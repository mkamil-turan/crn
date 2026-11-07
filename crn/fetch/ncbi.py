"""NCBI'dan DNA, RNA, protein ve gen çıkarma.

``Gencode`` ile aynı mantıkla çalışır: indirilen her kayıt yerel ``sys/``
klasöründe saklanır; bir dizi istendiğinde önce bu klasöre bakılır, kayıt
yoksa ``Ncbi.update()`` ile NCBI'dan indirilir. İstekler ``NcbiSession``
üzerinden (hız sınırı ve yeniden deneme ile) yapılır.

``sys/`` klasörünün yapısı::

    sys/ncbi/
        symbols/<organizma>.json        gen sembolü -> NCBI Gene ID
        genes/<gene_id>.json            genin özeti ve RNA/protein listesi
                                        (update tamamlandığında yazılır)
        nuccore/<erişim>.fa, .ft        transkript dizisi ve özellik tablosu
        protein/<erişim>.fa             protein dizisi
        region/<erişim>_<baş>_<son>_<zincir>.fa
                                        genomik bölge

Kullanım::

    from crn.fetch import Ncbi

    ncbi = Ncbi()                       # insan; NcbiSession kendiliğinden açılır
    ncbi.update("MEFV")                 # isteğe bağlı: indir (yoksa ilk istekte yapılır)
    mefv = ncbi.gene("MEFV")            # Gene: genomik DNA + tüm RefSeq RNA'lar + proteinler
    dna = ncbi.dna("MEFV")              # genin genomik DNA'sı (kodlayan zincir, 5'->3')
    rna = ncbi.rna("NM_000243")         # transkript (sürüm numarası olmadan da olur)
    prt = ncbi.protein("NM_000243")     # transkriptin proteini (NP_ erişimi de olur)
    bolge = ncbi.region("NC_000016.10", 3242028, 3256627, "-")
"""
from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import (Any, Dict, Iterator, List, Optional, Tuple)

from crn.constants import DNA_TO_RNA
from crn.seq import DNA, RNA, PRT, GenomicStrand
from crn.gene import Gene
from crn.fetch._store import LocalSource, atomic_write, split_fasta, read_single_fasta
from crn.fetch.session import NcbiSession


NCBI_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# organizma -> NCBI Taxonomy ID'leri
ORGANISMS = {
    "human":     ("9606",),
    "mouse":     ("10090",),
    "rat":       ("10116",),
    "zebrafish": ("7955",),
    "fruit_fly": ("7227",),
    "worm":      ("6239",),
}

# Tek istekte indirilecek en çok kayıt.
_BATCH = 100
_CDS = re.compile(r"^[<>]?(\d+)\t[<>]?(\d+)\tCDS\s*$")
_QUALIFIER = re.compile(r"^\t\t\t(\S+)\t(.*)$")


def _split_feature_table(text: str) -> Iterator[Tuple[str, str]]:
    """Çok kayıtlı özellik tablosunu (erişim, kaydın metni) çiftlerine böler."""
    for block in re.split(r"(?m)^(?=>Feature)", text):
        match = re.match(r">Feature \w+\|([^|]+)\|", block)
        if match:
            yield match.group(1), block.rstrip("\n") + "\n"


def _parse_feature_table(text: str) -> Dict[str, Any]:
    """Özellik tablosundan gen adını, GeneID'yi, ilk CDS'i ve proteinini çıkarır.

    Returns:
        ``{"gene_name", "gene_id", "cds": (başlangıç, bitiş) ya da None,
        "protein_id"}``. Koordinatlar transkript üzerinde 1 tabanlıdır;
        eksik uç işaretleri (``<``, ``>``) atılır.
    """
    result: Dict[str, Any] = {"gene_name": None, "gene_id": None, "cds": None, "protein_id": None}
    feature = None
    for line in text.splitlines():
        if not line.startswith("\t"):
            match = _CDS.match(line)
            if match and result["cds"] is None:
                feature = "CDS"
                result["cds"] = (int(match.group(1)), int(match.group(2)))
            elif "\t" in line and not line.startswith(">"):
                parts = line.split("\t")
                feature = parts[2] if len(parts) > 2 and parts[2] else feature
            continue
        match = _QUALIFIER.match(line)
        if not match:
            continue
        key, value = match.groups()
        if feature == "gene" and key == "gene" and result["gene_name"] is None:
            result["gene_name"] = value
        elif feature == "gene" and key == "db_xref" and value.startswith("GeneID:"):
            result["gene_id"] = result["gene_id"] or value.split(":", 1)[1]
        elif feature == "CDS" and key == "protein_id" and result["protein_id"] is None:
            result["protein_id"] = value.strip("|").split("|")[-1]
    return result


class Ncbi(LocalSource):
    """NCBI'dan DNA, RNA, PRT ve Gene nesneleri çıkarır.

    Args:
        organism: Gen sembolü aramasında kullanılacak organizma
            (``ORGANISMS`` anahtarlarından biri).
        session: HTTP oturumu; verilmezse API anahtarsız bir ``NcbiSession``
            açılır (saniyede 3 istek).
        sys_dir: Verilerin tutulacağı klasör; varsayılan proje kökündeki ``sys/``.
        verbose: True ise indirme ve okuma adımları ilerleme çubuklarıyla
            (tqdm) gösterilir.
    """

    SOURCE = "NCBI"

    def __init__(self, organism: str = "human", session: Optional[NcbiSession] = None,
                 sys_dir: Optional[str | Path] = None, verbose: bool = True):
        if organism not in ORGANISMS:
            raise ValueError(f"Desteklenmeyen organizma: {organism}; "
                             f"seçenekler: {', '.join(ORGANISMS)}")
        super().__init__(sys_dir, verbose)
        self.organism = organism
        self.session = session if session is not None else NcbiSession()

    # --- Klasörler ------------------------------------------------------------

    @property
    def ncbi_dir(self) -> Path:
        return self.sys_dir / "ncbi"

    def _gene_file(self, gene_id: str) -> Path:
        return self.ncbi_dir / "genes" / f"{gene_id}.json"

    def _symbols_file(self) -> Path:
        return self.ncbi_dir / "symbols" / f"{self.organism}.json"

    def _record_file(self, db: str, accession: str, suffix: str) -> Path:
        return self.ncbi_dir / db / f"{accession}.{suffix}"

    def _cached_accession(self, db: str, accession: str) -> Optional[str]:
        """sys/ içindeki kaydın tam erişimi; sürümsüz verilirse en yeni sürüm."""
        if "." in accession:
            return accession if self._record_file(db, accession, "fa").is_file() else None
        folder = self.ncbi_dir / db
        versions = [p.stem for p in folder.glob(f"{accession}.*.fa")] if folder.is_dir() else []
        return max(versions, key=lambda a: int(a.rsplit(".", 1)[1])) if versions else None

    # --- NCBI istekleri -------------------------------------------------------

    def _get(self, endpoint: str, **params: Any):
        return self.session.get(f"{NCBI_URL}/{endpoint}", params=params)

    def _efetch(self, db: str, ids: List[str], rettype: str, desc: str, **params: Any) -> str:
        """Kayıtları _BATCH'lik gruplar halinde indirip metinlerini birleştirir."""
        texts = []
        with self._bar(len(ids), desc, "kayıt") as bar:
            for i in range(0, len(ids), _BATCH):
                chunk = ids[i:i + _BATCH]
                response = self._get("efetch.fcgi", db=db, id=",".join(chunk),
                                     rettype=rettype, retmode="text", **params)
                texts.append(response.text)
                bar.update(len(chunk))
        return "\n".join(texts)

    def search_gene(self, gene: str) -> str:
        """Gen sembolünü (ör. ``"MEFV"``) NCBI Gene ID'sine çevirir.

        Rakamlardan oluşan değer zaten Gene ID kabul edilir. NCBI araması takma
        adları da eşleştirdiği için birden çok sonuç varsa resmi sembolü tam
        tutan gen seçilir (``MIP`` için TNPO1 değil MIP); yoksa ilk sonuç alınır.
        Bulunan eşleşmeler sys/ içinde saklanır; aynı sembol tekrar sorulduğunda
        ağa çıkılmaz.

        Raises:
            KeyError: Sembol bu organizmada bulunamazsa.
        """
        gene = str(gene)
        if gene.isdigit():
            return gene
        symbols_file = self._symbols_file()
        symbols = json.loads(symbols_file.read_text()) if symbols_file.is_file() else {}
        if gene in symbols:
            return symbols[gene]

        taxa = " OR ".join(f"{taxon}[Taxonomy ID]" for taxon in ORGANISMS[self.organism])
        term = f"{gene}[Gene Name] AND ({taxa}) AND alive[prop]"
        data = self._get("esearch.fcgi", db="gene", term=term, retmode="json").json()
        ids = data.get("esearchresult", {}).get("idlist", [])
        if not ids:
            raise KeyError(f"NCBI: '{gene}' geni {self.organism} organizmasında bulunamadı")
        found = ids[0]
        if len(ids) > 1:
            summaries = self._get("esummary.fcgi", db="gene", id=",".join(ids),
                                  retmode="json").json().get("result", {})
            found = next((i for i in ids if summaries.get(i, {}).get("name", "").casefold()
                          == gene.casefold()), found)
        symbols[gene] = found
        atomic_write(symbols_file, json.dumps(symbols, indent=1, sort_keys=True))
        return found

    def _summary(self, gene_id: str) -> Dict[str, Any]:
        data = self._get("esummary.fcgi", db="gene", id=gene_id, retmode="json").json()
        summary = data.get("result", {}).get(gene_id)
        if not summary or "error" in summary:
            raise KeyError(f"NCBI: Gene ID {gene_id} bulunamadı")
        return summary

    def _elink(self, gene_id: str, db: str, linkname: str) -> List[str]:
        data = self._get("elink.fcgi", dbfrom="gene", db=db, id=gene_id,
                         linkname=linkname, retmode="json").json()
        ids = []
        for linkset in data.get("linksets", []):
            for linksetdb in linkset.get("linksetdbs", []):
                if linksetdb.get("linkname") == linkname:
                    ids.extend(str(i) for i in linksetdb.get("links", []))
        return ids

    # --- update ---------------------------------------------------------------

    def update(self, gene: str, force: bool = False) -> str:
        """Genin tüm kayıtlarını NCBI'dan indirip sys/ içine yazar.

        İndirilenler: gen özeti, genin genomik bölgesi, gene bağlı tüm RefSeq
        RNA'ları (NM_, NR_, XM_, XR_) ile özellik tabloları ve bu RNA'ların
        kodladığı proteinler. sys/ içinde zaten bulunan RNA, protein ve bölge
        kayıtları yeniden indirilmez; ``force=True`` ise hepsi yenilenir.
        En sonda ``genes/<gene_id>.json`` yazılır; bu dosya yoksa gen hazır
        sayılmaz.

        Args:
            gene: Gen sembolü ya da NCBI Gene ID.
            force: True ise tüm kayıtlar yeniden indirilir.

        Returns:
            NCBI Gene ID.

        Raises:
            KeyError: Gen bulunamazsa.
        """
        with self._bar(4, f"NCBI {gene}", "adım") as bar:
            bar.set_postfix_str("arama ve özet")
            gene_id = self.search_gene(gene)
            summary = self._summary(gene_id)
            bar.update(1)

            bar.set_postfix_str("RNA'lar")
            rna_ids = self._elink(gene_id, "nuccore", "gene_nuccore_refseqrna")
            transcripts = self._download_rnas(rna_ids, force) if rna_ids else []
            bar.update(1)

            bar.set_postfix_str("proteinler")
            proteins = {}
            for accession in transcripts:
                ft = self._record_file("nuccore", accession, "ft")
                protein_id = _parse_feature_table(ft.read_text())["protein_id"] if ft.is_file() else None
                if protein_id:
                    proteins[protein_id] = accession
            self._download_proteins(list(proteins), force)
            bar.update(1)

            bar.set_postfix_str("genomik bölge")
            region = self._gene_region(summary)
            if region is not None:
                self._download_region(*region, force=force)
            bar.update(1)

        organism = summary.get("organism") or {}
        record = {
            "gene_id": gene_id,
            "gene_name": summary.get("name", ""),
            "description": summary.get("description", ""),
            "organism": organism.get("scientificname", ""),
            "chromosome": summary.get("chromosome", ""),
            "map_location": summary.get("maplocation", ""),
            "aliases": [a.strip() for a in summary.get("otheraliases", "").split(",") if a.strip()],
            "region": list(region) if region else None,
            "transcripts": transcripts,
            "proteins": proteins,
            "updated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        atomic_write(self._gene_file(gene_id), json.dumps(record, indent=1))
        return gene_id

    def _download_rnas(self, ids: List[str], force: bool) -> List[str]:
        """RNA'ların FASTA ve özellik tablolarını indirir; erişimlerini döndürür.

        NCBI UID ile çalıştığı için hangi kaydın sys/ içinde olduğu indirmeden
        bilinemez; bu yüzden FASTA her zaman indirilir, özellik tabloları
        yalnızca eksik olanlar için indirilir.
        """
        accessions = []
        for accession, text in split_fasta(self._efetch("nuccore", ids, "fasta", "RNA dizileri")):
            if force or not self._record_file("nuccore", accession, "fa").is_file():
                atomic_write(self._record_file("nuccore", accession, "fa"), text)
            accessions.append(accession)
        missing = [a for a in accessions
                   if force or not self._record_file("nuccore", a, "ft").is_file()]
        if missing:
            for accession, text in _split_feature_table(
                    self._efetch("nuccore", missing, "ft", "RNA özellik tabloları")):
                atomic_write(self._record_file("nuccore", accession, "ft"), text)
        # Önce doğrulanmış (NM_, NR_), sonra tahmini (XM_, XR_) kayıtlar
        return sorted(accessions, key=lambda a: (a.startswith("X"), a))

    def _download_proteins(self, accessions: List[str], force: bool) -> None:
        missing = [a for a in accessions
                   if force or not self._record_file("protein", a, "fa").is_file()]
        if missing:
            for accession, text in split_fasta(
                    self._efetch("protein", missing, "fasta", "Protein dizileri")):
                atomic_write(self._record_file("protein", accession, "fa"), text)

    def _region_file(self, accession: str, start: int, end: int, strand: str) -> Path:
        sign = "plus" if strand == "+" else "minus"
        return self.ncbi_dir / "region" / f"{accession}_{start}_{end}_{sign}.fa"

    def _download_region(self, accession: str, start: int, end: int, strand: str,
                         force: bool = False) -> Path:
        path = self._region_file(accession, start, end, strand)
        if force or not path.is_file():
            text = self._efetch("nuccore", [accession], "fasta", "Genomik bölge",
                                seq_start=start, seq_stop=end,
                                strand=1 if strand == "+" else 2)
            records = list(split_fasta(text))
            if not records:
                raise KeyError(f"NCBI: {accession}:{start}-{end} bölgesi indirilemedi")
            atomic_write(path, records[0][1])
        return path

    @staticmethod
    def _gene_region(summary: Dict[str, Any]) -> Optional[Tuple[str, int, int, str]]:
        """Özetteki genomik koordinatlardan (erişim, başlangıç, bitiş, zincir).

        NCBI 0 tabanlı koordinat verir; 1 tabanlıya çevrilir. ``chrstart >
        chrstop`` ise gen (-) zincirdedir.
        """
        info = (summary.get("genomicinfo") or [{}])[0]
        accession, chrstart, chrstop = info.get("chraccver"), info.get("chrstart"), info.get("chrstop")
        if accession is None or chrstart is None or chrstop is None:
            return None
        strand = "-" if chrstart > chrstop else "+"
        return accession, min(chrstart, chrstop) + 1, max(chrstart, chrstop) + 1, strand

    # --- sys/ içinden okuma ---------------------------------------------------

    def _gene_record(self, gene: str) -> Dict[str, Any]:
        """Genin sys/ içindeki kaydını döndürür; yoksa ya da eksikse update() çağırır."""
        gene_id = self.search_gene(gene)
        path = self._gene_file(gene_id)
        if path.is_file():
            record = json.loads(path.read_text())
            complete = (all(self._record_file("nuccore", a, "fa").is_file() for a in record["transcripts"])
                        and all(self._record_file("protein", a, "fa").is_file() for a in record["proteins"])
                        and (record["region"] is None
                             or self._region_file(*record["region"]).is_file()))
            if complete:
                return record
        self.update(gene_id)
        return json.loads(path.read_text())

    def _features(self, obj, **values) -> None:
        super()._features(obj, empty=(None, "", []), **values)

    # --- Çıkarma -------------------------------------------------------------

    def region(self, accession: str, start: int, end: int, strand: str = "+",
               name: str = "") -> DNA:
        """Bir genomik bölgeyi DNA olarak döndürür (sys/ içinde yoksa indirir).

        Args:
            accession: Kromozom ya da contig erişimi (ör. ``"NC_000016.10"``).
            start, end: 1 tabanlı, iki ucu dahil konumlar.
            strand: ``"+"`` referans zincir; ``"-"`` ise ters tümleyici alınır.
                Her iki durumda DNA 5'->3' yazılır.
            name: DNA'nın adı; boşsa ``NC_000016.10:3242028-3256627(-)`` olur.

        Raises:
            ValueError: Konumlar ya da strand geçersizse.
        """
        if strand not in ("+", "-"):
            raise ValueError(f"strand '+' ya da '-' olmalı, {strand!r} verildi")
        if not 1 <= start <= end:
            raise ValueError(f"Geçersiz bölge: {start}-{end}")
        _, seq = read_single_fasta(self._download_region(accession, start, end, strand))
        dna = DNA(seq, name=name or f"{accession}:{start}-{end}({strand})")
        self._features(dna, chromosome=accession, start=start, end=end, chromosome_strand=strand)
        return dna

    def dna(self, gene: str) -> DNA:
        """Genin genomik DNA'sını kodlayan (sense) zincirde, 5'->3' döndürür.

        Args:
            gene: Gen sembolü (ör. ``"MEFV"``) ya da NCBI Gene ID.

        Raises:
            KeyError: Gen ya da genin genomik koordinatı bulunamazsa.
        """
        record = self._gene_record(gene)
        if record["region"] is None:
            raise KeyError(f"NCBI: Gene ID {record['gene_id']} için genomik koordinat yok")
        dna = self.region(*record["region"], name=record["gene_name"] or record["gene_id"])
        self._features(dna, gene_id=record["gene_id"], gene_name=record["gene_name"],
                       description=record["description"])
        return dna

    def rna(self, accession: str,
            genomic_strand: GenomicStrand = GenomicStrand.Unknown) -> RNA:
        """RefSeq transkriptini RNA olarak döndürür (sys/ içinde yoksa indirir).

        Protein kodlayan transkriptlerde features'a ``utr5``, ``cds`` ve
        ``utr3`` bölgeleri ``(başlangıç, bitiş)`` olarak yazılır (transkript
        üzerinde 1 tabanlı, iki ucu dahil). CDS şöyle alınır::

            start, end = rna.features.first("cds")
            cds = rna[start - 1:end]

        Args:
            accession: RefSeq erişimi (ör. ``"NM_000243.3"``; sürümsüz de olur).
            genomic_strand: RNA'nın transkript edildiği genomik zincir; ``gene()``
                bunu genin zincirinden doldurur.

        Raises:
            KeyError: Kayıt NCBI'da bulunamazsa.
        """
        found = self._cached_accession("nuccore", accession)
        if found is None or not self._record_file("nuccore", found, "ft").is_file():
            downloaded = self._download_rnas([accession], force=False)
            if not downloaded:
                raise KeyError(f"NCBI: '{accession}' transkripti bulunamadı")
            found = downloaded[0]

        title, seq = read_single_fasta(self._record_file("nuccore", found, "fa"))
        ft = self._record_file("nuccore", found, "ft")
        table = _parse_feature_table(ft.read_text()) if ft.is_file() else {}
        rna = RNA(seq.translate(DNA_TO_RNA), name=found,
                  genomic_strand=genomic_strand)
        self._features(rna, transcript_id=found, title=title.split(None, 1)[-1],
                       predicted=found.startswith("X"), gene_id=table.get("gene_id"),
                       gene_name=table.get("gene_name"), protein_id=table.get("protein_id"))
        cds = table.get("cds")
        if cds:
            start, end = cds
            if start > 1:
                rna.features.set("utr5", (1, start - 1))
            rna.features.set("cds", (start, end))
            if end < len(rna):
                rna.features.set("utr3", (end + 1, len(rna)))
        return rna

    def protein(self, accession: str) -> PRT:
        """Proteini döndürür (sys/ içinde yoksa indirir).

        Args:
            accession: Protein erişimi (NP_, XP_) ya da proteini kodlayan
                transkriptin erişimi (NM_, XM_); sürümsüz de olur.

        Raises:
            KeyError: Protein bulunamazsa (ör. kodlamayan RNA).
        """
        transcript_id = None
        if accession.startswith(("NM_", "XM_", "NR_", "XR_")):
            rna = self.rna(accession)
            transcript_id = rna.name
            accession = rna.features.first("protein_id")
            if accession is None:
                raise KeyError(f"NCBI: {transcript_id} protein kodlamıyor")

        found = self._cached_accession("protein", accession)
        if found is None:
            self._download_proteins([accession], force=False)
            found = self._cached_accession("protein", accession)
            if found is None:
                raise KeyError(f"NCBI: '{accession}' proteini bulunamadı")

        title, seq = read_single_fasta(self._record_file("protein", found, "fa"))
        prt = PRT(seq, name=found)
        self._features(prt, protein_id=found, transcript_id=transcript_id,
                       title=title.split(None, 1)[-1])
        return prt

    def gene(self, gene: str) -> Gene:
        """Geni tüm içeriğiyle döndürür.

        Gene nesnesi şunları içerir: genomik DNA (``dna()``; koordinat yoksa
        boş), gene bağlı tüm RefSeq RNA'ları (kimlik = erişim) ve bunların
        kodladığı proteinler (kimlik = erişim). Proteinin hangi transkriptten
        geldiği proteinin ``transcript_id`` özelliğindedir.

        Args:
            gene: Gen sembolü (ör. ``"MEFV"``) ya da NCBI Gene ID.
        """
        record = self._gene_record(gene)
        name = record["gene_name"] or record["gene_id"]
        result = Gene(name, self.dna(record["gene_id"]) if record["region"] else None)

        region = record["region"]
        self._features(result, gene_id=record["gene_id"], gene_name=record["gene_name"],
                       description=record["description"], organism=record["organism"],
                       chromosome=record["chromosome"], map_location=record["map_location"])
        for alias in record["aliases"]:
            result.features.add("aliases", alias)
        if region:
            result.features.set("genomic_region", tuple(region))
        genomic_strand = (GenomicStrand.Minus if region and region[3] == "-"
                          else GenomicStrand.Plus if region else GenomicStrand.Unknown)

        protein_of = {transcript: protein for protein, transcript in record["proteins"].items()}
        for accession in self._progress(record["transcripts"], f"{name} transkriptleri", "transkript"):
            result.add_transcript(self.rna(accession, genomic_strand), accession)
            if accession in protein_of:
                prt = self.protein(protein_of[accession])
                prt.features.set("transcript_id", accession)
                result.add_protein(prt, prt.name)
        return result
