"""Gen listesi veren kaynakların (Hpo, Clinvar) ortak tabanı: listeden crn nesnelerine.

Alt sınıf bir sorgudan gen satırları üretir (``genes``); bu taban satırları tek DNA,
tek transkript ve tek proteinli ``Gene`` nesnelerine çevirir, dizilere ayırır ve
gruplar hâlinde tek CRN listesine yazar. Alt sınıfın sağlaması gerekenler:

==================  ============================================================
 ``genes``          filtrelere uyan gen satırları; her satırda ``ncbi_gene_id``
                    ve ``gene_symbol`` bulunur
 ``_describe``      oluşturulan ``Gene``'in features'ına kaynağa özgü bilgiyi yazar
 ``release``        kullanılan veri sürümü
 ``FEATURE_KEYS``   ``to_sequences``'ın genden diziye taşıdığı anahtarlar
 ``RELEASE_KEY``    sürümün yazıldığı features anahtarı
 ``GROUP_FILTER``   ``pack``'te grup değeri sözlük değilse aktarıldığı filtre
==================  ============================================================
"""
from __future__ import annotations

import csv
import re
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence as Seq, Tuple

from crn import fileio
from crn.fetch._store import LocalSource
from crn.fetch.gencode import AmbiguousGeneError, Gencode
from crn.fetch.hgnc import Hgnc
from crn.gene import Gene
from crn.seq import Sequence

SEQUENCE_KINDS = ("dna", "rna", "protein")
OVERLAP_POLICIES = ("exclude", "first", "error")


class GeneSetSource(LocalSource):
    """Gen listesini ``Gene`` nesnelerine, dizilere ve CRN gruplarına çeviren taban."""

    FEATURE_KEYS: Tuple[str, ...] = ()
    RELEASE_KEY = ""
    GROUP_FILTER = ""

    def __init__(self, sys_dir: Optional[str | Path] = None, verbose: bool = True):
        super().__init__(sys_dir, verbose)
        self.__gencode: Optional[Gencode] = None

    # --- alt sınıfın sağladıkları ----------------------------------------------

    @property
    def release(self) -> Optional[str]:
        raise NotImplementedError

    def genes(self, *args: Any, **filters: Any) -> List[Dict[str, Any]]:
        raise NotImplementedError

    def _describe(self, gene: Gene, row: Dict[str, Any], filters: Dict[str, Any]) -> None:
        raise NotImplementedError

    # --- ortak yardımcılar -----------------------------------------------------

    @staticmethod
    def _check_choices(name: str, values: List[Any], allowed: Seq[Any]) -> None:
        unknown = [value for value in values if value not in allowed]
        if unknown:
            raise ValueError(f"Bilinmeyen {name}: {', '.join(map(str, unknown))}; "
                             f"seçenekler: {', '.join(map(str, allowed))}")

    @staticmethod
    def _row_key(row: Dict[str, Any]) -> str:
        """Gen satırının kimliği: NCBI Gene ID, yoksa sembol."""
        return str(row["ncbi_gene_id"] or row["gene_symbol"])

    def _gencode(self) -> Gencode:
        if self.__gencode is None:
            self.__gencode = Gencode(sys_dir=self.sys_dir, verbose=False)
        return self.__gencode

    def _canonical(self, source: Gencode, gene_id: str, symbol: str) -> Gene:
        """Sembolle ana transkript; aynı adlı GENCODE kayıtları ayrılarak.

        Ad birden çok kayda uyarsa tek protein kodlayan kayıt alınır (CAST). İkisi
        de kodluyorsa (SHOX, PDE8B) NCBI kimliği HGNC'de Ensembl kimliğine
        çevrilir; kimlik yoksa ya da HGNC yerelde yoksa hata aynen geçer.
        """
        try:
            return source.canonical(symbol)
        except AmbiguousGeneError:
            coding = source.genes(gene_type="protein_coding", name=f"^{re.escape(symbol)}$")
            if len(coding) == 1:
                return source.canonical(coding[0]["gene_id"])
            hgnc = Hgnc(sys_dir=self.sys_dir, verbose=False)
            if not gene_id or not hgnc.is_ready():
                raise
            rows = hgnc.query(entrez_id=gene_id, fields=["ensembl_gene_id"])
            if not rows or not rows[0]["ensembl_gene_id"]:
                raise
            return source.canonical(rows[0]["ensembl_gene_id"])

    def _isoform(self, source: Any, gene_id: str, symbol: str) -> Gene:
        """Genin tek DNA, tek transkript ve tek proteinli hali (CUB'un istediği biçim)."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")     # izoform ayıklama bildirimleri
            if isinstance(source, Gencode):
                if not symbol:
                    raise ValueError("gen sembolü yok")
                canonical = self._canonical(source, gene_id, symbol)
                isoforms = canonical.isoforms()
                if not isoforms:        # PANK2 gibi ana transkriptin ORF'u eşleşmeyen genler
                    isoforms = source.gene(canonical.features.first("gene_id")).isoforms()
                    for gene in isoforms:
                        gene.features.set("canonical_tag", "other_transcript")
            else:
                isoforms = source.gene(gene_id or symbol).isoforms()
                preferred = [g for g in isoforms
                             if str(g.features.first("transcript_id")).startswith("NM_")]
                isoforms = preferred or isoforms
        if not isoforms:
            raise ValueError("protein ve transkript eşleşmedi")
        return isoforms[0]

    def _build(self, wanted: List[Dict[str, Any]], source: Any, group: Optional[str],
               key: str, filters: Dict[str, Any], desc: str
               ) -> Tuple[List[Gene], Dict[str, str]]:
        """Gen satırlarından ``Gene`` listesi; kurulamayanlar nedenleriyle ayrı döner."""
        source = source if source is not None else self._gencode()
        genes, skipped = [], {}
        for row in self._progress(wanted, desc, "gen"):
            label = row["gene_symbol"] or str(row["ncbi_gene_id"])
            try:
                gene = self._isoform(source, str(row["ncbi_gene_id"] or ""), row["gene_symbol"])
            except (KeyError, ValueError, RuntimeError) as error:
                skipped[label] = str(error).strip("'\"")
                continue
            self._describe(gene, row, filters)
            if group is not None:
                gene.features.set(key, group)
            genes.append(gene)
        return genes, skipped

    def _to_genes(self, filters: Dict[str, Any], source: Any, group: Optional[str],
                  key: str) -> List[Gene]:
        genes, skipped = self._build(self.genes(**filters), source, group, key, filters,
                                     "Genler çıkarılıyor")
        if skipped:
            warnings.warn(f"{len(skipped)} gen atlandı: {', '.join(sorted(skipped)[:10])}"
                          f"{' ...' if len(skipped) > 10 else ''}", UserWarning, stacklevel=3)
        return genes

    def _to_sequences(self, kind: str, filters: Dict[str, Any], source: Any) -> List[Sequence]:
        self._check_choices("kind", [kind], SEQUENCE_KINDS)
        sequences = []
        for gene in self._to_genes(filters, source, None, "group"):
            if kind == "dna":
                item = gene.dna
            elif kind == "rna":
                item = next(iter(gene.transcripts.values()))
            else:
                item = next(iter(gene.proteins.values()))
            item.features.set("gene_name", gene.features.first("gene_name"))
            for name in self.FEATURE_KEYS:
                item.features.remove(name)
                for value in gene.features.get(name, []):
                    item.features.add(name, value)
            item.features.set(self.RELEASE_KEY, self.release)
            sequences.append(item)
        return sequences

    def pack(self, groups: Dict[str, Any], file_name: str | Path, source: Any = None,
             key: str = "group", overlap: str = "exclude", **filters: Any) -> Dict[str, Any]:
        """Gen gruplarını tek CRN listesine yazar.

        Dosya ``DifferentialCodonUsage.from_crn`` ile doğrudan açılır; genin
        grubu features'ta ``key`` anahtarındadır.

        Args:
            groups: grup adı -> sorgu değeri (``GROUP_FILTER`` filtresine gider) ya
                da gruba özel filtre sözlüğü.
            file_name: Yazılacak CRN dosyası.
            source: ``to_genes`` ile aynı.
            key: Grup anahtarı.
            overlap: Birden çok gruba giren gen için ``"exclude"`` (hepsinden
                çıkar, varsayılan), ``"first"`` (ilk grupta kalır) ya da ``"error"``.
            **filters: Tüm gruplara uygulanan filtreler.

        Returns:
            ``file, release, groups`` (grup -> gen sayısı), ``overlap`` (çakışan
            semboller) ve ``skipped`` (atlanan gen -> neden).

        Raises:
            ValueError: overlap geçersizse, gruplar çakışıyorsa (``"error"``) ya
                da bir grup boş kalırsa.
        """
        self._check_choices("overlap", [overlap], OVERLAP_POLICIES)
        if not groups:
            raise ValueError("En az bir grup gerekli")
        wanted: Dict[str, List[Dict[str, Any]]] = {}
        settings: Dict[str, Dict[str, Any]] = {}
        owners: Dict[str, List[str]] = defaultdict(list)
        labels: Dict[str, str] = {}
        for name, spec in groups.items():
            settings[name] = {**filters, **(spec if isinstance(spec, dict) else {self.GROUP_FILTER: spec})}
            wanted[name] = self.genes(**settings[name])
            for row in wanted[name]:
                owners[self._row_key(row)].append(name)
                labels[self._row_key(row)] = row["gene_symbol"] or str(row["ncbi_gene_id"])

        shared = {gene: names for gene, names in owners.items() if len(names) > 1}
        if shared and overlap == "error":
            raise ValueError(f"{len(shared)} gen birden çok grupta: "
                             f"{', '.join(sorted(labels[g] for g in shared)[:10])}")
        for name in wanted:
            wanted[name] = [row for row in wanted[name] if self._row_key(row) not in shared
                            or (overlap == "first" and shared[self._row_key(row)][0] == name)]
        if shared:
            warnings.warn(f"{len(shared)} gen birden çok grupta; overlap='{overlap}' uygulandı",
                          UserWarning, stacklevel=2)

        items: List[Gene] = []
        skipped: Dict[str, str] = {}
        counts: Dict[str, int] = {}
        for name in wanted:
            genes, missed = self._build(wanted[name], source, name, key, settings[name],
                                        f"{name} genleri")
            if not genes:
                raise ValueError(f"'{name}' grubunda gen kalmadı")
            items += genes
            skipped.update(missed)
            counts[name] = len(genes)
        fileio.write_crn_list(items, file_name)
        return {"file": str(file_name), "release": self.release, "groups": counts,
                "overlap": sorted(labels[g] for g in shared), "skipped": skipped}

    @staticmethod
    def write_tsv(rows: Seq[Dict[str, Any]], file_name: str | Path) -> None:
        """Sorgu sonucunu sekmeyle ayrılmış dosyaya yazar."""
        if not rows:
            raise ValueError("Yazılacak satır yok")
        with open(file_name, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
