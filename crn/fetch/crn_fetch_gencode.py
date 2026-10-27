"""GENCODE verilerinden DNA, RNA, protein ve gen çıkarma.

Veriler yerel ``sys/`` klasöründe tutulur. Bir dizi istendiğinde önce bu
klasöre bakılır; veri yoksa ``Gencode.update()`` ile indirilir ve hazırlanır.

``sys/`` klasörünün yapısı::

    sys/gencode/<tür>/<sürüm>/          ör. sys/gencode/human/v50/
        raw/                            GENCODE'dan indirilen .gz dosyaları
        genome.fa, genome.fai           açılmış genom ve konum dizini
        transcripts.fa, transcripts.idx transkript dizileri ve dizini
        translations.fa, translations.idx
                                        protein dizileri ve dizini
        genes.tsv, transcripts.tsv      GTF'den çıkarılan gen ve transkript
                                        koordinatları
        regions.tsv                     protein kodlayan transkriptlerin
                                        UTR5/CDS/UTR3 koordinatları
        manifest.json                   hazırlık tamamlandığında yazılır

Kullanım::

    from crn.fetch.crn_fetch_gencode import Gencode

    gc = Gencode()                  # insan, sys/ içindeki en yeni sürüm
    gc.update()                     # isteğe bağlı: indir/hazırla (yoksa ilk istekte yapılır)
    mefv = gc.gene("MEFV")          # Gene: genomik DNA + tüm RNA'lar + proteinler
    dna = gc.dna("MEFV")            # genin genomik DNA'sı (kodlayan zincir, 5'->3')
    rna = gc.rna("ENST00000219596") # transkript (sürüm numarası olmadan da olur)
    prt = gc.protein("ENST00000219596")   # transkriptin proteini (ENSP de olur)
    bolge = gc.region("chr16", 3242028, 3256627, "-")
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import (Dict, Iterator, List, Optional, Tuple)

from tqdm import tqdm

from crn.io import (DNA, RNA, PRT, Gene, GenomicStrand, Maps)


GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode"

# tür -> (GENCODE klasörü, sürüm öneki)
SPECIES = {
    "human": ("Gencode_human", "v"),
    "mouse": ("Gencode_mouse", "vM"),
}

# Proje kökündeki sys/ klasörü.
DEFAULT_SYS_DIR = Path(__file__).resolve().parents[2] / "sys"

_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')
_VERSION = re.compile(r"\.\d+")
_GENES_COLUMNS = ["gene_id", "gene_name", "gene_type", "chromosome", "start", "end", "strand"]
_TRANSCRIPTS_COLUMNS = ["transcript_id", "gene_id", "transcript_name", "transcript_type",
                        "chromosome", "start", "end", "strand"]
_REGION = re.compile(r"(UTR5|CDS|UTR3):(\d+)-(\d+)")
# Hazırlık sonunda bulunması gereken dosyalar; biri eksikse sürüm hazır sayılmaz.
_OUTPUTS = ["genome.fa", "genome.fai", "transcripts.fa", "transcripts.idx",
            "translations.fa", "translations.idx", "genes.tsv", "transcripts.tsv",
            "regions.tsv", "manifest.json"]


def _unversioned(identifier: str) -> str:
    """ENSG00000103313.13 -> ENSG00000103313 (``_PAR_Y`` gibi ekler korunur)."""
    return _VERSION.sub("", identifier)


def _reverse_complement(seq: str) -> str:
    return seq.translate(Maps.DEFAULT_COMPLEMENT_MAP_FOR_DNA.value)[::-1]


class Gencode:
    """GENCODE sürümünden DNA, RNA, PRT ve Gene nesneleri çıkarır.

    Args:
        species: ``"human"`` ya da ``"mouse"``.
        release: ``"v50"``, ``"50"`` ya da fare için ``"vM39"`` gibi bir sürüm.
            None ise ``sys/`` içindeki en yeni hazır sürüm kullanılır; hiç
            yoksa ``update()`` GENCODE'daki en son sürümü indirir.
        sys_dir: Verilerin tutulacağı klasör; varsayılan proje kökündeki ``sys/``.
        base_url: GENCODE sunucusunun adresi.
        verbose: True ise indirme, hazırlık ve arama adımları ilerleme
            çubuklarıyla (tqdm) gösterilir.
    """

    def __init__(self, species: str = "human", release: Optional[str] = None,
                 sys_dir: Optional[str | Path] = None,
                 base_url: str = GENCODE_URL, verbose: bool = True):
        if species not in SPECIES:
            raise ValueError(f"Desteklenmeyen tür: {species}; seçenekler: {', '.join(SPECIES)}")
        self.species = species
        self.sys_dir = Path(sys_dir) if sys_dir is not None else DEFAULT_SYS_DIR
        self.base_url = base_url.rstrip("/")
        self.verbose = verbose
        self.__requested = self._normalize_release(release) if release else None
        self.__release = self.__requested
        self.__tables: Optional[dict] = None

    # --- Sürüm ve klasörler -------------------------------------------------

    def _normalize_release(self, release: str) -> str:
        prefix = SPECIES[self.species][1]
        release = str(release)
        if release.startswith(prefix):
            return release
        return prefix + release.lstrip("vM")

    @property
    def species_dir(self) -> Path:
        return self.sys_dir / "gencode" / self.species

    @property
    def release(self) -> Optional[str]:
        """Kullanılan sürüm. Belirtilmediyse sys/ içindeki en yeni hazır sürüm."""
        if self.__release is None:
            installed = self.installed_releases()
            if installed:
                self.__release = installed[-1]
        return self.__release

    @property
    def release_dir(self) -> Path:
        if self.release is None:
            raise ValueError("Sürüm belirli değil; önce update() çalıştırın")
        return self.species_dir / self.release

    def installed_releases(self) -> List[str]:
        """sys/ içinde hazırlığı tamamlanmış sürümler, eskiden yeniye."""
        if not self.species_dir.is_dir():
            return []
        ready = [p.name for p in self.species_dir.iterdir() if (p / "manifest.json").is_file()]
        return sorted(ready, key=lambda r: int(re.sub(r"\D", "", r) or 0))

    def is_ready(self) -> bool:
        """Kullanılan sürüm sys/ içinde hazır mı (tüm hazırlık dosyaları var mı)?"""
        return (self.release is not None
                and all((self.release_dir / name).is_file() for name in _OUTPUTS))

    def latest_release(self) -> str:
        """GENCODE sunucusundaki en son sürümü döndürür (ağ bağlantısı gerekir)."""
        folder, prefix = SPECIES[self.species]
        listing = self._read_url(f"{self.base_url}/{folder}/latest_release/").decode("utf-8", "replace")
        match = re.search(rf"gencode\.({prefix}\d+)\.annotation\.gtf\.gz", listing)
        if not match:
            raise RuntimeError(f"GENCODE {self.species} için son sürüm bulunamadı")
        return match.group(1)

    def _remote_dir(self) -> str:
        folder, prefix = SPECIES[self.species]
        number = self.release[len(prefix):]
        return f"{self.base_url}/{folder}/release_{number if prefix == 'v' else 'M' + number}"

    # --- update ----------------------------------------------------------------

    def update(self, force: bool = False) -> str:
        """Seçili sürümün dosyalarını GENCODE'dan sys/ içine indirir ve hazırlar.

        Nesne bir sürümle oluşturulduysa o sürüm, oluşturulmadıysa sunucudaki
        en son sürüm hazırlanır (sys/ içinde daha eski bir sürüm olsa bile).
        İndirilen her dosya MD5SUMS ile doğrulanır. Daha önce indirilmiş ve
        doğrulanmış dosyalar yeniden indirilmez; daha önce hazırlanmış
        dosyalar yeniden hazırlanmaz, yalnızca eksikler tamamlanır. Hazırlık
        bitince ``manifest.json`` yazılır.

        Args:
            force: True ise dosyalar yeniden indirilir ve yeniden hazırlanır.

        Returns:
            Hazırlanan sürüm (ör. ``"v50"``).
        """
        return self._prepare(self.__requested or self.latest_release(), force)

    def _prepare(self, release: str, force: bool = False) -> str:
        """Verilen sürümü indirir ve eksik hazırlık dosyalarını üretir."""
        self.__release = release
        self.__tables = None
        if self.is_ready() and not force:
            self._log(f"GENCODE {self.species} {self.release} zaten hazır: {self.release_dir}")
            return self.release

        raw = self.release_dir / "raw"
        raw.mkdir(parents=True, exist_ok=True)
        (self.release_dir / "manifest.json").unlink(missing_ok=True)

        remote = self._remote_dir()
        checksums = self._read_checksums(remote)
        files = self._select_files(checksums)
        for name in files.values():
            self._download(f"{remote}/{name}", raw / name, checksums[name], force)

        # Her adım çıktısını önce .tmp olarak yazar ve sonda yerine koyar; bu
        # yüzden var olan bir çıktı tamamlanmış demektir ve atlanabilir.
        steps = [
            ("genome.fai", lambda: self._build_genome(raw / files["genome"])),
            ("transcripts.idx", lambda: self._build_fasta(
                raw / files["transcripts"], self.release_dir / "transcripts",
                "Transkriptler dizinleniyor")),
            ("translations.idx", lambda: self._build_fasta(
                raw / files["translations"], self.release_dir / "translations",
                "Proteinler dizinleniyor")),
            ("transcripts.tsv", lambda: self._build_tables(raw / files["annotation"])),
            ("regions.tsv", lambda: self._build_regions(raw / files["pc_transcripts"])),
        ]
        for output, build in steps:
            if force or not (self.release_dir / output).is_file():
                build()

        manifest = {
            "species": self.species,
            "release": self.release,
            "source": remote,
            "files": files,
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        (self.release_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
        self.__tables = None
        self._log(f"GENCODE {self.species} {self.release} hazır: {self.release_dir}")
        return self.release

    def _ensure(self) -> None:
        """Veri sys/ içinde hazır değilse update() ile indirir.

        Sürüm verilmediyse ve sys/ içinde hazır bir sürüm varsa ağa çıkılmaz.
        """
        if not self.is_ready():
            self._prepare(self.release or self.latest_release())

    def _select_files(self, checksums: Dict[str, str]) -> Dict[str, str]:
        release = self.release
        wanted = {
            "annotation": f"gencode.{release}.primary_assembly.annotation.gtf.gz",
            "transcripts": f"gencode.{release}.transcripts.fa.gz",
            "translations": f"gencode.{release}.pc_translations.fa.gz",
            "pc_transcripts": f"gencode.{release}.pc_transcripts.fa.gz",
        }
        genome = [n for n in checksums if re.fullmatch(r"GRC\w+\.primary_assembly\.genome\.fa\.gz", n)]
        if not genome:
            raise RuntimeError(f"{release}: MD5SUMS içinde primary_assembly genomu yok")
        wanted["genome"] = genome[0]
        missing = [n for n in wanted.values() if n not in checksums]
        if missing:
            raise RuntimeError(f"{release}: sunucuda bulunamayan dosyalar: {', '.join(missing)}")
        return wanted

    def _read_checksums(self, remote: str) -> Dict[str, str]:
        text = self._read_url(f"{remote}/MD5SUMS").decode("utf-8")
        checksums = {}
        for line in text.splitlines():
            parts = line.split()
            if len(parts) == 2:
                checksums[parts[1].lstrip("*")] = parts[0]
        return checksums

    def _read_url(self, url: str) -> bytes:
        with urllib.request.urlopen(url, timeout=60) as response:
            return response.read()

    def _download(self, url: str, target: Path, md5: str, force: bool) -> None:
        if target.exists() and not force:
            if self._md5(target, f"MD5 {target.name}") == md5:
                self._log(f"Var, atlandı: {target.name}")
                return
            self._log(f"MD5 tutmuyor, yeniden indiriliyor: {target.name}")

        part = target.with_name(target.name + ".part")
        digest = hashlib.md5()
        with urllib.request.urlopen(url, timeout=60) as response, open(part, "wb") as out:
            total = int(response.headers.get("Content-Length") or 0)
            with self._bar(total, f"İndiriliyor {target.name}") as bar:
                while chunk := response.read(1 << 20):
                    out.write(chunk)
                    digest.update(chunk)
                    bar.update(len(chunk))
        if digest.hexdigest() != md5:
            part.unlink(missing_ok=True)
            raise RuntimeError(f"{target.name}: MD5 doğrulaması başarısız")
        part.replace(target)

    def _md5(self, path: Path, desc: str) -> str:
        digest = hashlib.md5()
        with open(path, "rb") as f, self._bar(path.stat().st_size, desc) as bar:
            while chunk := f.read(1 << 20):
                digest.update(chunk)
                bar.update(len(chunk))
        return digest.hexdigest()

    def _log(self, message: str) -> None:
        if self.verbose:
            tqdm.write(message)

    def _bar(self, total: Optional[int], desc: str, unit: str = "B") -> tqdm:
        """İlerleme çubuğu; verbose False ise gösterilmez."""
        in_bytes = unit == "B"
        return tqdm(total=total or None, desc=desc, unit=unit, unit_scale=in_bytes,
                    unit_divisor=1024 if in_bytes else 1000, dynamic_ncols=True,
                    disable=not self.verbose)

    def _gzip_lines(self, source: Path, desc: str) -> Iterator[bytes]:
        """gzip dosyasının satırlarını verir; çubuk okunan sıkıştırılmış baytlara göre ilerler."""
        with open(source, "rb") as raw, gzip.GzipFile(fileobj=raw) as src, \
                self._bar(source.stat().st_size, desc) as bar:
            for n, line in enumerate(src):
                yield line
                if n % 50000 == 0:
                    bar.update(raw.tell() - bar.n)
            bar.update(raw.tell() - bar.n)

    def _file_lines(self, path: Path, desc: str) -> Iterator[str]:
        """Metin dosyasının satırlarını verir; çubuk okunan baytlara göre ilerler."""
        with open(path, "rb") as f, self._bar(path.stat().st_size, desc) as bar:
            done = 0
            for n, line in enumerate(f):
                done += len(line)
                yield line.decode()
                if n % 50000 == 0:
                    bar.update(done - bar.n)
            bar.update(done - bar.n)

    # --- Hazırlık (dizinleme) -------------------------------------------------

    def _build_genome(self, source: Path) -> None:
        """Genomu genome.fa olarak açar ve samtools faidx biçiminde genome.fai yazar.

        .fai satırı: ad, uzunluk, dizinin bayt konumu, satırdaki baz, satırdaki bayt.
        """
        target = self.release_dir / "genome.fa"
        rows = []
        current = None
        offset = 0
        with open(self._tmp(target), "wb") as out:
            for line in self._gzip_lines(source, "Genom açılıyor"):
                out.write(line)
                if line.startswith(b">"):
                    if current:
                        rows.append(current)
                    name = line[1:].split()[0].decode()
                    current = [name, 0, offset + len(line), 0, 0]
                elif current is not None:
                    bases = len(line.rstrip(b"\r\n"))
                    if current[3] == 0:
                        current[3], current[4] = bases, len(line)
                    current[1] += bases
                offset += len(line)
        if current:
            rows.append(current)
        fai = self.release_dir / "genome.fai"
        with open(self._tmp(fai), "w") as f:
            for row in rows:
                f.write("\t".join(map(str, row)) + "\n")
        self._commit(target, fai)

    def _build_fasta(self, source: Path, stem: Path, desc: str) -> None:
        """FASTA'yı <stem>.fa olarak açar ve <stem>.idx yazar.

        .idx satırı: kimlik, dizinin bayt konumu, dizinin bayt uzunluğu, başlık.
        """
        rows = []
        current = None
        offset = 0
        with open(self._tmp(stem.with_suffix(".fa")), "wb") as out:
            for line in self._gzip_lines(source, desc):
                out.write(line)
                if line.startswith(b">"):
                    if current:
                        current[2] = offset - current[1]
                        rows.append(current)
                    header = line[1:].rstrip(b"\r\n").decode()
                    current = [header.split("|")[0], offset + len(line), 0, header]
                offset += len(line)
        if current:
            current[2] = offset - current[1]
            rows.append(current)
        with open(self._tmp(stem.with_suffix(".idx")), "w") as f:
            for row in rows:
                f.write("\t".join(map(str, row)) + "\n")
        self._commit(stem.with_suffix(".fa"), stem.with_suffix(".idx"))

    def _build_tables(self, source: Path) -> None:
        """GTF'deki gene ve transcript satırlarından genes.tsv ve transcripts.tsv yazar."""
        genes_path = self.release_dir / "genes.tsv"
        transcripts_path = self.release_dir / "transcripts.tsv"
        genes = open(self._tmp(genes_path), "w")
        transcripts = open(self._tmp(transcripts_path), "w")
        with genes, transcripts:
            genes.write("\t".join(_GENES_COLUMNS) + "\n")
            transcripts.write("\t".join(_TRANSCRIPTS_COLUMNS) + "\n")
            for raw_line in self._gzip_lines(source, "GTF süzülüyor"):
                line = raw_line.decode()
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] not in ("gene", "transcript"):
                    continue
                attributes = dict(_ATTRIBUTE.findall(fields[8]))
                location = [fields[0], fields[3], fields[4], fields[6]]
                if fields[2] == "gene":
                    row = [attributes.get("gene_id", ""), attributes.get("gene_name", ""),
                           attributes.get("gene_type", "")] + location
                    genes.write("\t".join(row) + "\n")
                else:
                    row = [attributes.get("transcript_id", ""), attributes.get("gene_id", ""),
                           attributes.get("transcript_name", ""),
                           attributes.get("transcript_type", "")] + location
                    transcripts.write("\t".join(row) + "\n")
        self._commit(genes_path, transcripts_path)

    def _build_regions(self, source: Path) -> None:
        """pc_transcripts başlıklarındaki UTR5/CDS/UTR3 koordinatlarından regions.tsv yazar.

        Satır: transcript_id, bölgeler (ör. ``UTR5:1-367,CDS:368-2173,UTR3:2174-4922``).
        Koordinatlar transkript üzerinde 1 tabanlı ve iki ucu dahildir.
        """
        target = self.release_dir / "regions.tsv"
        with open(self._tmp(target), "w") as out:
            for line in self._gzip_lines(source, "CDS koordinatları çıkarılıyor"):
                if not line.startswith(b">"):
                    continue
                fields = line[1:].rstrip(b"\r\n").decode().split("|")
                regions = [f for f in fields if _REGION.fullmatch(f)]
                if regions:
                    out.write(f"{fields[0]}\t{','.join(regions)}\n")
        self._commit(target)

    @staticmethod
    def _tmp(path: Path) -> Path:
        return path.with_name(path.name + ".tmp")

    def _commit(self, *paths: Path) -> None:
        """.tmp dosyalarını sırayla gerçek adlarına taşır."""
        for path in paths:
            self._tmp(path).replace(path)

    # --- Dizinleri yükleme ----------------------------------------------------

    def _load(self) -> dict:
        """sys/ içindeki tabloları ve dizinleri belleğe bir kez yükler."""
        self._ensure()
        if self.__tables is not None:
            return self.__tables
        d = self.release_dir

        def read_tsv(name: str) -> List[Dict[str, str]]:
            lines = self._file_lines(d / name, f"Yükleniyor {name}")
            columns = next(lines).rstrip("\n").split("\t")
            return [dict(zip(columns, line.rstrip("\n").split("\t"))) for line in lines]

        genes = read_tsv("genes.tsv")
        transcripts = read_tsv("transcripts.tsv")

        def read_idx(name: str) -> Dict[str, Tuple[int, int, str]]:
            index = {}
            for line in self._file_lines(d / name, f"Yükleniyor {name}"):
                key, offset, size, header = line.rstrip("\n").split("\t", 3)
                index[key] = (int(offset), int(size), header)
            return index

        genome = {}
        with open(d / "genome.fai") as f:
            for line in f:
                name, length, offset, line_bases, line_bytes = line.split("\t")
                genome[name] = (int(length), int(offset), int(line_bases), int(line_bytes))

        tables = {
            "genes": {}, "gene_names": {}, "transcripts": {}, "transcript_names": {},
            "gene_transcripts": {}, "genome": genome,
            "transcript_seqs": read_idx("transcripts.idx"),
            "protein_seqs": read_idx("translations.idx"),
            "protein_of_transcript": {},
            "regions": {},
        }
        for line in self._file_lines(d / "regions.tsv", "Yükleniyor regions.tsv"):
            transcript_id, regions = line.rstrip("\n").split("\t")
            for key in {transcript_id, _unversioned(transcript_id)}:
                tables["regions"][key] = regions
        for g in genes:
            for key in {g["gene_id"], _unversioned(g["gene_id"])}:
                tables["genes"][key] = g
            tables["gene_names"].setdefault(g["gene_name"], []).append(g)
        for t in transcripts:
            for key in {t["transcript_id"], _unversioned(t["transcript_id"])}:
                tables["transcripts"][key] = t
            tables["transcript_names"][t["transcript_name"]] = t
            tables["gene_transcripts"].setdefault(t["gene_id"], []).append(t)
        for protein_id, (_, _, header) in tables["protein_seqs"].items():
            transcript_id = header.split("|")[1]
            for key in {transcript_id, _unversioned(transcript_id)}:
                tables["protein_of_transcript"][key] = protein_id
        for key in list(tables["transcript_seqs"]):
            tables["transcript_seqs"].setdefault(_unversioned(key), tables["transcript_seqs"][key])
        for key in list(tables["protein_seqs"]):
            tables["protein_seqs"].setdefault(_unversioned(key), tables["protein_seqs"][key])

        self.__tables = tables
        return tables

    def _find_gene(self, gene: str) -> Dict[str, str]:
        tables = self._load()
        if gene in tables["genes"]:
            return tables["genes"][gene]
        matches = tables["gene_names"].get(gene, [])
        if len(matches) > 1:
            matches = [g for g in matches if "_PAR_Y" not in g["gene_id"]] or matches
        if len(matches) == 1:
            return matches[0]
        if matches:
            ids = ", ".join(g["gene_id"] for g in matches)
            raise ValueError(f"'{gene}' adında birden çok gen var: {ids}; gene_id kullanın")
        raise KeyError(f"GENCODE {self.release}: '{gene}' geni bulunamadı")

    def _find_transcript(self, transcript: str) -> Dict[str, str]:
        tables = self._load()
        found = tables["transcripts"].get(transcript) or tables["transcript_names"].get(transcript)
        if found is None:
            raise KeyError(f"GENCODE {self.release}: '{transcript}' transkripti bulunamadı")
        return found

    def _read_record(self, file_name: str, entry: Tuple[int, int, str]) -> str:
        offset, size, _ = entry
        with open(self.release_dir / file_name, "rb") as f:
            f.seek(offset)
            return f.read(size).decode().replace("\n", "").replace("\r", "")

    def _features(self, obj, **values) -> None:
        obj.features.set("source", "GENCODE")
        obj.features.set("release", self.release)
        for key, value in values.items():
            if value not in (None, ""):
                obj.features.set(key, value)

    # --- Çıkarma -------------------------------------------------------------

    def region(self, chromosome: str, start: int, end: int, strand: str = "+",
               name: str = "") -> DNA:
        """Genomdan bir bölgeyi DNA olarak çıkarır.

        Args:
            chromosome: Kromozom adı (ör. ``"chr16"``).
            start, end: 1 tabanlı, iki ucu dahil konumlar (GTF ile aynı).
            strand: ``"+"`` referans zincir; ``"-"`` ise ters tümleyici alınır.
                Her iki durumda DNA 5'->3' yazılır.
            name: DNA'nın adı; boşsa ``chr16:3242028-3256627(-)`` biçiminde olur.

        Raises:
            KeyError: Kromozom yoksa.
            ValueError: Konumlar kromozomun dışındaysa ya da strand geçersizse.
        """
        genome = self._load()["genome"]
        if chromosome not in genome:
            raise KeyError(f"GENCODE {self.release}: '{chromosome}' kromozomu yok")
        if strand not in ("+", "-"):
            raise ValueError(f"strand '+' ya da '-' olmalı, {strand!r} verildi")
        length, offset, line_bases, line_bytes = genome[chromosome]
        if not 1 <= start <= end <= length:
            raise ValueError(f"{chromosome}:{start}-{end} kromozomun dışında (uzunluk {length})")

        def position(base: int) -> int:
            return offset + base // line_bases * line_bytes + base % line_bases

        first, last = position(start - 1), position(end - 1) + 1
        with open(self.release_dir / "genome.fa", "rb") as f:
            f.seek(first)
            seq = f.read(last - first).decode().replace("\n", "").replace("\r", "")
        if strand == "-":
            seq = _reverse_complement(seq)

        dna = DNA(seq, name=name or f"{chromosome}:{start}-{end}({strand})")
        self._features(dna, chromosome=chromosome, start=start, end=end,
                       chromosome_strand=strand)
        return dna

    def dna(self, gene: str) -> DNA:
        """Genin genomik DNA'sını kodlayan (sense) zincirde, 5'->3' döndürür.

        (-) zincirdeki genler için referans bölgenin ters tümleyicisi alınır.

        Args:
            gene: Gen adı (ör. ``"MEFV"``) ya da gene_id (sürümlü ya da sürümsüz).
        """
        g = self._find_gene(gene)
        dna = self.region(g["chromosome"], int(g["start"]), int(g["end"]), g["strand"],
                          name=g["gene_name"] or g["gene_id"])
        self._features(dna, gene_id=g["gene_id"], gene_name=g["gene_name"],
                       gene_type=g["gene_type"])
        return dna

    def rna(self, transcript: str) -> RNA:
        """Transkriptin (olgun, splice edilmiş) dizisini RNA olarak döndürür.

        Protein kodlayan transkriptlerde features'a ``utr5``, ``cds`` ve
        ``utr3`` bölgeleri ``(başlangıç, bitiş)`` olarak yazılır (transkript
        üzerinde 1 tabanlı, iki ucu dahil). CDS şöyle alınır::

            start, end = rna.features.first("cds")
            cds = rna[start - 1:end]

        Args:
            transcript: transcript_id (sürümlü ya da sürümsüz) ya da
                transkript adı (ör. ``"MEFV-201"``).
        """
        t = self._find_transcript(transcript)
        entry = self._load()["transcript_seqs"].get(t["transcript_id"])
        if entry is None:
            raise KeyError(f"GENCODE {self.release}: {t['transcript_id']} için dizi yok "
                           f"(yalnızca referans kromozomlardaki transkriptlerin dizisi vardır)")
        seq = self._read_record("transcripts.fa", entry).translate(Maps.CODING_TO_RNA_MAP.value)
        genomic_strand = GenomicStrand.Plus if t["strand"] == "+" else GenomicStrand.Minus
        rna = RNA(seq, name=t["transcript_id"], genomic_strand=genomic_strand)
        gene_name = self._load()["genes"].get(t["gene_id"], {}).get("gene_name")
        self._features(rna, transcript_id=t["transcript_id"], transcript_name=t["transcript_name"],
                       transcript_type=t["transcript_type"], gene_id=t["gene_id"],
                       gene_name=gene_name, chromosome=t["chromosome"],
                       start=int(t["start"]), end=int(t["end"]))
        for kind, start, end in _REGION.findall(self._load()["regions"].get(t["transcript_id"], "")):
            rna.features.add(kind.lower(), (int(start), int(end)))
        return rna

    def protein(self, identifier: str) -> PRT:
        """Proteini döndürür.

        Args:
            identifier: protein_id (ENSP), transcript_id (ENST) ya da transkript
                adı; sürüm numarası olmadan da verilebilir.

        Raises:
            KeyError: Protein yoksa (ör. kodlamayan transkript).
        """
        tables = self._load()
        entry = tables["protein_seqs"].get(identifier)
        if entry is None:
            transcript_id = identifier
            if identifier in tables["transcript_names"]:
                transcript_id = tables["transcript_names"][identifier]["transcript_id"]
            protein_id = tables["protein_of_transcript"].get(transcript_id)
            entry = tables["protein_seqs"].get(protein_id) if protein_id else None
        if entry is None:
            raise KeyError(f"GENCODE {self.release}: '{identifier}' için protein bulunamadı")

        fields = entry[2].split("|")
        prt = PRT(self._read_record("translations.fa", entry), name=fields[0])
        self._features(prt, protein_id=fields[0], transcript_id=fields[1], gene_id=fields[2],
                       transcript_name=fields[5], gene_name=fields[6])
        return prt

    def gene(self, gene: str) -> Gene:
        """Geni tüm içeriğiyle döndürür.

        Gene nesnesi şunları içerir: genomik DNA (``dna()``), genin dizisi olan
        tüm transkriptleri (RNA, kimlik = transcript_id) ve protein kodlayan
        transkriptlerin proteinleri (PRT, kimlik = protein_id). Proteinin hangi
        transkriptten geldiği proteinin ``transcript_id`` özelliğindedir.

        Args:
            gene: Gen adı (ör. ``"MEFV"``) ya da gene_id.
        """
        g = self._find_gene(gene)
        tables = self._load()
        result = Gene(g["gene_name"] or g["gene_id"], self.dna(g["gene_id"]))
        self._features(result, gene_id=g["gene_id"], gene_name=g["gene_name"],
                       gene_type=g["gene_type"], chromosome=g["chromosome"],
                       start=int(g["start"]), end=int(g["end"]), chromosome_strand=g["strand"])

        transcripts = tables["gene_transcripts"].get(g["gene_id"], [])
        for t in tqdm(transcripts, desc=f"{result.name} transkriptleri", unit="transkript",
                      dynamic_ncols=True, disable=not self.verbose):
            transcript_id = t["transcript_id"]
            if transcript_id not in tables["transcript_seqs"]:
                continue
            result.add_transcript(self.rna(transcript_id), transcript_id)
            if transcript_id in tables["protein_of_transcript"]:
                prt = self.protein(transcript_id)
                result.add_protein(prt, prt.name)
        return result
