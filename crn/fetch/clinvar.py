"""ClinVar: klinik varyantlardan gene ve hastalığa.

ClinVar'ın aylık ``variant_summary`` arşivini indirir, tek bir SQLite veritabanına
hazırlar ve varyantları klinik anlam, inceleme düzeyi, hastalık, varyant tipi,
protein etkisi ve konuma göre süzer. Süzülen varyantlardan gen listesi ve hastalık
listesi çıkarır; genleri ``Hpo`` ile aynı yoldan crn nesnelerine (``Gene``, ``DNA``,
``RNA``, ``PRT``) dönüştürür.

Sürüm, arşivin ayıdır (``"2026-11"``). Arşivler aylık tutulduğu için aynı sürüm
her zaman aynı veriyi verir; haftalık güncel dosya kullanılmaz.

``sys/`` klasörünün yapısı::

    sys/clinvar/<sürüm>/                       ör. sys/clinvar/2026-11/
        raw/variant_summary_<sürüm>.txt.gz     NCBI arşiv dosyası (yaklaşık 450 MB)
        clinvar_<derleme>.sqlite               hazırlanmış veritabanı (GRCh38 ya da GRCh37)
        manifest.json                          hazırlık tamamlandığında yazılır

Kullanım::

    from crn.fetch import Clinvar

    cv = Clinvar()                        # sys/ içindeki en yeni sürüm, GRCh38
    cv.update()                           # isteğe bağlı: indir ve hazırla (yoksa ilk sorguda yapılır)

    cv.variants("MEFV", significance="pathogenic", min_stars=2)
    cv.genes("Familial Mediterranean fever", significance=["pathogenic", "likely_pathogenic"])
    cv.symbols(condition="cardiomyopathy", effect="nonsense", min_stars=1)
    cv.conditions("mediterranean")        # hastalık ara
    cv.conditions(gene="MEFV")            # genin hastalıkları

    genler = cv.to_genes("Familial Mediterranean fever", significance="pathogenic")
    cv.annotate(Gencode().gene("MEFV"))   # features'a ClinVar özetini yazar
    cv.pack({"kardiyo": "cardiomyopathy", "epilepsi": "epilepsy"}, "gruplar.crn",
            significance=["pathogenic", "likely_pathogenic"], min_stars=1)

Filtreler (``variants`` ve onu kullanan tüm metodlar):

==================  =========================================================
 ``gene``           sembol ya da NCBI Gene ID (liste olabilir)
 ``condition``      hastalık adının bir parçası, ``re.compile`` deseni ya da
                    kimlik (``OMIM:249100``, ``MedGen:C0031069``, ``MONDO:0011122``)
 ``significance``   ``SIGNIFICANCES`` sınıfları; biri tutması yeter
 ``min_stars``      en az bu inceleme yıldızı (0-4)
 ``type``           ClinVar varyant tipi, ör. ``"single nucleotide variant"``
 ``origin``         ``germline``, ``somatic``, ``de novo`` ...; parça eşleşmesi
 ``effect``         ``EFFECTS`` protein etkileri (HGVS adından çıkarılır)
 ``chromosome``     kromozom; ``start`` ve ``end`` ile bölgeyle örtüşme
 ``max_genes``      en çok bu kadar gen içeren varyantlar (varsayılan 1); None sınırsız
==================  =========================================================

Yorumlarken bilinmesi gerekenler:

* Klinik anlam ClinVar'ın germline toplam sınıflandırmasıdır. ``Pathogenic/Likely
  pathogenic`` iki sınıfa birden girer; ``Conflicting classifications`` yalnızca
  ``conflicting`` sınıfına girer, patojen sayılmaz.
* Varyantın gene atanması ClinVar'ın gen listesine dayanır. Büyük kopya sayısı
  değişimleri onlarca gen listeler; ``max_genes=1`` varsayılanı bunları gen
  düzeyindeki sayımlardan çıkarır. ``LOC`` ile başlayan tanımlanmamış lokuslar
  sayılmaz.
* ``effect`` varyantın adındaki HGVS protein yazımından çıkarılır; protein değişimi
  yazılmamış varyantlarda boştur ve bir ardışık işlem hattının (VEP gibi) yerine geçmez.
* Bir varyant iki derlemede de konumlanmışsa yalnızca seçilen derlemenin konumu
  tutulur; yalnızca öteki derlemede konumu olan varyantların konumu o derlemedendir
  (``assembly`` sütunu). Bölge filtresi yalnızca seçilen derlemedeki konumlara uygulanır.
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import sqlite3
import urllib.error
import urllib.request
import zlib
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple

from crn.fetch._geneset import GeneSetSource
from crn.gene import Gene
from crn.seq import Sequence

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited"
ASSEMBLIES = ("GRCh38", "GRCh37")
SIGNIFICANCES = ("pathogenic", "likely_pathogenic", "uncertain", "likely_benign", "benign",
                 "conflicting", "risk", "drug_response", "other", "not_provided")
EFFECTS = ("missense", "nonsense", "synonymous", "frameshift", "inframe_indel",
           "start_lost", "stop_lost", "other")
SCHEMA_VERSION = 2

# ClinVar inceleme durumu -> yıldız (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
REVIEW_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the single variant": 0,
    "no classification for the individual variant": 0,
    "-": 0,
}

# Toplam sınıflandırmadaki terimler; değer "/" ve ";" ile birleşik gelir.
_SIGNIFICANCE_TERMS = {
    "pathogenic": "pathogenic",
    "pathogenic, low penetrance": "pathogenic",
    "likely pathogenic": "likely_pathogenic",
    "likely pathogenic, low penetrance": "likely_pathogenic",
    "uncertain significance": "uncertain",
    "vus-high": "uncertain",
    "vus-mid": "uncertain",
    "vus-low": "uncertain",
    "likely benign": "likely_benign",
    "benign": "benign",
    "conflicting classifications of pathogenicity": "conflicting",
    "conflicting interpretations of pathogenicity": "conflicting",
    "conflicting data from submitters": "conflicting",
    "risk factor": "risk",
    "established risk allele": "risk",
    "likely risk allele": "risk",
    "uncertain risk allele": "risk",
    "drug response": "drug_response",
    "association": "other",
    "association not found": "other",
    "protective": "other",
    "affects": "other",
    "confers sensitivity": "other",
    "other": "other",
    "not provided": "not_provided",
    "-": "not_provided",
    "no classification for the single variant": "not_provided",
    "no classification for the individual variant": "not_provided",
    "no classifications from unflagged records": "not_provided",
}
_PATHOGENIC = ("pathogenic", "likely_pathogenic")
_BENIGN = ("benign", "likely_benign")
# Hastalık sayılmayan koşul adları.
UNINFORMATIVE_CONDITIONS = frozenset({"not provided", "not specified", "see cases"})
_FEATURE_KEYS = ("clinvar_condition", "clinvar_variants", "clinvar_pathogenic")
_ASSEMBLY_RANK = {"GRCh38": 0, "GRCh37": 1, "na": 3, "NCBI36": 4}
_MONTHS = {name: index for index, name in enumerate(
    ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), 1)}
_NAME = re.compile(r"^(?P<tx>N[MR]_\d+(?:\.\d+)?)(?:\([^)]*\))?:(?P<c>[cn]\.\S+)"
                   r"(?:\s+\(p\.(?P<p>.+)\))?$")
_AA = r"[A-Z][a-z]{2}"
_CONDITION_ID = re.compile(r"^[A-Za-z][\w.]*:\S+$")
_REQUIRED = ("#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
             "LastEvaluated", "RS# (dbSNP)", "PhenotypeIDS", "PhenotypeList", "OriginSimple",
             "Assembly", "Chromosome", "Start", "Stop", "ReviewStatus", "NumberSubmitters",
             "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")
# Metin sütunları (tip, klinik anlam, inceleme durumu) labels tablosunda bir kez tutulur,
# koşullar conditions tablosunda; varyant satırı yalnızca kodları taşır. Sınıflar bit maskesidir.
_SCHEMA = """
CREATE TABLE variants (
    allele_id INTEGER PRIMARY KEY, variation_id INTEGER, type_code INTEGER, name TEXT,
    effect TEXT, gene_id TEXT, gene_symbol TEXT, n_genes INTEGER,
    significance_code INTEGER, classes INTEGER, stars INTEGER, review_code INTEGER,
    last_evaluated TEXT, origin TEXT, rs TEXT, submitters INTEGER,
    assembly TEXT, rank INTEGER, chromosome TEXT, start INTEGER, stop INTEGER,
    position_vcf INTEGER, ref TEXT, alt TEXT, oncogenicity TEXT, somatic_impact TEXT);
CREATE TABLE labels (code INTEGER PRIMARY KEY, text TEXT);
CREATE TABLE conditions (condition_id INTEGER PRIMARY KEY, name TEXT, ids TEXT);
CREATE TABLE variant_genes (allele_id INTEGER, symbol TEXT,
                            PRIMARY KEY (allele_id, symbol)) WITHOUT ROWID;
CREATE TABLE variant_conditions (allele_id INTEGER, condition_id INTEGER,
                                 PRIMARY KEY (allele_id, condition_id)) WITHOUT ROWID;
"""
_INDEXES = """
CREATE INDEX variant_genes_symbol ON variant_genes (symbol);
CREATE INDEX variant_conditions_condition ON variant_conditions (condition_id);
CREATE TABLE symbols AS
    SELECT g.symbol AS symbol, max(v.gene_id) AS gene_id, count(*) AS variants
    FROM variant_genes g JOIN variants v ON v.allele_id = g.allele_id GROUP BY g.symbol;
CREATE UNIQUE INDEX symbols_symbol ON symbols (symbol);
CREATE INDEX symbols_symbol_nocase ON symbols (symbol COLLATE NOCASE);
CREATE INDEX symbols_gene_id ON symbols (gene_id);
"""
_UPSERT = (f"INSERT INTO variants VALUES ({', '.join('?' * 26)}) ON CONFLICT (allele_id) DO UPDATE SET "
           "assembly = excluded.assembly, rank = excluded.rank, chromosome = excluded.chromosome, "
           "start = excluded.start, stop = excluded.stop, position_vcf = excluded.position_vcf, "
           "ref = excluded.ref, alt = excluded.alt WHERE excluded.rank < variants.rank")
_SELECT = """
SELECT v.allele_id, v.variation_id, v.gene_symbol, v.gene_id, v.n_genes, v.name, v.effect,
       lt.text AS type, ls.text AS significance, v.classes, v.stars, lr.text AS review_status,
       v.last_evaluated, v.origin, v.rs, v.submitters, v.assembly, v.chromosome, v.start, v.stop,
       v.position_vcf, v.ref, v.alt,
       (SELECT group_concat(c.name, '|') FROM variant_conditions vc
        JOIN conditions c ON c.condition_id = vc.condition_id WHERE vc.allele_id = v.allele_id) AS conditions,
       (SELECT group_concat(c.ids, '|') FROM variant_conditions vc
        JOIN conditions c ON c.condition_id = vc.condition_id WHERE vc.allele_id = v.allele_id) AS condition_ids
FROM variants v
LEFT JOIN labels lt ON lt.code = v.type_code
LEFT JOIN labels ls ON ls.code = v.significance_code
LEFT JOIN labels lr ON lr.code = v.review_code
"""
_ROW_ORDER = ("allele_id", "variation_id", "gene_symbol", "gene_id", "n_genes", "name", "transcript", "hgvs_c",
              "hgvs_p", "effect", "type", "significance", "classes", "stars", "review_status",
              "last_evaluated", "origin", "rs", "submitters", "assembly", "chromosome", "start",
              "stop", "position_vcf", "ref", "alt", "conditions", "condition_ids")

ConditionQuery = Any    # str | re.Pattern | Iterable
Filter = Any            # str | int | Iterable


def _as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple, set, frozenset)) else [value]


def significance_classes(text: str) -> List[str]:
    """ClinVar toplam sınıflandırmasını ``SIGNIFICANCES`` sınıflarına çevirir.

    ``"Pathogenic/Likely pathogenic; risk factor"`` -> ``["pathogenic",
    "likely_pathogenic", "risk"]``. Tanınmayan terim ``other`` olur; sıra ilk
    görülme sırasıdır.
    """
    classes: List[str] = []
    for term in re.split(r"/|;", text or "-"):
        term = term.strip().lower()
        if term:
            value = _SIGNIFICANCE_TERMS.get(term, "other")
            if value not in classes:
                classes.append(value)
    return classes or ["not_provided"]


def _mask(classes: Iterable[str]) -> int:
    """Sınıf adlarını bit maskesine çevirir (bit sırası ``SIGNIFICANCES``)."""
    return sum(1 << SIGNIFICANCES.index(name) for name in set(classes))


def _unmask(mask: int) -> List[str]:
    return [name for bit, name in enumerate(SIGNIFICANCES) if mask >> bit & 1]


def review_stars(status: str) -> int:
    """İnceleme durumunu ClinVar yıldızına (0-4) çevirir; tanınmayan durum 0."""
    return REVIEW_STARS.get((status or "-").strip().lower(), 0)


def protein_effect(hgvs_p: str) -> str:
    """HGVS protein yazımından (``p.`` sonrası) protein etkisi; yazım yoksa boş.

    Kurallar sırayla: ``=`` eşanlamlı, ``fs`` çerçeve kayması, ``ext`` ya da
    ``Ter``'in başka amino aside dönmesi stop kaybı, ``Ter``'le biten tek değişim ya
    da ``Ter`` içeren ekleme anlamsız (nonsense), ``Met1`` değişimi başlangıç kaybı,
    tek amino asit değişimi yanlış anlamlı (missense), ``del``/``ins``/``dup`` çerçeve
    içi indel; geri kalanı ``other``.
    """
    p = (hgvs_p or "").strip().strip("()")
    if not p or p == "?":
        return ""
    if p.endswith("="):
        return "synonymous"
    if "fs" in p:
        return "frameshift"
    if "ext" in p or re.fullmatch(rf"Ter\d+{_AA}", p):
        return "stop_lost"
    if re.fullmatch(rf"{_AA}\d+Ter", p) or re.search(r"ins(?:[A-Z][a-z]{2})*Ter", p):
        return "nonsense"
    if re.match(r"Met1(?!\d)", p):
        return "start_lost"
    if re.fullmatch(rf"{_AA}\d+{_AA}", p):
        return "missense"
    if re.search(r"del|ins|dup", p):
        return "inframe_indel"
    return "other"


def _date(text: str) -> str:
    """``"Dec 17, 2024"`` -> ``"2024-12-17"``; yerel ayardan bağımsız."""
    match = re.fullmatch(r"(\w{3}) (\d{1,2}), (\d{4})", text or "")
    if not match or match.group(1) not in _MONTHS:
        return ""
    month, day, year = _MONTHS[match.group(1)], int(match.group(2)), match.group(3)
    return f"{year}-{month:02d}-{day:02d}"


def _condition_ids(text: str) -> str:
    """Kimlikleri virgülle birleştirir; ``MONDO:MONDO:x`` yazımı ``MONDO:x`` olur."""
    ids = [item.strip().replace("MONDO:MONDO:", "MONDO:") for item in text.split(",")]
    return ",".join(item for item in ids if item and item != "-")


def _conditions(ids_text: str, names_text: str) -> Iterator[Tuple[str, str]]:
    """``PhenotypeIDS`` ve ``PhenotypeList``'ten (ad, kimlikler) çiftleri.

    RCV'ler ``|``, aynı RCV'deki koşullar ``;`` ile ayrılır. Koşul sayıları
    tutmazsa adlar kimliksiz döner; beşten çok koşul yerine yazılan sayılar atlanır.
    """
    id_groups = ids_text.split("|") if ids_text and ids_text != "-" else []
    for index, names in enumerate(names_text.split("|") if names_text else []):
        name_list = [name.strip() for name in names.split(";")]
        ids = id_groups[index].split(";") if index < len(id_groups) else []
        if len(ids) != len(name_list):
            ids = [""] * len(name_list)
        for name, identifiers in zip(name_list, ids):
            if name and name != "-" and not name.isdigit():
                yield name, _condition_ids(identifiers)


# Büyük CNV'lerde gen alanı serbest metindir: "subset of 102 genes: PRRT2:TBX1" yalnızca
# doz duyarlı genleri sayar, "covers 17 genes, none of which ..." hiçbir gen saymaz.
_CNV_GENES = re.compile(r"(?:subset of|covers) (?P<total>\d+) genes?(?:: (?P<listed>.*)|,.*)?$")


def _gene_symbols(text: str) -> List[str]:
    text = (text or "").strip()
    cnv = _CNV_GENES.match(text)
    parts = (cnv.group("listed") or "").split(":") if cnv else text.split(";")
    return [s for s in dict.fromkeys(part.strip() for part in parts) if s and s != "-"]


def _gene_count(symbols: List[str], text: str = "") -> int:
    """Varyantın kapsadığı gen sayısı.

    CNV metninde ClinVar'ın verdiği toplam kullanılır; öbür durumlarda ``LOC``
    lokusları sayılmaz, yalnızca onlar varsa sayıları döner.
    """
    cnv = _CNV_GENES.match((text or "").strip())
    if cnv:
        return int(cnv.group("total"))
    named = [s for s in symbols if not re.fullmatch(r"LOC\d+", s)]
    return len(named) or len(symbols)


class Clinvar(GeneSetSource):
    """ClinVar varyantları, gen ve hastalık sorguları, crn nesnelerine dönüşüm.

    Args:
        release: ``"2026-11"`` biçiminde arşiv ayı. None ise ``sys/`` içindeki en
            yeni hazır sürüm kullanılır; hiç yoksa ``update()`` en son arşivi indirir.
        assembly: Konumların alınacağı derleme, ``"GRCh38"`` ya da ``"GRCh37"``.
        sys_dir: Verilerin tutulacağı klasör; varsayılan proje kökündeki ``sys/``.
        verbose: True ise indirme ve hazırlık ilerleme çubuklarıyla (tqdm) gösterilir.
    """

    SOURCE = "ClinVar"
    FEATURE_KEYS = _FEATURE_KEYS
    RELEASE_KEY = "clinvar_release"
    GROUP_FILTER = "condition"

    def __init__(self, release: Optional[str] = None, assembly: str = "GRCh38",
                 sys_dir: Optional[str | Path] = None, verbose: bool = True):
        super().__init__(sys_dir, verbose)
        self._check_choices("assembly", [assembly], ASSEMBLIES)
        self.assembly = assembly
        self.__requested = self._normalize_release(release) if release else None
        self.__release = self.__requested
        self.__db: Optional[sqlite3.Connection] = None
        self.__patterns: Dict[str, re.Pattern] = {}

    # --- Sürüm ve klasörler -------------------------------------------------

    @staticmethod
    def _normalize_release(release: str) -> str:
        if not re.fullmatch(r"\d{4}-\d{2}", release):
            raise ValueError(f"ClinVar sürümü YYYY-AA biçiminde olmalı, {release!r} verildi")
        return release

    @property
    def clinvar_dir(self) -> Path:
        return self.sys_dir / "clinvar"

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
        return self.clinvar_dir / self.release

    @property
    def database(self) -> Path:
        """Seçili sürüm ve derlemenin SQLite dosyası."""
        return self.release_dir / f"clinvar_{self.assembly}.sqlite"

    def _manifest(self, folder: Path) -> Dict[str, Any]:
        path = folder / "manifest.json"
        return json.loads(path.read_text()) if path.is_file() else {}

    def _prepared(self, folder: Path) -> bool:
        info = self._manifest(folder).get("assemblies", {}).get(self.assembly, {})
        return (info.get("schema") == SCHEMA_VERSION
                and (folder / f"clinvar_{self.assembly}.sqlite").is_file())

    def installed_releases(self) -> List[str]:
        """sys/ içinde bu derleme için hazırlanmış sürümler, eskiden yeniye."""
        if not self.clinvar_dir.is_dir():
            return []
        return sorted(p.name for p in self.clinvar_dir.iterdir() if p.is_dir() and self._prepared(p))

    def is_ready(self) -> bool:
        """Kullanılan sürümün veritabanı bu derleme için hazır mı?"""
        return self.release is not None and self._prepared(self.release_dir)

    def _open(self, url: str, method: str = "GET"):
        request = urllib.request.Request(url, method=method, headers={"User-Agent": "crn"})
        return urllib.request.urlopen(request, timeout=60)

    def latest_release(self) -> str:
        """Sunucudaki en son aylık arşiv (ağ bağlantısı gerekir)."""
        with self._open(f"{CLINVAR_URL}/archive/") as response:
            html = response.read().decode("utf-8", "replace")
        releases = re.findall(r'variant_summary_(\d{4}-\d{2})\.txt\.gz"', html)
        if not releases:
            raise RuntimeError("ClinVar arşiv listesinde variant_summary bulunamadı")
        return max(releases)

    def _locate(self, release: str) -> Tuple[str, int]:
        """Arşiv dosyasının adresi ve boyutu; eski yıllar alt klasördedir."""
        name = f"variant_summary_{release}.txt.gz"
        for url in (f"{CLINVAR_URL}/archive/{name}", f"{CLINVAR_URL}/archive/{release[:4]}/{name}"):
            try:
                with self._open(url, "HEAD") as response:
                    return url, int(response.headers.get("Content-Length") or 0)
            except urllib.error.HTTPError as error:
                if error.code != 404:
                    raise RuntimeError(f"ClinVar arşivine ulaşılamadı (HTTP {error.code})") from error
        raise ValueError(f"ClinVar sürümü bulunamadı: {release}")

    # --- update ----------------------------------------------------------------

    def update(self, force: bool = False) -> str:
        """Sürümün arşivini indirir ve seçili derleme için veritabanını hazırlar.

        Nesne bir sürümle oluşturulduysa o sürüm (hazırsa ağa çıkılmaz),
        oluşturulmadıysa sunucudaki en son aylık arşiv alınır. NCBI arşiv dosyaları
        için özet (MD5) yayımlamaz; indirme boyutla, içerik gzip CRC32 denetimiyle
        doğrulanır ve dosyanın SHA-256'sı manifest'e yazılır. İndirilmiş dosya
        yeniden indirilmez; öteki derleme aynı dosyadan hazırlanır.

        Args:
            force: True ise dosya yeniden indirilir ve veritabanı yeniden kurulur.

        Returns:
            Hazırlanan sürüm (ör. ``"2026-11"``).
        """
        if self.__requested and self.is_ready() and not force:
            self._log(f"ClinVar {self.release} ({self.assembly}) zaten hazır: {self.database}")
            return self.release
        release = self.__requested or self.latest_release()
        self.__release = release
        self._close()
        if self.is_ready() and not force:
            self._log(f"ClinVar {release} ({self.assembly}) zaten hazır: {self.database}")
            return release

        url, size = self._locate(release)
        raw = self.release_dir / "raw" / url.rsplit("/", 1)[1]
        manifest = self._manifest(self.release_dir)
        if raw.is_file() and raw.stat().st_size == size and manifest.get("sha256") and not force:
            self._log(f"Var, atlandı: {raw.name}")
            digest = manifest["sha256"]
        else:
            digest = self._download(url, raw, size)
        try:
            stats = self._prepare(raw)
        except (EOFError, gzip.BadGzipFile, zlib.error) as error:
            raw.unlink(missing_ok=True)
            raise RuntimeError(f"{raw.name} bozuk ({error}); dosya silindi, update() "
                               f"yeniden çalıştırılmalı") from error

        manifest.update(release=release, source=url, size=size, sha256=digest)
        manifest.setdefault("assemblies", {})[self.assembly] = {
            "schema": SCHEMA_VERSION, **stats,
            "created": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        (self.release_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        self._log(f"ClinVar {release} ({self.assembly}) hazır: {self.database}")
        return release

    def _download(self, url: str, target: Path, size: int) -> str:
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + ".part")
        digest = hashlib.sha256()
        with self._open(url) as response, open(part, "wb") as out, \
                self._bar(size, f"İndiriliyor {target.name}") as bar:
            while chunk := response.read(1 << 20):
                out.write(chunk)
                digest.update(chunk)
                bar.update(len(chunk))
        if size and part.stat().st_size != size:
            part.unlink(missing_ok=True)
            raise RuntimeError(f"{target.name}: eksik indirme ({size} bayt bekleniyordu)")
        part.replace(target)
        return digest.hexdigest()

    def _records(self, raw: Path, bar: Any) -> Iterator[Dict[str, str]]:
        """Arşivdeki satırlar; ilerleme sıkıştırılmış baytlarla izlenir."""
        with open(raw, "rb") as handle, gzip.GzipFile(fileobj=handle) as unzipped:
            text = io.TextIOWrapper(unzipped, encoding="utf-8", newline="")
            header = text.readline().rstrip("\n").split("\t")
            missing = [column for column in _REQUIRED if column not in header]
            if missing:
                raise ValueError(f"{raw.name}: beklenen sütunlar yok: {', '.join(missing)}")
            optional = {"Oncogenicity", "SomaticClinicalImpact"}
            for count, line in enumerate(text, 1):
                values = line.rstrip("\n").split("\t")
                row = dict(zip(header, values))
                for column in optional:
                    row.setdefault(column, "-")
                yield row
                if count % 20000 == 0:
                    bar.update(handle.tell() - bar.n)
            bar.update(handle.tell() - bar.n)

    def _prepare(self, raw: Path) -> Dict[str, Any]:
        """Arşivi seçili derleme için SQLite veritabanına çevirir.

        Her varyant (AlleleID) bir kez yazılır; konum önce seçili derlemeden, yoksa
        öteki derlemeden alınır. Gen ve koşul satırları tekrarsızdır.
        """
        target = self.database
        part = target.with_name(target.name + ".part")
        part.unlink(missing_ok=True)
        db = sqlite3.connect(part)
        db.executescript("PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF; PRAGMA cache_size = -524288; PRAGMA temp_store = MEMORY;" + _SCHEMA)
        rank = {**_ASSEMBLY_RANK, self.assembly: -1}
        unknown_significance: Counter = Counter()
        unknown_review: Counter = Counter()
        labels: Dict[str, int] = {}
        condition_codes: Dict[Tuple[str, str], int] = {}
        variants, genes, links = [], [], []
        previous = None

        def label(text: str) -> int:
            return labels.setdefault(text, len(labels) + 1)

        def condition(name: str, ids: str) -> int:
            return condition_codes.setdefault((name, ids), len(condition_codes) + 1)

        def flush() -> None:
            db.executemany(_UPSERT, variants)
            db.executemany("INSERT OR IGNORE INTO variant_genes VALUES (?, ?)", genes)
            db.executemany("INSERT OR IGNORE INTO variant_conditions VALUES (?, ?)", links)
            variants.clear(), genes.clear(), links.clear()

        with self._bar(raw.stat().st_size, f"Hazırlanıyor {self.assembly}") as bar:
            for row in self._records(raw, bar):
                allele = int(row["#AlleleID"])
                significance = row["ClinicalSignificance"]
                for term in re.split(r"/|;", significance):
                    if term.strip() and term.strip().lower() not in _SIGNIFICANCE_TERMS:
                        unknown_significance[term.strip()] += 1
                status = row["ReviewStatus"]
                if status.strip().lower() not in REVIEW_STARS:
                    unknown_review[status] += 1
                match = _NAME.match(row["Name"])
                symbols = _gene_symbols(row["GeneSymbol"])
                variants.append((
                    allele, int(row["VariationID"] or 0), label(row["Type"]), row["Name"],
                    protein_effect((match.group("p") or "") if match else ""),
                    "" if row["GeneID"] in ("-1", "", "-") else row["GeneID"],
                    ";".join(symbols), _gene_count(symbols, row["GeneSymbol"]),
                    label(significance), _mask(significance_classes(significance)),
                    review_stars(status), label(status), _date(row["LastEvaluated"]),
                    row["OriginSimple"], "" if row["RS# (dbSNP)"] in ("-1", "-") else row["RS# (dbSNP)"],
                    int(row["NumberSubmitters"] or 0),
                    row["Assembly"], rank.get(row["Assembly"], 5), row["Chromosome"],
                    _int(row["Start"]), _int(row["Stop"]), _int(row["PositionVCF"]),
                    _allele(row["ReferenceAlleleVCF"]), _allele(row["AlternateAlleleVCF"]),
                    _dash(row["Oncogenicity"]), _dash(row["SomaticClinicalImpact"])))
                if allele != previous:      # aynı varyantın derleme satırları ardışıktır
                    genes.extend((allele, symbol) for symbol in symbols)
                    links.extend((allele, condition(name, ids)) for name, ids in
                                 _conditions(row["PhenotypeIDS"], row["PhenotypeList"]))
                    previous = allele
                if len(variants) >= 20000:
                    flush()
            flush()
        db.executemany("INSERT INTO labels VALUES (?, ?)", [(c, t) for t, c in labels.items()])
        db.executemany("INSERT INTO conditions VALUES (?, ?, ?)",
                       [(c, name, ids) for (name, ids), c in condition_codes.items()])
        db.executescript(_INDEXES)
        counts = dict(db.execute("SELECT 'variants', count(*) FROM variants UNION ALL "
                                 "SELECT 'genes', count(*) FROM symbols UNION ALL "
                                 "SELECT 'conditions', count(DISTINCT name) FROM conditions"))
        classes = {name: db.execute("SELECT count(*) FROM variants WHERE classes & ? != 0",
                                    (_mask([name]),)).fetchone()[0] for name in SIGNIFICANCES}
        db.commit()
        db.execute("VACUUM")
        db.close()
        part.replace(target)
        return {**counts, "classes": classes,
                "unknown_significance": dict(unknown_significance.most_common()),
                "unknown_review_status": dict(unknown_review.most_common())}

    def _ensure(self) -> None:
        """Veri hazır değilse indirir; sys/ içinde hazır sürüm varsa ağa çıkılmaz."""
        if not self.is_ready():
            self.update()

    def _connection(self) -> sqlite3.Connection:
        if self.__db is None:
            self._ensure()
            db = sqlite3.connect(f"file:{self.database}?mode=ro", uri=True)
            db.row_factory = sqlite3.Row
            db.create_function("crn_match", 2, self._regex_match, deterministic=True)
            self.__db = db
        return self.__db

    def _close(self) -> None:
        if self.__db is not None:
            self.__db.close()
            self.__db = None

    def _regex_match(self, key: str, text: Optional[str]) -> bool:
        return self.__patterns[key].search(text or "") is not None

    # --- Süzgeçler -------------------------------------------------------------

    def _symbols_for(self, gene: Filter) -> List[str]:
        """Sembol ya da NCBI Gene ID'leri ClinVar sembollerine çevirir."""
        db = self._connection()
        found: List[str] = []
        missing: List[str] = []
        for item in _as_list(gene):
            text = str(item).strip()
            column = "gene_id" if text.isdigit() else "symbol"
            rows = db.execute(f"SELECT symbol FROM symbols WHERE {column} = ? COLLATE NOCASE",
                              (text,)).fetchall()
            found += [row["symbol"] for row in rows]
            if not rows:
                missing.append(text)
        if missing:
            raise KeyError(f"ClinVar {self.release}: gen bulunamadı: {', '.join(missing)}")
        return list(dict.fromkeys(found))

    def _condition_clause(self, condition: ConditionQuery) -> Tuple[str, List[Any]]:
        """Koşul süzgecinin SQL parçası (``conditions`` tablosunun takma adı ``c``)."""
        parts, params = [], []
        for item in _as_list(condition):
            if isinstance(item, re.Pattern):
                key = str(len(self.__patterns))
                self.__patterns[key] = item
                parts.append("crn_match(?, c.name)")
                params.append(key)
            elif _CONDITION_ID.fullmatch(str(item)) and " " not in str(item):
                parts.append("instr(',' || c.ids || ',', ',' || ? || ',') > 0")
                params.append(_condition_ids(str(item)))
            else:
                parts.append("instr(lower(c.name), lower(?)) > 0")
                params.append(str(item))
        return "(" + " OR ".join(parts) + ")", params

    def _where(self, gene: Optional[Filter] = None, condition: Optional[ConditionQuery] = None,
               significance: Optional[str | Iterable[str]] = None, min_stars: Optional[int] = None,
               type: Optional[str | Iterable[str]] = None, origin: Optional[str | Iterable[str]] = None,
               effect: Optional[str | Iterable[str]] = None, chromosome: Optional[str] = None,
               start: Optional[int] = None, end: Optional[int] = None,
               max_genes: Optional[int] = 1) -> Tuple[str, List[Any], List[str]]:
        """Süzgeçlerin ``variants v`` üzerindeki WHERE ifadesi, parametreleri ve semboller."""
        clauses: List[str] = []
        params: List[Any] = []
        symbols: List[str] = []
        if gene is not None:
            symbols = self._symbols_for(gene)
            clauses.append("v.allele_id IN (SELECT allele_id FROM variant_genes WHERE symbol IN "
                           f"({', '.join('?' * len(symbols))}))")
            params += symbols
        if condition is not None:
            text, values = self._condition_clause(condition)
            clauses.append("v.allele_id IN (SELECT vc.allele_id FROM variant_conditions vc JOIN conditions c "
                           f"ON c.condition_id = vc.condition_id WHERE {text})")
            params += values
        classes = _as_list(significance)
        if classes:
            self._check_choices("significance", classes, SIGNIFICANCES)
            clauses.append("v.classes & ? != 0")
            params.append(_mask(classes))
        if min_stars is not None:
            if min_stars not in range(5):
                raise ValueError(f"min_stars 0 ile 4 arasında olmalı, {min_stars} verildi")
            clauses.append("v.stars >= ?")
            params.append(min_stars)
        types = _as_list(type)
        if types:
            clauses.append("v.type_code IN (SELECT code FROM labels WHERE lower(text) IN "
                           f"({', '.join('?' * len(types))}))")
            params += [str(item).lower() for item in types]
        origins = _as_list(origin)
        if origins:
            clauses.append("(" + " OR ".join("instr(lower(v.origin), lower(?)) > 0" for _ in origins) + ")")
            params += origins
        effects = _as_list(effect)
        if effects:
            self._check_choices("effect", effects, EFFECTS)
            clauses.append(f"v.effect IN ({', '.join('?' * len(effects))})")
            params += effects
        if chromosome is not None or start is not None or end is not None:
            if chromosome is None:
                raise ValueError("start ve end için chromosome gerekli")
            clauses.append("v.assembly = ? AND v.chromosome = ?")
            params += [self.assembly, str(chromosome).removeprefix("chr")]
            if start is not None:
                clauses.append("v.stop >= ?")
                params.append(int(start))
            if end is not None:
                clauses.append("v.start <= ?")
                params.append(int(end))
        if max_genes is not None:
            clauses.append("v.n_genes <= ?")
            params.append(int(max_genes))
        return (" AND ".join(clauses) or "1"), params, symbols

    # --- Sorgular --------------------------------------------------------------

    def variants(self, gene: Optional[Filter] = None, *, limit: Optional[int] = None,
                 **filters: Any) -> List[Dict[str, Any]]:
        """Filtrelere uyan varyantlar, AlleleID sırasıyla.

        Satır: ``allele_id, variation_id, gene_symbol, gene_id, n_genes, name,
        transcript, hgvs_c, hgvs_p, effect, type, significance, classes, stars, review_status,
        last_evaluated, origin, rs, submitters, assembly, chromosome, start, stop,
        position_vcf, ref, alt, conditions, condition_ids``. ``classes`` ve
        ``conditions`` ``|`` ile ayrılmıştır; ``gene_symbol`` çok genli varyantta
        ``;`` ile ayrılmıştır. ``n_genes`` varyantın kapsadığı gen sayısıdır; büyük
        CNV'lerde ClinVar yalnızca doz duyarlı genleri listelediğinden sembol
        sayısından büyük olabilir.

        Args:
            gene: Sembol ya da NCBI Gene ID (liste olabilir).
            limit: En çok bu kadar satır.
            **filters: Modül açıklamasındaki filtreler.

        Raises:
            KeyError: Gen ClinVar'da yoksa.
            ValueError: Filtre değeri geçersizse.
        """
        where, params, _ = self._where(gene, **filters)
        sql = f"{_SELECT} WHERE {where} ORDER BY v.allele_id"
        if limit is not None:
            sql += f" LIMIT {int(limit)}"
        rows = []
        for record in self._connection().execute(sql, params):
            row = dict(record)
            match = _NAME.match(row["name"])
            row.update(transcript=match.group("tx") if match else "",
                       hgvs_c=match.group("c") if match else "",
                       hgvs_p=(match.group("p") or "") if match else "",
                       classes="|".join(_unmask(row["classes"])),
                       conditions=row["conditions"] or "", condition_ids=row["condition_ids"] or "")
            rows.append({key: row[key] for key in _ROW_ORDER})
        return rows

    def genes(self, condition: Optional[ConditionQuery] = None, **filters: Any) -> List[Dict[str, Any]]:
        """Filtrelere uyan varyantların genleri; en çok varyantı olan önce.

        Satır: ``ncbi_gene_id, gene_symbol, variants, pathogenic, uncertain,
        benign, conflicting, max_stars``. ``pathogenic`` patojen ya da olası
        patojen, ``benign`` benign ya da olası benign sınıfına giren varyant
        sayısıdır; sayılar filtreden geçen varyantlar üzerindendir. ``gene``
        filtresi verilirse çok genli varyantların öteki genleri listelenmez.
        """
        where, params, symbols = self._where(condition=condition, **filters)
        if symbols:
            where += f" AND g.symbol IN ({', '.join('?' * len(symbols))})"
            params = params + symbols

        def has(*names: str) -> str:
            return f"sum(v.classes & {_mask(names)} != 0)"

        sql = (f"SELECT coalesce(s.gene_id, '') AS ncbi_gene_id, g.symbol AS gene_symbol, "
               f"count(*) AS variants, {has(*_PATHOGENIC)} AS pathogenic, {has('uncertain')} AS uncertain, "
               f"{has(*_BENIGN)} AS benign, {has('conflicting')} AS conflicting, max(v.stars) AS max_stars "
               "FROM variant_genes g JOIN variants v ON v.allele_id = g.allele_id "
               "LEFT JOIN symbols s ON s.symbol = g.symbol "
               f"WHERE {where} GROUP BY g.symbol ORDER BY variants DESC, g.symbol")
        return [dict(row) for row in self._connection().execute(sql, params)]

    def symbols(self, condition: Optional[ConditionQuery] = None, **filters: Any) -> List[str]:
        """Filtrelere uyan genlerin sembolleri, alfabetik."""
        return sorted(row["gene_symbol"] for row in self.genes(condition, **filters))

    def conditions(self, name: Optional[str | re.Pattern] = None, *, informative: bool = True,
                   limit: Optional[int] = None, **filters: Any) -> List[Dict[str, Any]]:
        """Koşullar (hastalıklar): adında ``name`` geçenler ve/veya filtrelere uyanlar.

        Satır: ``condition_name, condition_ids, variants, genes``; sayılar
        filtreden geçen varyantlar üzerindendir, en çok varyantı olan önce.

        Args:
            name: Adın bir parçası (büyük/küçük harf duyarsız), ``re.compile`` deseni
                ya da kimlik.
            informative: True ise ``not provided``, ``not specified`` ve ``see cases``
                çıkarılır.
            limit: En çok bu kadar satır.
            **filters: ``variants`` filtreleri, ör. ``gene="MEFV"``.
        """
        where, params, symbols = self._where(**filters)
        if name is not None:
            text, values = self._condition_clause(name)
            where += f" AND {text}"
            params = params + values
        if informative:
            where += f" AND lower(c.name) NOT IN ({', '.join('?' * len(UNINFORMATIVE_CONDITIONS))})"
            params = params + sorted(UNINFORMATIVE_CONDITIONS)
        gene_join = f" AND g.symbol IN ({', '.join('?' * len(symbols))})" if symbols else ""
        sql = ("SELECT c.name AS condition_name, max(c.ids) AS condition_ids, "
               "count(DISTINCT v.allele_id) AS variants, count(DISTINCT g.symbol) AS genes "
               "FROM variant_conditions vc JOIN conditions c ON c.condition_id = vc.condition_id "
               "JOIN variants v ON v.allele_id = vc.allele_id "
               f"LEFT JOIN variant_genes g ON g.allele_id = v.allele_id{gene_join} "
               f"WHERE {where} GROUP BY c.name ORDER BY variants DESC, c.name")
        if limit is not None:
            sql += f" LIMIT {int(limit)}"
        return [dict(row) for row in self._connection().execute(sql, symbols + params)]

    def summary(self) -> Dict[str, Any]:
        """Sürüm, derleme, kayıt sayıları, sınıf dağılımı ve tanınmayan terimler."""
        self._ensure()
        info = self._manifest(self.release_dir)["assemblies"][self.assembly]
        return {"release": self.release, "assembly": self.assembly, "variants": info["variants"],
                "genes": info["genes"], "conditions": info["conditions"], **info["classes"],
                "unknown_significance": sum(info["unknown_significance"].values()),
                "unknown_review_status": sum(info["unknown_review_status"].values())}

    # --- crn nesneleri ---------------------------------------------------------

    def annotate(self, obj: Gene | Sequence, gene: Optional[str] = None, **filters: Any) -> Gene | Sequence:
        """``Gene`` ya da diziye genin ClinVar özetini yazar.

        features'a ``clinvar_condition`` (bilgilendirici koşul adları, tekrarsız),
        ``clinvar_variants`` (süzgeçten geçen varyant sayısı), ``clinvar_pathogenic``
        (bunlardan patojen ya da olası patojen olanlar) ve ``clinvar_release``
        yazılır; aynı anahtarların eski değerleri silinir.

        Args:
            obj: ``Gene``, ``DNA``, ``RNA`` ya da ``PRT``.
            gene: Gen sembolü; verilmezse features'taki ``gene_name``, o da yoksa
                ``Gene`` nesnesinin adı kullanılır.
            **filters: ``variants`` filtreleri.

        Returns:
            Aynı nesne.

        Raises:
            ValueError: Gen sembolü belirlenemezse.
            KeyError: Gen ClinVar'da yoksa.
        """
        symbol = gene or obj.features.first("gene_name") or (obj.name if isinstance(obj, Gene) else None)
        if not symbol:
            raise ValueError("Gen sembolü belirlenemedi; gene= ile verin")
        rows = self.variants(symbol, **{k: v for k, v in filters.items() if k != "gene"})
        names = dict.fromkeys(name for row in rows for name in row["conditions"].split("|")
                              if name and name.lower() not in UNINFORMATIVE_CONDITIONS)
        for key in _FEATURE_KEYS:
            obj.features.remove(key)
        for name in names:
            obj.features.add("clinvar_condition", name)
        obj.features.set("clinvar_variants", len(rows))
        obj.features.set("clinvar_pathogenic",
                         sum(any(c in row["classes"].split("|") for c in _PATHOGENIC) for row in rows))
        obj.features.set(self.RELEASE_KEY, self.release)
        return obj

    def _describe(self, gene: Gene, row: Dict[str, Any], filters: Dict[str, Any]) -> None:
        """``_build`` kancası: genin süzgeçten geçen ClinVar özetini yazar."""
        # Gen filtresi satırın kendi genine daralır; özet yalnız o genin varyantlarıdır.
        rest = {k: v for k, v in filters.items() if k != "gene"}
        self.annotate(gene, gene=row["gene_symbol"], **rest)

    def to_genes(self, condition: Optional[ConditionQuery] = None, source: Any = None,
                 group: Optional[str] = None, key: str = "group", **filters: Any) -> List[Gene]:
        """Filtrelere uyan genleri ``Gene`` nesneleri olarak çıkarır.

        Her gen bir DNA, bir transkript ve bir protein içerir, yani doğrudan
        ``CUB``'a verilebilir. Genin ClinVar özeti features'a yazılır (bkz.
        ``annotate``). Protein kodlamayan ya da proteini transkriptle eşleşmeyen
        genler atlanır ve ``UserWarning`` ile bildirilir.

        Args:
            condition: Koşul adı parçası, desen ya da kimlik.
            source: ``Gencode`` (ana transkript) ya da ``Ncbi`` (NM_ öncelikli ilk
                izoform). None ise sys/ içindeki GENCODE.
            group: Verilirse her genin features'ına ``key`` anahtarıyla yazılır.
            key: Grup anahtarı.
            **filters: ``variants`` filtreleri.
        """
        return self._to_genes({"condition": condition, **filters}, source, group, key)

    def to_sequences(self, kind: str, condition: Optional[ConditionQuery] = None,
                     source: Any = None, **filters: Any) -> List[Sequence]:
        """Filtrelere uyan genlerin dizileri; ``write_fasta`` ile yazılabilir.

        Args:
            kind: ``"dna"`` (genomik), ``"rna"`` (ana transkript) ya da ``"protein"``.
            condition, source, **filters: ``to_genes`` ile aynı.

        Raises:
            ValueError: kind geçersizse.
        """
        return self._to_sequences(kind, {"condition": condition, **filters}, source)

    def __repr__(self) -> str:
        return f"Clinvar({self.release!r}, {self.assembly!r}, {str(self.clinvar_dir)!r})"


def _int(text: str) -> Optional[int]:
    return int(text) if text and text.lstrip("-").isdigit() and text != "-1" else None


def _allele(text: str) -> str:
    return "" if text in ("na", "-") else text


def _dash(text: str) -> str:
    return "" if text in ("-", "") else text

