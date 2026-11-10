"""HGNC gen tablosu ve HCOP ortolog tablosundan sorgulama.

Veriler yerel ``sys/`` klasöründe tutulur. Bir sorgu yapıldığında önce bu
klasöre bakılır; veri yoksa ``Hgnc.update()`` ile indirilir ve hazırlanır.

``sys/`` klasörünün yapısı::

    sys/hgnc/
        hgnc_complete_set.txt           HGNC'nin tüm genleri (sekmeyle ayrılmış)
        raw/human_<tür>_hcop_fifteen_column.txt.gz
                                        HCOP'tan indirilen tür başına ortologlar
        orthologs.tsv                   hazırlanan türlerin ortologları (varsayılan
                                        insana yakın 10 tür), tek tabloda;
                                        ``species`` ve ``taxon_id`` sütunlarıyla ayrılır
        manifest.json                   hazırlık tamamlandığında yazılır

Sorgularda her anahtar kelime bir sütundur; değer şunlardan biri olabilir:

* metin ya da sayı: hücrenin ``|`` ile ayrılmış parçalarından biri eşitse
  (büyük/küçük harf duyarsız), ör. ``gene_group="Aquaporins"``
* liste/küme: öğelerinden biri tutarsa, ör. ``symbol=["AQP1", "AQP2"]``
* ``re.compile(...)``: parçalardan birinde desen bulunursa
* fonksiyon: hücrenin tüm metniyle çağrılır, True dönerse

Kullanım::

    from crn.fetch import Hgnc

    hgnc = Hgnc()
    hgnc.update()                                   # isteğe bağlı: indir/hazırla
    hgnc.query(gene_group="Aquaporins", fields=["symbol", "hgnc_id", "entrez_id"])
    kodlayan = hgnc.symbols(locus_type="gene with protein product", status="Approved")
    hgnc.orthologs(human_symbol="AQP1", species=["mouse", "rat"])
    hgnc.orthologs(human_symbol="AQP1", taxid=9031)   # tavuk; ilk istekte eklenir
    hgnc.values("locus_type")                       # sütundaki değerler ve sayıları
"""
from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import io
import json
import re
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence

from crn.fetch._store import LocalSource, atomic_write, parse_taxid


HGNC_URL = "https://storage.googleapis.com/public-download-files"

# HCOP'un sunduğu bütün türler: dosya adındaki tür -> NCBI Taxonomy ID. Taxid'ler
# HCOP'un birleşik dosyasındaki (human_all_hcop_sixteen_column) değerlerdir.
HCOP_SPECIES = {
    "chimpanzee":   "9598",
    "macaque":      "9544",
    "mouse":        "10090",
    "rat":          "10116",
    "dog":          "9615",
    "cat":          "9685",
    "cattle":       "9913",
    "pig":          "9823",
    "horse":        "9796",
    "opossum":      "13616",
    "platypus":     "9258",
    "chicken":      "9031",
    "anole_lizard": "28377",
    "xenopus":      "8364",
    "zebrafish":    "7955",
    "fruitfly":     "7227",
    "c.elegans":    "6239",
    "s.cerevisiae": "4932",
    "s.pombe":      "284812",
}

# HCOP'un suş düzeyinde kullandığı türlerin tür düzeyi taxid'leri ve NCBI Gene
# kayıtlarının bağlı olduğu suşlar; bu taxid'ler de aynı HCOP türüne çıkar.
HCOP_TAXID_ALIASES = {"4896": "s.pombe", "559292": "s.cerevisiae"}

# update()'in varsayılan olarak hazırladığı, insana yakın 10 memeli. Öteki HCOP
# türleri istendiğinde (orthologs(taxid=...) ya da update(species=...)) eklenir.
ORTHOLOG_SPECIES = {
    "chimpanzee": "9598",
    "macaque":    "9544",
    "mouse":      "10090",
    "rat":        "10116",
    "dog":        "9615",
    "cat":        "9685",
    "cattle":     "9913",
    "pig":        "9823",
    "horse":      "9796",
    "opossum":    "13616",
}

# HCOP 15 sütunlu dosyalarda 8-14. sütunlar türün adıyla başlar (ör.
# mouse_symbol, mgi_id); tek tabloda hepsi aynı adı taşısın diye yeniden adlandırılır.
_HCOP_HUMAN = ["human_entrez_gene", "human_ensembl_gene", "hgnc_id", "human_name",
               "human_symbol", "human_chr", "human_assert_ids"]
_HCOP_ORTHOLOG = ["ortholog_entrez_gene", "ortholog_ensembl_gene", "ortholog_db_id",
                  "ortholog_name", "ortholog_symbol", "ortholog_chr", "ortholog_assert_ids"]
ORTHOLOG_COLUMNS = ["species", "taxon_id", *_HCOP_HUMAN, *_HCOP_ORTHOLOG, "support"]

_GENES = "hgnc_complete_set.txt"
_ORTHOLOGS = "orthologs.tsv"
_OUTPUTS = [_GENES, _ORTHOLOGS, "manifest.json"]

Filter = Any  # str | int | Iterable | re.Pattern | Callable[[str], bool]


def _values(value: Any) -> List[Any]:
    """Tek değeri listeye sarar; liste, küme ve demetleri listeye çevirir."""
    return [value] if isinstance(value, (str, int)) else list(value)


def _matches(cell: str, wanted: Filter) -> bool:
    """Hücre, sorgu değerine uyuyor mu (kurallar modül açıklamasında)."""
    if callable(wanted) and not isinstance(wanted, re.Pattern):
        return bool(wanted(cell))
    if isinstance(wanted, (list, tuple, set, frozenset)):
        return any(_matches(cell, w) for w in wanted)
    parts = cell.split("|") if cell else [""]
    if isinstance(wanted, re.Pattern):
        return any(wanted.search(p) for p in parts)
    target = str(wanted).casefold()
    return any(p.casefold() == target for p in parts)


class Hgnc(LocalSource):
    """HGNC gen tablosunda ve HCOP ortolog tablosunda sütun sütun sorgu.

    Args:
        sys_dir: Verilerin tutulacağı klasör; varsayılan proje kökündeki ``sys/``.
        base_url: HGNC dosyalarının sunucusu.
        verbose: True ise indirme ve hazırlık ilerleme çubuklarıyla gösterilir.
    """

    SOURCE = "HGNC"

    def __init__(self, sys_dir: Optional[str | Path] = None,
                 base_url: str = HGNC_URL, verbose: bool = True):
        super().__init__(sys_dir, verbose)
        self.base_url = base_url.rstrip("/")
        self.__columns: Optional[List[str]] = None

    @property
    def hgnc_dir(self) -> Path:
        return self.sys_dir / "hgnc"

    def is_ready(self) -> bool:
        """sys/ içinde gen ve ortolog tabloları hazır mı?"""
        return all((self.hgnc_dir / name).is_file() for name in _OUTPUTS)

    # --- update ----------------------------------------------------------------

    def prepared_species(self) -> Dict[str, str]:
        """Ortolog tablosunda hazır olan türler (ad -> taxid); hazır değilse boş."""
        manifest = self.hgnc_dir / "manifest.json"
        if not (self.is_ready() and manifest.is_file()):
            return {}
        return dict(json.loads(manifest.read_text()).get("species", {}))

    @staticmethod
    def species_name(value: str | int) -> str:
        """HCOP tür adını ya da taxid'i HCOP tür adına çevirir (``9031`` -> ``"chicken"``).

        Raises:
            ValueError: Tür HCOP'ta yoksa; HCOP yalnızca ``HCOP_SPECIES`` türlerini içerir.
        """
        number = parse_taxid(value)
        if number is None and value in HCOP_SPECIES:
            return str(value)
        found = next((name for name, taxid in HCOP_SPECIES.items() if taxid == number),
                     HCOP_TAXID_ALIASES.get(number or ""))
        if found is None:
            raise ValueError(f"HCOP bu türü içermiyor: {value}; desteklenen türler (taxid): "
                             + ", ".join(f"{n} ({t})" for n, t in HCOP_SPECIES.items()))
        return found

    def update(self, force: bool = False, species: Optional[Iterable[str | int]] = None) -> Path:
        """HGNC gen tablosunu ve HCOP ortologlarını indirip hazırlar.

        Varsayılan olarak ``ORTHOLOG_SPECIES`` içindeki 10 tür hazırlanır;
        daha önce eklenmiş türler korunur. İndirilen her dosya sunucunun MD5'i
        ile doğrulanır. Var olan dosyalar yeniden indirilmez; ``orthologs.tsv``
        yalnızca tür kümesi değiştiğinde yeniden üretilir.

        Args:
            force: True ise her şey yeniden indirilir ve hazırlanır.
            species: Ayrıca eklenecek HCOP türleri; ad ya da taxid
                (``["chicken", 7955]``).

        Returns:
            Verilerin bulunduğu klasör.

        Raises:
            ValueError: İstenen tür HCOP'ta yoksa.
        """
        wanted = {name: HCOP_SPECIES[name] for name in
                  [*ORTHOLOG_SPECIES, *self.prepared_species(),
                   *(self.species_name(s) for s in (species or []))]}
        if self.is_ready() and not force and set(wanted) <= set(self.prepared_species()):
            self._log(f"HGNC zaten hazır: {self.hgnc_dir}")
            return self.hgnc_dir
        changed = force or set(wanted) != set(self.prepared_species())
        (self.hgnc_dir / "manifest.json").unlink(missing_ok=True)

        genes = self.hgnc_dir / _GENES
        self._download(f"{self.base_url}/hgnc/tsv/tsv/{_GENES}", genes, force)
        raws = []
        for name in wanted:
            file_name = f"human_{name}_hcop_fifteen_column.txt.gz"
            raws.append((name, self.hgnc_dir / "raw" / file_name))
            self._download(f"{self.base_url}/hcop/{file_name}", raws[-1][1], force)

        orthologs = self.hgnc_dir / _ORTHOLOGS
        if changed or not orthologs.is_file():
            self._merge_orthologs(raws, orthologs)

        manifest = {"source": self.base_url,
                    "species": wanted,
                    "prepared": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        atomic_write(self.hgnc_dir / "manifest.json", json.dumps(manifest, indent=2) + "\n")
        self.__columns = None
        return self.hgnc_dir

    def _download(self, url: str, target: Path, force: bool) -> None:
        if target.exists() and not force:
            self._log(f"Var, atlandı: {target.name}")
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        part = target.with_name(target.name + ".part")
        digest = hashlib.md5()
        with urllib.request.urlopen(url, timeout=60) as response, open(part, "wb") as out:
            hashes = ",".join(response.headers.get_all("x-goog-hash") or [])
            expected = re.search(r"md5=([A-Za-z0-9+/=]+)", hashes)
            total = int(response.headers.get("Content-Length") or 0)
            with self._bar(total, f"İndiriliyor {target.name}") as bar:
                while chunk := response.read(1 << 20):
                    out.write(chunk)
                    digest.update(chunk)
                    bar.update(len(chunk))
        if expected and base64.b64encode(digest.digest()).decode() != expected.group(1):
            part.unlink(missing_ok=True)
            raise RuntimeError(f"{target.name}: MD5 doğrulaması başarısız")
        part.replace(target)

    def _merge_orthologs(self, raws: Sequence[tuple], target: Path) -> None:
        """Tür başına HCOP dosyalarını ortak sütun adlarıyla tek tabloda birleştirir."""
        buffer = io.StringIO()
        writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
        writer.writerow(ORTHOLOG_COLUMNS)
        for species, raw in self._progress(raws, "Ortologlar birleştiriliyor", "tür"):
            with gzip.open(raw, "rt", newline="") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader)
                if len(header) != 15 or header[:7] != _HCOP_HUMAN:
                    raise RuntimeError(f"{raw.name}: beklenmeyen HCOP sütunları: {header}")
                for row in reader:
                    writer.writerow([species, HCOP_SPECIES[species], *row])
        atomic_write(target, buffer.getvalue())

    # --- Sorgu -----------------------------------------------------------------------

    def _rows(self, name: str) -> Iterator[Dict[str, str]]:
        if not self.is_ready():
            self.update()
        with open(self.hgnc_dir / name, newline="") as f:
            yield from csv.DictReader(f, delimiter="\t")

    @property
    def columns(self) -> List[str]:
        """HGNC gen tablosunun sütunları."""
        if self.__columns is None:
            self.__columns = list(next(self._rows(_GENES)))
        return self.__columns

    @staticmethod
    def _select(rows: Iterable[Dict[str, str]], columns: Sequence[str],
                fields: Optional[Sequence[str]], filters: Dict[str, Filter]
                ) -> List[Dict[str, str]]:
        unknown = [c for c in [*filters, *(fields or [])] if c not in columns]
        if unknown:
            raise ValueError(f"Bilinmeyen sütun: {', '.join(unknown)}; "
                             f"sütunlar: {', '.join(columns)}")
        result = []
        for row in rows:
            if all(_matches(row[c], w) for c, w in filters.items() if w is not None):
                result.append({c: row[c] for c in fields} if fields else row)
        return result

    def query(self, fields: Optional[Sequence[str]] = None, **filters: Filter
              ) -> List[Dict[str, str]]:
        """HGNC gen tablosunda filtrelere uyan satırlar.

        Args:
            fields: Döndürülecek sütunlar; None ise tüm sütunlar.
            **filters: sütun=değer çiftleri; hepsi birlikte sağlanmalı (VE).

        Örn. ``query(gene_group="Aquaporins", fields=["symbol", "entrez_id"])``.
        """
        return self._select(self._rows(_GENES), self.columns, fields, filters)

    def symbols(self, **filters: Filter) -> List[str]:
        """Filtrelere uyan genlerin onaylı sembolleri."""
        return [row["symbol"] for row in self.query(fields=["symbol"], **filters)]

    def values(self, column: str, **filters: Filter) -> Dict[str, int]:
        """Bir sütundaki (``|`` ile ayrılmış) değerler ve gen sayıları, çoktan aza.

        Hangi değerlerle sorgu yapılabileceğini görmek için kullanılır.
        """
        counts = Counter(part for row in self.query(fields=[column], **filters)
                         for part in row[column].split("|") if part)
        return dict(counts.most_common())

    def orthologs(self, fields: Optional[Sequence[str]] = None, **filters: Filter
                  ) -> List[Dict[str, str]]:
        """HCOP ortolog tablosunda filtrelere uyan satırlar.

        Sütunlar ``ORTHOLOG_COLUMNS``. Tür ``species`` (HCOP adı) ya da ``taxid``
        (NCBI Taxonomy ID, ``taxon_id`` sütunu) ile seçilir; ikisi de tek değer
        ya da liste olabilir. İstenen HCOP türü henüz hazırlanmamışsa dosyası
        indirilip tabloya eklenir.

        Örn. ``orthologs(human_symbol="AQP1", species=["mouse", "rat"])`` ya da
        ``orthologs(human_symbol="AQP1", taxid=[9031, 7955])``.

        Raises:
            ValueError: Tür HCOP'ta yoksa ya da sütun bilinmiyorsa.
        """
        taxid = filters.pop("taxid", None)
        if taxid is not None:
            names = [self.species_name(v) for v in _values(taxid)]
            filters["species"] = [*_values(filters.get("species")), *names] \
                if filters.get("species") is not None else names
        species = filters.get("species")
        if species is not None and not callable(species) and not isinstance(species, re.Pattern):
            names = [self.species_name(s) for s in _values(species)]
            missing = [n for n in names if n not in self.prepared_species()]
            if missing:
                self.update(species=missing)
            filters["species"] = names
        return self._select(self._rows(_ORTHOLOGS), ORTHOLOG_COLUMNS, fields, filters)

    @staticmethod
    def write_tsv(rows: Sequence[Dict[str, str]], file_name: str | Path) -> None:
        """Sorgu sonucunu sekmeyle ayrılmış dosyaya yazar."""
        if not rows:
            raise ValueError("Yazılacak satır yok")
        with open(file_name, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def __repr__(self) -> str:
        return f"Hgnc({str(self.hgnc_dir)!r})"
