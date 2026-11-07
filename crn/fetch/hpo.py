"""Human Phenotype Ontology (HPO): semptomdan gen grubuna.

Bir semptom ya da fenotip teriminden gen listesi çıkarır, listeyi filtreler ve
crn nesnelerine (``Gene``, ``DNA``, ``RNA``, ``PRT``) dönüştürür. Ontoloji
sayesinde bir terim alt terimleriyle birlikte sorgulanır: "Seizure" istendiğinde
"Focal seizure" ile işaretlenmiş genler de gelir.

Veriler yerel ``sys/`` klasöründe tutulur. Bir sorgu yapıldığında önce bu
klasöre bakılır; veri yoksa ``Hpo.update()`` ile indirilir.

``sys/`` klasörünün yapısı::

    sys/hpo/<sürüm>/                    ör. sys/hpo/v2026-09-01/
        raw/hp.obo                      ontoloji: terimler, eş anlamlılar, hiyerarşi
        raw/genes_to_phenotype.txt      gen, terim, sıklık, hastalık
        raw/genes_to_disease.txt        gen, hastalık, ilişki tipi
        raw/phenotype.hpoa              hastalık notları: kanıt, başlangıç, kalıtım
        manifest.json                   hazırlık tamamlandığında yazılır

Kullanım::

    from crn.fetch import Hpo

    hpo = Hpo()                         # sys/ içindeki en yeni sürüm
    hpo.update()                        # isteğe bağlı: indir (yoksa ilk sorguda yapılır)

    hpo.search("seizure")               # terim ara (ad ve eş anlamlılar)
    hpo.term("Seizure")                 # tek terim: tanım, üst ve alt terimler
    hpo.symbols("HP:0001250", min_frequency=0.3, association="MENDELIAN")
    hpo.genes("Seizure", exclude="Autism")          # gen başına özet
    hpo.diseases("parkinson")           # hastalık ara
    hpo.phenotypes("SCN1A")             # genin fenotipleri

    genler = hpo.to_genes("Ataxia", min_frequency=0.8)     # List[Gene], CUB'a hazır
    proteinler = hpo.to_sequences("protein", "Ataxia")     # List[PRT]
    hpo.annotate(Gencode().gene("SCN1A"))                  # features'a terimleri yazar

    rapor = hpo.pack({"nobet": "Seizure", "ataksi": "Ataxia"}, "gruplar.crn")
    DifferentialCodonUsage.from_crn("gruplar.crn", "nobet", "ataksi")

Filtreler (``annotations`` ve onu kullanan tüm metodlar):

==================  =========================================================
 ``term``           terim ya da terimler (birleşim); alt terimler dahildir
 ``exclude``        bu terimlerle ilişkili genleri sonuçtan çıkarır
 ``gene``           sembol ya da NCBI Gene ID
 ``disease``        OMIM/ORPHA kimliği, hastalığın tam adı ya da ``re.compile``
 ``min_frequency``  en az bu sıklık (0-1); sıklığı bilinmeyen satırlar elenir
 ``association``    ``MENDELIAN``, ``POLYGENIC``, ``UNKNOWN``
 ``evidence``       ``PCS`` (yayın), ``TAS`` (uzman), ``IEA`` (otomatik)
 ``onset``          hastalığın başlangıç terimi, alt terimleriyle
 ``inheritance``    hastalığın kalıtım terimi, alt terimleriyle
 ``descendants``    False ise yalnızca terimin kendisiyle işaretli satırlar
==================  =========================================================

Yorumlarken bilinmesi gerekenler:

* HPO'da gen-fenotip ilişkisi hastalık üzerindendir: gen bir hastalığa,
  hastalık fenotipe bağlıdır. Genin o fenotipi tek başına doğrudan yaptığı
  anlamına gelmez; hastalığın başka bir özelliği olabilir.
* ``genes_to_phenotype.txt`` hastalıkta bulunmadığı açıkça yazılmış (``NOT``)
  notları ve hiç gözlenmemiş (``0/N``) satırları da içerir. Bunlar yüklenirken
  çıkarılır; sayıları ``summary()`` içindedir.
* Sıklık kategorilerinde alt sınır kullanılır: "Frequent" (%30-79) 0.30 sayılır,
  yani ``min_frequency`` geçilen satırın o sıklığa ulaştığını garanti eder.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.error
import urllib.request
import warnings
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence as Seq, Set, Tuple

from crn import fileio
from crn.fetch._store import LocalSource
from crn.fetch.gencode import AmbiguousGeneError, Gencode
from crn.fetch.hgnc import Hgnc
from crn.gene import Gene
from crn.seq import Sequence


HPO_RELEASES = "https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases"

PHENOTYPE_ROOT = "HP:0000118"       # Phenotypic abnormality
ONSET_ROOT = "HP:0003674"           # Onset
ASSOCIATIONS = ("MENDELIAN", "POLYGENIC", "UNKNOWN")
EVIDENCES = ("PCS", "TAS", "IEA")
SEQUENCE_KINDS = ("dna", "rna", "protein")

# HPO sıklık terimlerinin tanımındaki alt sınırlar.
FREQUENCY_TERMS = {"HP:0040280": 1.0, "HP:0040281": 0.80, "HP:0040282": 0.30,
                   "HP:0040283": 0.05, "HP:0040284": 0.01, "HP:0040285": 0.0}

_OBO, _GENE_PHENOTYPE, _GENE_DISEASE, _HPOA = ("hp.obo", "genes_to_phenotype.txt",
                                               "genes_to_disease.txt", "phenotype.hpoa")
_FILES = (_OBO, _GENE_PHENOTYPE, _GENE_DISEASE, _HPOA)
_SYNONYM = re.compile(r'^"(.*)" (\w+)')
_DEFINITION = re.compile(r'^"(.*)"')
_NO_SYMBOL = "-"
_FEATURE_KEYS = ("hpo_id", "hpo_name", "disease_id", "disease_name")

TermQuery = Any         # str | Iterable[str]
Filter = Any            # str | int | Iterable | re.Pattern


def _as_list(value: Any) -> List[Any]:
    if value is None:
        return []
    return list(value) if isinstance(value, (list, tuple, set, frozenset)) else [value]


def frequency_value(frequency: str) -> Optional[float]:
    """HPO sıklık yazımını 0-1 arasına çevirir; bilinmiyorsa None.

    ``"7/12"``, ``"45%"`` ve ``"HP:0040282"`` yazımlarını tanır. Kategorilerde
    alt sınır döner.
    """
    if frequency in FREQUENCY_TERMS:
        return FREQUENCY_TERMS[frequency]
    if "/" in frequency:
        observed, total = frequency.split("/")
        return int(observed) / int(total) if int(total) else None
    if frequency.endswith("%"):
        return float(frequency[:-1]) / 100
    return None


class Hpo(LocalSource):
    """HPO terimleri, gen-fenotip notları ve crn nesnelerine dönüşüm.

    Args:
        release: ``"v2026-09-01"`` ya da ``"2026-09-01"``. None ise ``sys/``
            içindeki en yeni hazır sürüm kullanılır; hiç yoksa ``update()`` en
            son sürümü indirir.
        sys_dir: Verilerin tutulacağı klasör; varsayılan proje kökündeki ``sys/``.
        verbose: True ise indirme ve yükleme ilerleme çubuklarıyla (tqdm) gösterilir.
    """

    SOURCE = "HPO"

    def __init__(self, release: Optional[str] = None, sys_dir: Optional[str | Path] = None,
                 verbose: bool = True):
        super().__init__(sys_dir, verbose)
        self.__requested = self._normalize_release(release) if release else None
        self.__release = self.__requested
        self.__data: Optional[dict] = None
        self.__gencode: Optional[Gencode] = None

    # --- Sürüm ve klasörler -------------------------------------------------

    @staticmethod
    def _normalize_release(release: str) -> str:
        return release if release.startswith("v") else f"v{release}"

    @property
    def hpo_dir(self) -> Path:
        return self.sys_dir / "hpo"

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
        return self.hpo_dir / self.release

    def installed_releases(self) -> List[str]:
        """sys/ içinde hazırlığı tamamlanmış sürümler, eskiden yeniye."""
        if not self.hpo_dir.is_dir():
            return []
        return sorted(p.name for p in self.hpo_dir.iterdir() if (p / "manifest.json").is_file())

    def is_ready(self) -> bool:
        """Kullanılan sürümün tüm dosyaları sys/ içinde hazır mı?"""
        return (self.release is not None
                and (self.release_dir / "manifest.json").is_file()
                and all((self.release_dir / "raw" / name).is_file() for name in _FILES))

    def _release_info(self, release: Optional[str]) -> dict:
        url = f"{HPO_RELEASES}/tags/{release}" if release else f"{HPO_RELEASES}/latest"
        request = urllib.request.Request(url, headers={"User-Agent": "crn",
                                                       "Accept": "application/vnd.github+json"})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise ValueError(f"HPO sürümü bulunamadı: {release}") from error
            raise RuntimeError(f"HPO sürüm bilgisi alınamadı (HTTP {error.code}); "
                               f"GitHub API saatlik sınırı dolmuş olabilir") from error

    def latest_release(self) -> str:
        """Sunucudaki en son sürümü döndürür (ağ bağlantısı gerekir)."""
        return self._release_info(None)["tag_name"]

    # --- update ----------------------------------------------------------------

    def update(self, force: bool = False) -> str:
        """Sürümün dosyalarını sys/ içine indirir.

        Nesne bir sürümle oluşturulduysa o sürüm (hazırsa ağa çıkılmaz),
        oluşturulmadıysa sunucudaki en son sürüm alınır. Her dosya GitHub'ın
        SHA-256 özetiyle doğrulanır; doğrulanmış dosyalar yeniden indirilmez.

        Args:
            force: True ise dosyalar yeniden indirilir.

        Returns:
            Hazırlanan sürüm (ör. ``"v2026-09-01"``).
        """
        if self.__requested and self.is_ready() and not force:
            self._log(f"HPO {self.release} zaten hazır: {self.release_dir}")
            return self.release
        info = self._release_info(self.__requested)
        self.__release = info["tag_name"]
        self.__data = None
        if self.is_ready() and not force:
            self._log(f"HPO {self.release} zaten hazır: {self.release_dir}")
            return self.release

        assets = {asset["name"]: asset for asset in info["assets"]}
        missing = [name for name in _FILES if name not in assets]
        if missing:
            raise RuntimeError(f"{self.release}: sürümde bulunamayan dosyalar: {', '.join(missing)}")
        (self.release_dir / "manifest.json").unlink(missing_ok=True)
        raw = self.release_dir / "raw"
        raw.mkdir(parents=True, exist_ok=True)
        digests = {}
        for name in _FILES:
            digests[name] = self._download(assets[name], raw / name, force)

        manifest = {"release": self.release, "source": info.get("html_url", ""),
                    "published": info.get("published_at", ""), "sha256": digests,
                    "created": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        (self.release_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        self._log(f"HPO {self.release} hazır: {self.release_dir}")
        return self.release

    def _download(self, asset: dict, target: Path, force: bool) -> str:
        expected = (asset.get("digest") or "").removeprefix("sha256:")
        if target.exists() and not force and expected:
            if self._sha256(target) == expected:
                self._log(f"Var, atlandı: {target.name}")
                return expected
            self._log(f"SHA-256 tutmuyor, yeniden indiriliyor: {target.name}")

        part = target.with_name(target.name + ".part")
        digest = hashlib.sha256()
        with urllib.request.urlopen(asset["browser_download_url"], timeout=60) as response, \
                open(part, "wb") as out, \
                self._bar(asset.get("size"), f"İndiriliyor {target.name}") as bar:
            while chunk := response.read(1 << 20):
                out.write(chunk)
                digest.update(chunk)
                bar.update(len(chunk))
        if expected and digest.hexdigest() != expected:
            part.unlink(missing_ok=True)
            raise RuntimeError(f"{target.name}: SHA-256 doğrulaması başarısız")
        part.replace(target)
        return digest.hexdigest()

    def _sha256(self, path: Path) -> str:
        digest = hashlib.sha256()
        with open(path, "rb") as f, self._bar(path.stat().st_size, f"SHA-256 {path.name}") as bar:
            while chunk := f.read(1 << 20):
                digest.update(chunk)
                bar.update(len(chunk))
        return digest.hexdigest()

    def _ensure(self) -> None:
        """Veri hazır değilse indirir; sys/ içinde hazır sürüm varsa ağa çıkılmaz."""
        if not self.is_ready():
            self.update()

    # --- Yükleme ---------------------------------------------------------------

    def _rows(self, name: str) -> csv.DictReader:
        path = self.release_dir / "raw" / name
        with open(path, newline="", encoding="utf-8") as f:
            lines = [line for line in f if not line.startswith("#")]
        return csv.DictReader(lines, delimiter="\t", quoting=csv.QUOTE_NONE)

    def _parse_ontology(self) -> dict:
        terms: Dict[str, Dict[str, Any]] = {}
        alternatives: Dict[str, str] = {}
        current: Optional[Dict[str, Any]] = None
        with open(self.release_dir / "raw" / _OBO, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("["):
                    current = None
                    if line == "[Term]":
                        current = {"id": "", "name": "", "definition": "", "synonyms": [],
                                   "layperson": [], "parents": [], "obsolete": False,
                                   "replaced_by": ""}
                    continue
                if current is None or ": " not in line:
                    continue
                key, value = line.split(": ", 1)
                if key == "id":
                    current["id"] = value.split()[0]
                    terms[current["id"]] = current
                elif key == "name":
                    current["name"] = value.removeprefix("obsolete ")
                elif key == "def":
                    current["definition"] = _DEFINITION.match(value).group(1)
                elif key == "synonym":
                    text, _ = _SYNONYM.match(value).groups()
                    current["synonyms"].append(text)
                    if " layperson " in value:
                        current["layperson"].append(text)
                elif key == "is_a":
                    current["parents"].append(value.split()[0])
                elif key == "alt_id":
                    alternatives[value.split()[0]] = current["id"]
                elif key == "is_obsolete":
                    current["obsolete"] = value == "true"
                elif key == "replaced_by":
                    current["replaced_by"] = value.split()[0]

        children: Dict[str, List[str]] = defaultdict(list)
        names: Dict[str, Set[str]] = defaultdict(set)
        synonyms: Dict[str, Set[str]] = defaultdict(set)
        for term in terms.values():
            if term["obsolete"]:
                continue
            names[term["name"].casefold()].add(term["id"])
            for synonym in term["synonyms"]:
                synonyms[synonym.casefold()].add(term["id"])
            for parent in term["parents"]:
                children[parent].append(term["id"])
        return {"terms": terms, "alternatives": alternatives, "children": dict(children),
                "names": dict(names), "synonyms": dict(synonyms)}

    def _load(self) -> dict:
        if self.__data is not None:
            return self.__data
        self._ensure()
        data = self._parse_ontology()

        # hastalık -> ad, kalıtım, başlangıç; (hastalık, terim) -> kanıtlar
        diseases: Dict[str, Dict[str, Any]] = {}
        evidence: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
        negated: Set[Tuple[str, str]] = set()
        onset_terms = self._closure([ONSET_ROOT], data["children"])
        for row in self._rows(_HPOA):
            disease_id, hpo_id = row["database_id"], row["hpo_id"]
            disease = diseases.setdefault(disease_id, {"name": row["disease_name"],
                                                       "inheritance": set(), "onset": set()})
            if row["qualifier"] == "NOT":
                negated.add((disease_id, hpo_id))
                continue
            evidence[disease_id, hpo_id].add(row["evidence"])
            if row["aspect"] == "I":
                disease["inheritance"].add(hpo_id)
            elif hpo_id in onset_terms:
                disease["onset"].add(hpo_id)
            if row["onset"]:
                disease["onset"].add(row["onset"])

        association = {}
        for row in self._rows(_GENE_DISEASE):
            gene_id = row["ncbi_gene_id"].removeprefix("NCBIGene:")
            association[gene_id, row["disease_id"]] = row["association_type"]

        annotations: List[Dict[str, Any]] = []
        removed_not = removed_zero = 0
        for row in self._rows(_GENE_PHENOTYPE):
            disease_id, hpo_id = row["disease_id"], row["hpo_id"]
            value = frequency_value(row["frequency"])
            if (disease_id, hpo_id) in negated and (disease_id, hpo_id) not in evidence:
                removed_not += 1
                continue
            if value == 0:
                removed_zero += 1
                continue
            symbol = row["gene_symbol"]
            annotations.append({
                "ncbi_gene_id": row["ncbi_gene_id"],
                "gene_symbol": "" if symbol == _NO_SYMBOL else symbol,
                "hpo_id": hpo_id,
                "hpo_name": row["hpo_name"],
                "frequency": "" if row["frequency"] == "-" else row["frequency"],
                "frequency_value": value,
                "disease_id": disease_id,
                "disease_name": diseases.get(disease_id, {}).get("name", ""),
                "association_type": association.get((row["ncbi_gene_id"], disease_id), "UNKNOWN"),
                "evidence": "|".join(sorted(evidence.get((disease_id, hpo_id), ()))),
            })

        by_term: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        by_gene: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        symbols: Dict[str, str] = {}
        for row in annotations:
            by_term[row["hpo_id"]].append(row)
            by_gene[row["ncbi_gene_id"]].append(row)
            if row["gene_symbol"]:
                symbols[row["gene_symbol"].casefold()] = row["ncbi_gene_id"]
        data.update(diseases=diseases, annotations=annotations, by_term=dict(by_term),
                    by_gene=dict(by_gene), symbols=symbols, removed_not=removed_not,
                    removed_zero=removed_zero)
        self.__data = data
        return data

    # --- Ontoloji --------------------------------------------------------------

    @staticmethod
    def _closure(seeds: Iterable[str], links: Dict[str, List[str]]) -> Set[str]:
        seen: Set[str] = set()
        stack = list(seeds)
        while stack:
            current = stack.pop()
            if current not in seen:
                seen.add(current)
                stack.extend(links.get(current, ()))
        return seen

    def _resolve(self, term: str) -> str:
        """Kimlik, eski kimlik, terim adı ya da eş anlamlıdan güncel HP kimliği."""
        data = self._load()
        identifier = str(term).strip()
        identifier = data["alternatives"].get(identifier, identifier)
        record = data["terms"].get(identifier)
        if record is not None:
            if record["obsolete"]:
                if record["replaced_by"]:
                    return record["replaced_by"]
                raise ValueError(f"{identifier} ({record['name']}) kullanımdan kalkmış "
                                 f"ve yerine geçen terim yok")
            return identifier
        key = identifier.casefold()
        matches = data["names"].get(key) or data["synonyms"].get(key)
        if matches and len(matches) == 1:
            return next(iter(matches))
        if matches:
            options = ", ".join(f"{i} ({data['terms'][i]['name']})" for i in sorted(matches))
            raise ValueError(f"'{term}' birden çok terime uyuyor: {options}")
        close = [f"{row['id']} {row['name']}" for row in self.search(identifier, limit=5)]
        hint = f"; yakın terimler: {', '.join(close)}" if close else ""
        raise ValueError(f"HPO terimi bulunamadı: '{term}'{hint}")

    def _terms(self, query: TermQuery, descendants: bool = True) -> Set[str]:
        ids = [self._resolve(item) for item in _as_list(query)]
        if not descendants:
            return set(ids)
        data = self._load()
        return self._closure(ids, data["children"])

    def term(self, term: str) -> Dict[str, Any]:
        """Terimin kaydı: ``id, name, definition, synonyms, layperson, parents, children``.

        Args:
            term: HP kimliği, eski kimlik, terim adı ya da eş anlamlısı.

        Raises:
            ValueError: Terim bulunamazsa, belirsizse ya da yerine geçeni yoksa.
        """
        data = self._load()
        record = data["terms"][self._resolve(term)]
        return {"id": record["id"], "name": record["name"], "definition": record["definition"],
                "synonyms": list(record["synonyms"]), "layperson": list(record["layperson"]),
                "parents": list(record["parents"]),
                "children": list(data["children"].get(record["id"], ()))}

    def ancestors(self, term: str) -> List[str]:
        """Terimin tüm üst terimleri (kendisi hariç)."""
        data = self._load()
        identifier = self._resolve(term)
        parents = {key: value["parents"] for key, value in data["terms"].items()}
        return sorted(self._closure([identifier], parents) - {identifier})

    def descendants(self, term: str) -> List[str]:
        """Terimin tüm alt terimleri (kendisi hariç)."""
        identifier = self._resolve(term)
        return sorted(self._terms(identifier) - {identifier})

    def search(self, text: str | re.Pattern, limit: Optional[int] = None) -> List[Dict[str, Any]]:
        """Adında ya da eş anlamlılarında metin geçen terimler.

        Satır: ``id, name, matched, genes``. ``matched`` eşleşen ad ya da eş
        anlamlıdır, ``genes`` alt terimler dahil gen sayısıdır. Tam eşleşmeler
        önce, sonra gen sayısına göre sıralanır.

        Args:
            text: Büyük/küçük harf duyarsız metin ya da ``re.compile`` deseni.
            limit: En çok bu kadar satır.
        """
        data = self._load()
        pattern = text if isinstance(text, re.Pattern) else re.compile(re.escape(text), re.I)
        found = []
        for term in data["terms"].values():
            if term["obsolete"]:
                continue
            matched = next((n for n in (term["name"], *term["synonyms"]) if pattern.search(n)), None)
            if matched is not None:
                exact = not isinstance(text, re.Pattern) and matched.casefold() == text.casefold()
                found.append((not exact, term, matched))
        rows = []
        for rank, term, matched in found:
            genes = {row["ncbi_gene_id"] for t in self._terms(term["id"])
                     for row in data["by_term"].get(t, ())}
            rows.append((rank, -len(genes), term["name"], {
                "id": term["id"], "name": term["name"], "matched": matched, "genes": len(genes)}))
        rows.sort(key=lambda item: item[:3])
        return [item[3] for item in rows[:limit]]

    # --- Sorgular --------------------------------------------------------------

    def _gene_ids(self, query: Filter) -> Set[str]:
        data = self._load()
        ids = set()
        for gene in _as_list(query):
            text = str(gene)
            ids.add(text if text in data["by_gene"] else data["symbols"].get(text.casefold(), ""))
        return ids - {""}

    def _disease_ids(self, query: Filter) -> Set[str]:
        diseases = self._load()["diseases"]
        ids: Set[str] = set()
        for item in _as_list(query):
            if isinstance(item, re.Pattern):
                ids |= {key for key, value in diseases.items() if item.search(value["name"])}
            elif item in diseases:
                ids.add(item)
            else:
                wanted = str(item).casefold()
                ids |= {key for key, value in diseases.items() if value["name"].casefold() == wanted}
        if not ids:
            raise ValueError(f"Hastalık bulunamadı: {query}; aramak için diseases() kullanın")
        return ids

    @staticmethod
    def _check_choices(name: str, values: List[str], allowed: Seq[str]) -> None:
        unknown = [value for value in values if value not in allowed]
        if unknown:
            raise ValueError(f"Bilinmeyen {name}: {', '.join(map(str, unknown))}; "
                             f"seçenekler: {', '.join(allowed)}")

    def annotations(self, term: Optional[TermQuery] = None, *, exclude: Optional[TermQuery] = None,
                    gene: Optional[Filter] = None, disease: Optional[Filter] = None,
                    min_frequency: Optional[float] = None,
                    association: Optional[str | Iterable[str]] = None,
                    evidence: Optional[str | Iterable[str]] = None,
                    onset: Optional[TermQuery] = None, inheritance: Optional[TermQuery] = None,
                    descendants: bool = True) -> List[Dict[str, Any]]:
        """Filtrelere uyan gen-terim-hastalık satırları.

        Satır: ``ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency,
        frequency_value, disease_id, disease_name, association_type, evidence``.
        Filtreler modül açıklamasındadır; hepsi birlikte sağlanmalıdır (VE).

        Raises:
            ValueError: Terim bulunamazsa ya da filtre değeri geçersizse.
        """
        data = self._load()
        associations, evidences = _as_list(association), set(_as_list(evidence))
        self._check_choices("association", associations, ASSOCIATIONS)
        self._check_choices("evidence", list(evidences), EVIDENCES)
        if min_frequency is not None and not 0 <= min_frequency <= 1:
            raise ValueError(f"min_frequency 0 ile 1 arasında olmalı, {min_frequency} verildi")

        terms = None if term is None else self._terms(term, descendants)
        if gene is not None:
            rows: Iterable = [row for gene_id in self._gene_ids(gene)
                              for row in data["by_gene"][gene_id]]
        elif terms is not None:
            rows = [row for hpo_id in terms for row in data["by_term"].get(hpo_id, ())]
        else:
            rows = data["annotations"]
        disease_ids = None if disease is None else self._disease_ids(disease)
        onsets = None if onset is None else self._terms(onset)
        modes = None if inheritance is None else self._terms(inheritance)
        excluded = set()
        if exclude is not None:
            excluded = {row["ncbi_gene_id"] for hpo_id in self._terms(exclude)
                        for row in data["by_term"].get(hpo_id, ())}

        result = []
        for row in rows:
            if row["ncbi_gene_id"] in excluded:
                continue
            if terms is not None and row["hpo_id"] not in terms:
                continue
            if disease_ids is not None and row["disease_id"] not in disease_ids:
                continue
            if min_frequency is not None and (row["frequency_value"] is None
                                              or row["frequency_value"] < min_frequency):
                continue
            if associations and row["association_type"] not in associations:
                continue
            if evidences and not evidences & set(row["evidence"].split("|")):
                continue
            if onsets is not None or modes is not None:
                info = data["diseases"].get(row["disease_id"], {})
                if onsets is not None and not onsets & info.get("onset", set()):
                    continue
                if modes is not None and not modes & info.get("inheritance", set()):
                    continue
            result.append(dict(row))
        return result

    def genes(self, term: Optional[TermQuery] = None, **filters: Any) -> List[Dict[str, Any]]:
        """Filtrelere uyan genlerin özeti; en çok hastalıkla ilişkili olan önce.

        Satır: ``ncbi_gene_id, gene_symbol, diseases, terms, max_frequency``.
        ``diseases`` ve ``terms`` filtreden geçen satırlardaki sayılardır.
        """
        grouped: Dict[str, Dict[str, Any]] = {}
        for row in self.annotations(term, **filters):
            entry = grouped.setdefault(row["ncbi_gene_id"], {
                "ncbi_gene_id": row["ncbi_gene_id"], "gene_symbol": row["gene_symbol"],
                "_diseases": set(), "_terms": set(), "max_frequency": None})
            entry["_diseases"].add(row["disease_id"])
            entry["_terms"].add(row["hpo_id"])
            value = row["frequency_value"]
            if value is not None and (entry["max_frequency"] is None or value > entry["max_frequency"]):
                entry["max_frequency"] = value
        rows = []
        for entry in grouped.values():
            rows.append({"ncbi_gene_id": entry["ncbi_gene_id"], "gene_symbol": entry["gene_symbol"],
                         "diseases": len(entry["_diseases"]), "terms": len(entry["_terms"]),
                         "max_frequency": entry["max_frequency"]})
        return sorted(rows, key=lambda r: (-r["diseases"], r["gene_symbol"]))

    def symbols(self, term: Optional[TermQuery] = None, **filters: Any) -> List[str]:
        """Filtrelere uyan genlerin sembolleri, alfabetik (sembolü olmayanlar hariç)."""
        return sorted({row["gene_symbol"] for row in self.genes(term, **filters) if row["gene_symbol"]})

    def phenotypes(self, gene: Filter, **filters: Any) -> List[Dict[str, Any]]:
        """Genin fenotip terimleri; kalıtım ve başlangıç gibi terimler dahil değildir.

        Satır: ``hpo_id, hpo_name, diseases, max_frequency``; en çok hastalıkta
        geçen önce.

        Raises:
            KeyError: Gen HPO'da yoksa.
        """
        if not self._gene_ids(gene):
            raise KeyError(f"HPO {self.release}: '{gene}' geni bulunamadı")
        grouped: Dict[str, Dict[str, Any]] = {}
        for row in self.annotations(filters.pop("term", PHENOTYPE_ROOT), gene=gene, **filters):
            entry = grouped.setdefault(row["hpo_id"], {"hpo_id": row["hpo_id"], "hpo_name": row["hpo_name"],
                                                       "_diseases": set(), "max_frequency": None})
            entry["_diseases"].add(row["disease_id"])
            value = row["frequency_value"]
            if value is not None and (entry["max_frequency"] is None or value > entry["max_frequency"]):
                entry["max_frequency"] = value
        rows = [{"hpo_id": e["hpo_id"], "hpo_name": e["hpo_name"], "diseases": len(e["_diseases"]),
                 "max_frequency": e["max_frequency"]} for e in grouped.values()]
        return sorted(rows, key=lambda r: (-r["diseases"], r["hpo_name"]))

    def diseases(self, name: Optional[str | re.Pattern] = None, **filters: Any) -> List[Dict[str, Any]]:
        """Hastalıklar: adında ``name`` geçenler ve/veya filtrelere uyanlar.

        Satır: ``disease_id, disease_name, genes, terms, inheritance, onset``;
        ``inheritance`` ve ``onset`` terim adlarıdır (``|`` ile ayrılmış). Gen
        sayısı filtreden geçen satırlardan hesaplanır; en çok gen önce.

        Args:
            name: Büyük/küçük harf duyarsız metin ya da ``re.compile`` deseni.
            **filters: ``annotations`` filtreleri, ör. ``term="Ataxia"``.
        """
        data = self._load()
        pattern = name if isinstance(name, re.Pattern) else re.compile(re.escape(name or ""), re.I)
        wanted = {key for key, value in data["diseases"].items() if pattern.search(value["name"])}
        genes: Dict[str, Set[str]] = defaultdict(set)
        terms: Dict[str, Set[str]] = defaultdict(set)
        for row in self.annotations(**filters):
            genes[row["disease_id"]].add(row["ncbi_gene_id"])
            terms[row["disease_id"]].add(row["hpo_id"])
        if filters:
            wanted &= set(genes)

        def names(ids: Iterable[str]) -> str:
            return "|".join(sorted(data["terms"][i]["name"] for i in ids if i in data["terms"]))

        rows = [{"disease_id": key, "disease_name": data["diseases"][key]["name"],
                 "genes": len(genes[key]), "terms": len(terms[key]),
                 "inheritance": names(data["diseases"][key]["inheritance"]),
                 "onset": names(data["diseases"][key]["onset"])} for key in wanted]
        return sorted(rows, key=lambda r: (-r["genes"], r["disease_name"]))

    def summary(self) -> Dict[str, Any]:
        """Sürüm, kayıt sayıları ve yükleme sırasında çıkarılan satırlar."""
        data = self._load()
        return {"release": self.release,
                "terms": sum(not t["obsolete"] for t in data["terms"].values()),
                "genes": len(data["by_gene"]), "diseases": len(data["diseases"]),
                "annotations": len(data["annotations"]),
                "removed_not": data["removed_not"], "removed_zero_frequency": data["removed_zero"]}

    # --- crn nesneleri ---------------------------------------------------------

    def _write_features(self, obj: Gene | Sequence, rows: List[Dict[str, Any]]) -> None:
        """Satırlardaki terim ve hastalıkları obj.features'a tekrarsız yazar."""
        for key in _FEATURE_KEYS:
            obj.features.remove(key)
            for value in dict.fromkeys(row[key] for row in rows if row[key]):
                obj.features.add(key, value)
        obj.features.set("hpo_release", self.release)

    def annotate(self, obj: Gene | Sequence, gene: Optional[str] = None,
                 term: Optional[TermQuery] = PHENOTYPE_ROOT, **filters: Any) -> Gene | Sequence:
        """``Gene`` ya da diziye HPO terimlerini ve hastalıklarını yazar.

        features'a ``hpo_id``, ``hpo_name``, ``disease_id``, ``disease_name``
        listeleri ve ``hpo_release`` eklenir; aynı anahtarların eski değerleri
        silinir. Gen HPO'da yoksa bu listeler boş kalır. Yazılanlar CRN
        dosyalarında korunur.

        Args:
            obj: ``Gene``, ``DNA``, ``RNA`` ya da ``PRT``.
            gene: Gen sembolü; verilmezse features'taki ``gene_name``, o da yoksa
                ``Gene`` nesnesinin adı kullanılır.
            term: Yazılacak terimleri sınırlar; varsayılan tüm fenotipler.
            **filters: ``annotations`` filtreleri.

        Returns:
            Aynı nesne.

        Raises:
            ValueError: Gen sembolü belirlenemezse.
        """
        symbol = gene or obj.features.first("gene_name") or (obj.name if isinstance(obj, Gene) else None)
        if not symbol:
            raise ValueError("Gen sembolü belirlenemedi; gene= ile verin")
        self._write_features(obj, self.annotations(term, gene=symbol, **filters))
        return obj

    def _gencode(self) -> Gencode:
        if self.__gencode is None:
            self.__gencode = Gencode(sys_dir=self.sys_dir, verbose=False)
        return self.__gencode

    def _canonical(self, source: Gencode, gene_id: str, symbol: str) -> Gene:
        """Sembolle ana transkript; aynı adlı GENCODE kayıtları ayrılarak.

        Ad birden çok kayda uyarsa tek protein kodlayan kayıt alınır (CAST). İkisi
        de kodluyorsa (SHOX, PDE8B) HPO'nun NCBI kimliği HGNC'de Ensembl kimliğine
        çevrilir; HGNC yerelde yoksa hata aynen geçer.
        """
        try:
            return source.canonical(symbol)
        except AmbiguousGeneError:
            coding = source.genes(gene_type="protein_coding", name=f"^{re.escape(symbol)}$")
            if len(coding) == 1:
                return source.canonical(coding[0]["gene_id"])
            hgnc = Hgnc(sys_dir=self.sys_dir, verbose=False)
            if not hgnc.is_ready():
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
                isoforms = source.gene(gene_id).isoforms()
                preferred = [g for g in isoforms
                             if str(g.features.first("transcript_id")).startswith("NM_")]
                isoforms = preferred or isoforms
        if not isoforms:
            raise ValueError("protein ve transkript eşleşmedi")
        return isoforms[0]

    def _build(self, wanted: List[Dict[str, Any]], source: Any, group: Optional[str],
               key: str, filters: Dict[str, Any], desc: str
               ) -> Tuple[List[Gene], Dict[str, str]]:
        source = source if source is not None else self._gencode()
        genes, skipped = [], {}
        for row in self._progress(wanted, desc, "gen"):
            label = row["gene_symbol"] or row["ncbi_gene_id"]
            try:
                gene = self._isoform(source, row["ncbi_gene_id"], row["gene_symbol"])
            except (KeyError, ValueError, RuntimeError) as error:
                skipped[label] = str(error).strip("'\"")
                continue
            self._write_features(gene, self.annotations(**{**filters, "gene": row["ncbi_gene_id"]}))
            if group is not None:
                gene.features.set(key, group)
            genes.append(gene)
        return genes, skipped

    def to_genes(self, term: Optional[TermQuery] = None, source: Any = None,
                 group: Optional[str] = None, key: str = "group", **filters: Any) -> List[Gene]:
        """Filtrelere uyan genleri ``Gene`` nesneleri olarak çıkarır.

        Her gen bir DNA, bir transkript ve bir protein içerir, yani doğrudan
        ``CUB``'a verilebilir. Eşleşen HPO terimleri ve hastalıklar features'a
        yazılır. Protein kodlamayan ya da standart kodla çevrilemeyen genler
        (mitokondri gibi) atlanır ve ``UserWarning`` ile bildirilir.

        Args:
            term: Terim ya da terimler.
            source: ``Gencode`` (MANE_Select, yoksa Ensembl_canonical; ORF'u
                proteine eşleşmezse ilk uyan transkript, ``canonical_tag`` değeri
                ``other_transcript``) ya da ``Ncbi`` (NM_ öncelikli ilk izoform).
                None ise sys/ içindeki GENCODE.
            group: Verilirse her genin features'ına ``key`` anahtarıyla yazılır.
            key: Grup anahtarı.
            **filters: ``annotations`` filtreleri.
        """
        filters = {"term": term, **filters}
        genes, skipped = self._build(self.genes(**filters), source, group, key, filters,
                                     "Genler çıkarılıyor")
        if skipped:
            warnings.warn(f"{len(skipped)} gen atlandı: {', '.join(sorted(skipped)[:10])}"
                          f"{' ...' if len(skipped) > 10 else ''}", UserWarning, stacklevel=2)
        return genes

    def to_sequences(self, kind: str, term: Optional[TermQuery] = None, source: Any = None,
                     **filters: Any) -> List[Sequence]:
        """Filtrelere uyan genlerin dizileri; ``write_fasta`` ile yazılabilir.

        Args:
            kind: ``"dna"`` (genomik), ``"rna"`` (ana transkript) ya da ``"protein"``.
            term, source, **filters: ``to_genes`` ile aynı.

        Returns:
            ``DNA``, ``RNA`` ya da ``PRT`` listesi. Her dizinin features'ında
            gen adı, HPO terimleri ve hastalıklar vardır.

        Raises:
            ValueError: kind geçersizse.
        """
        self._check_choices("kind", [kind], SEQUENCE_KINDS)
        sequences = []
        for gene in self.to_genes(term, source, **filters):
            if kind == "dna":
                item = gene.dna
            elif kind == "rna":
                item = next(iter(gene.transcripts.values()))
            else:
                item = next(iter(gene.proteins.values()))
            item.features.set("gene_name", gene.features.first("gene_name"))
            for name in _FEATURE_KEYS:
                item.features.remove(name)
                for value in gene.features.get(name, []):
                    item.features.add(name, value)
            item.features.set("hpo_release", self.release)
            sequences.append(item)
        return sequences

    def pack(self, groups: Dict[str, Any], file_name: str | Path, source: Any = None,
             key: str = "group", overlap: str = "exclude", **filters: Any) -> Dict[str, Any]:
        """Gen gruplarını tek CRN listesine yazar.

        Dosya ``DifferentialCodonUsage.from_crn`` ile doğrudan açılır; genin
        grubu features'ta ``key`` anahtarındadır.

        Args:
            groups: grup adı -> terim, terim listesi ya da gruba özel filtreler
                (ör. ``{"term": "Ataxia", "min_frequency": 0.3}``).
            file_name: Yazılacak CRN dosyası.
            source: ``to_genes`` ile aynı.
            key: Grup anahtarı.
            overlap: Birden çok gruba giren gen için ``"exclude"`` (hepsinden
                çıkar, varsayılan), ``"first"`` (ilk grupta kalır) ya da ``"error"``.
            **filters: Tüm gruplara uygulanan ``annotations`` filtreleri.

        Returns:
            ``file, release, groups`` (grup -> gen sayısı), ``overlap`` (çakışan
            semboller) ve ``skipped`` (atlanan gen -> neden).

        Raises:
            ValueError: overlap geçersizse, gruplar çakışıyorsa (``"error"``) ya
                da bir grup boş kalırsa.
        """
        self._check_choices("overlap", [overlap], ("exclude", "first", "error"))
        if not groups:
            raise ValueError("En az bir grup gerekli")
        wanted: Dict[str, List[Dict[str, Any]]] = {}
        settings: Dict[str, Dict[str, Any]] = {}
        owners: Dict[str, List[str]] = defaultdict(list)
        labels: Dict[str, str] = {}
        for name, spec in groups.items():
            settings[name] = {**filters, **(spec if isinstance(spec, dict) else {"term": spec})}
            wanted[name] = self.genes(**settings[name])
            for row in wanted[name]:
                owners[row["ncbi_gene_id"]].append(name)
                labels[row["ncbi_gene_id"]] = row["gene_symbol"] or row["ncbi_gene_id"]

        shared = {gene_id: names for gene_id, names in owners.items() if len(names) > 1}
        if shared and overlap == "error":
            raise ValueError(f"{len(shared)} gen birden çok grupta: "
                             f"{', '.join(sorted(labels[g] for g in shared)[:10])}")
        for name in wanted:
            wanted[name] = [row for row in wanted[name] if row["ncbi_gene_id"] not in shared
                            or (overlap == "first" and shared[row["ncbi_gene_id"]][0] == name)]
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

    def __repr__(self) -> str:
        return f"Hpo({self.release!r}, {str(self.hpo_dir)!r})"
