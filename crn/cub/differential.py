"""Kodon kullanım yanlılığı (codon usage bias) — J. Diferansiyel kodon kullanımı.

İki gen grubunu kodon kodon karşılaştırır: bir kodon, karşılaştırma havuzu
içinde bir grupta ötekine göre anlamlı biçimde daha çok mu seçiliyor? Her kodon
için şu 2×2 tablo kurulur::

                 kodon      geri kalan
    Grup A         a            b
    Grup B         c            d

"Geri kalan" karşılaştırma düzeyine bağlıdır:

==========  ==================================  =================================
 Düzey      Geri kalan                          Ne ölçer
==========  ==================================  =================================
 global     tüm öteki anlamlı kodonlar          genel bolluk; amino asit bileşimi
                                                farkıyla karışır
 family     aynı amino asidin öteki kodonları   saf sinonim tercih (esas düzey)
 subfamily  aynı amino asit ve ilk iki baz      wobble konumundaki tercih
==========  ==================================  =================================

Stop kodonları hiçbir düzeye girmez; tek kodonlu amino asitler (standart kodda
Met, Trp) aile ve alt aile düzeylerinde test edilmez.

Her satırda şunlar raporlanır:

* **Olasılık oranı:** Haldane-Anscombe düzeltmesiyle
  ``OR = ((a+0.5)(d+0.5)) / ((b+0.5)(c+0.5))``. ``log2_or`` pozitifse kodon A
  grubunda, negatifse B grubunda göreli daha çok seçilir.
* **p:** iki yönlü Fisher kesin testi (``scipy.stats.fisher_exact``). Olasılığı
  gözlenen tablonunkinden büyük olmayan tüm tabloların hipergeometrik
  olasılıkları toplanır.
* **q:** o düzeydeki testler üzerinde Benjamini-Hochberg düzeltmesi,
  ``q_(i) = min_{k >= i} (m / k) p_(k)``.

Genom ölçeğinde p değerleri çift duyarlıklı sayının altına inebilir.
``fisher_exact`` yaklaşık 1e-295'in altında kuyruklardan birini kaybeder, daha
aşağıda 0 döndürür. Bu yüzden p ``_SCIPY_MIN_P`` altına indiğinde aynı toplam
log uzayında yeniden hesaplanır: ``p`` ve ``q`` 0'a yuvarlansa bile ``log10_p``
ve ``log10_q`` sonludur; yanardağ grafiği bunları kullanır.

Fisher testi sayımları grup içinde havuzlar ve her kodonu bağımsız bir gözlem
sayar. Oysa genler kodon kullanımında birbirinden zaten farklıdır, bu yüzden
gerçek fark olmayan rastgele gen gruplarında bile onlarca kodon anlamlı çıkabilir.
``permutation`` bunu sınar: genlerin grup etiketlerini karıştırır ve log2 OR'un
gerçek değerini rastgele gruplarla karşılaştırır. Fisher'da da permütasyonda da
geçen kodonlar ``robust`` olarak işaretlenir; ``interpret`` yalnızca bunları yorumlar.

Permütasyon genler arası çeşitliliği hesaba katar, ama iki grubun uzunluk ya da
bölgesel GC farkını gidermez. Kontrol grubu bu özelliklere göre seçilmelidir.

Gruplar tek bir CRN dosyasındaki gen listesinden okunur. Her genin grubu
features'ta tutulur::

    from crn import write_crn_list
    from crn.cub import DifferentialCodonUsage

    for gene in yuksek:
        gene.features.set("group", "yuksek")
    for gene in dusuk:
        gene.features.set("group", "dusuk")
    write_crn_list(yuksek + dusuk, "gruplar.crn")

    dcu = DifferentialCodonUsage.from_crn("gruplar.crn", "yuksek", "dusuk")
    dcu.table("family")                 # tüm kodonlar
    dcu.significant("family")           # Fisher q < 0.05, etkiye göre sıralı
    dcu.robust("family")                # gen düzeyinde permütasyondan da geçenler
    dcu.compare_levels(permutations=5000)   # düzeyler arası okuma
    print(dcu.interpret())              # Türkçe yorum (permütasyona dayanır)
    dcu.plot_volcano("family", file_name="yanardag.png")

Her grup ayrı da verilebilir: ``CUB``, ``Gene`` listesi ya da
``write_crn_list`` ile yazılmış dosya yolu::

    dcu = DifferentialCodonUsage(cub_a, "grup_b.crn", name_a="hasta", name_b="kontrol")

Genler ``CUB`` kuralına uymalıdır: bir DNA, bir transkript, bir protein. Çok
izoformlu genler önce ``Gene.isoforms()`` ile ayrılır.
"""
from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Sequence, Tuple

from crn import fileio
from crn.constants import GeneticCode
from crn.cub import _plotting
from crn.cub.composition import Composition
from crn.cub.core import CUB, CubWarning
from crn.cub.enc import EffectiveCodons
from crn.cub.synonymous import check_choice, check_cub
from crn.gene import Gene

if TYPE_CHECKING:
    from matplotlib.axes import Axes

LEVELS = ("global", "family", "subfamily")
_LOG10 = math.log(10)
# fisher_exact bunun altında log uzayıyla değiştirilir. Bozulma p≈1e-295'te
# başladığı için sınır güvenli tarafta seçildi.
_SCIPY_MIN_P = 1e-200
# Yorumda kullanılan uzlaşımsal sınırlar; testin kendisini etkilemez.
_LARGE_EFFECT = 1.0       # |log2 OR| >= 1: iki kat ve üzeri
_DOMINANT_GENE = 0.25     # bir genin grubun kodonlarındaki payı
_MIN_THIRD_BASE = 3       # 3. konum karşılaştırması için grup başına en az anlamlı kodon
_GC3S_GAP = 0.05          # GC3s farkı bu değerin üstündeyse bileşim açıklaması güçlenir
PERMUTATIONS = 5000       # gen etiketi karıştırma sayısı


def _log_comb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _log_add(x: float, y: float) -> float:
    """log(e^x + e^y); çok küçük değerlerde taşmaz."""
    if x < y:
        x, y = y, x
    return x + math.log1p(math.exp(y - x))


def fisher_log10_p(a: int, b: int, c: int, d: int) -> float:
    """[[a, b], [c, d]] tablosunun iki yönlü Fisher kesin testi; log10(p).

    p ``scipy.stats.fisher_exact`` ile hesaplanır; ``_SCIPY_MIN_P`` altındaysa
    sonuç log uzayında yeniden bulunur, böylece değer her zaman sonludur.
    scipy ilk çağrıda yüklenir (``crn.cub`` importunu yavaşlatmaz).
    """
    from scipy.stats import fisher_exact

    p = float(fisher_exact([[a, b], [c, d]]).pvalue)
    if p >= _SCIPY_MIN_P:
        return min(0.0, math.log10(p))
    return _fisher_log10_tails(a, b, c, d)


def _fisher_log10_tails(a: int, b: int, c: int, d: int) -> float:
    """İki yönlü Fisher log10(p) değeri, tümüyle log uzayında.

    Kuyruklar en olası tablodan uzaklaşarak ardışık oranlarla toplanır; hesap
    kuyruktaki tablo sayısına göre değil, olasılığın küçülmesine göre biter, bu
    yüzden milyonluk sayımlarda da hızlıdır.
    """
    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2
    low, high = max(0, col1 - row2), min(row1, col1)
    if low == high:
        return 0.0

    def log_p(x: int) -> float:
        return _log_comb(row1, x) + _log_comb(row2, col1 - x) - _log_comb(total, col1)

    def tail(start: int, step: int) -> float:
        """log Σ p(x), x = start'tan ``step`` yönünde uca kadar."""
        term, summed, x = 1.0, 1.0, start
        while low <= x + step <= high:
            if step > 0:
                term *= (row1 - x) * (col1 - x) / ((x + 1) * (row2 - col1 + x + 1))
            else:
                term *= x * (row2 - col1 + x) / ((row1 - x + 1) * (col1 - x + 1))
            summed += term
            x += step
            if term < summed * 1e-17:
                break
        return log_p(start) + math.log(summed)

    mode = min(high, max(low, (row1 + 1) * (col1 + 1) // (total + 2)))
    if log_p(mode + 1 if mode < high else mode) > log_p(mode):
        mode += 1
    limit = log_p(a) + 1e-7     # R'daki göreli tolerans; eş olasılıklı tablolar dahil edilir

    def first_at_most(first: int, last: int, step: int) -> Optional[int]:
        """Mod'dan ``step`` yönünde olasılığı limit altına ilk inen x (ikili arama)."""
        if not low <= last <= high or log_p(last) > limit:
            return None
        while first != last:
            middle = (first + last) // 2 if step > 0 else (first + last + 1) // 2
            if log_p(middle) <= limit:
                last = middle
            else:
                first = middle + step
        return first

    left = first_at_most(mode, low, -1)
    right = first_at_most(mode, high, 1)
    if right is not None and right == left:     # gözlenen tablo en olası tablo: iki kuyruk bitişik
        right = right + 1 if right < high else None
    tails = []
    if left is not None:
        tails.append(tail(left, -1))
    if right is not None:
        tails.append(tail(right, 1))
    if not tails:       # hiçbir tablo eşiğin altında değil
        return 0.0
    log_sum = tails[0] if len(tails) == 1 else _log_add(*tails)
    return min(0.0, log_sum / _LOG10)


def benjamini_hochberg_log10(log10_p: Sequence[float]) -> List[float]:
    """log10(p) değerlerinden Benjamini-Hochberg log10(q) değerleri, aynı sırayla."""
    m = len(log10_p)
    order = sorted(range(m), key=log10_p.__getitem__)
    result = [0.0] * m
    running = 0.0
    for rank in range(m, 0, -1):
        index = order[rank - 1]
        running = min(running, math.log10(m / rank) + log10_p[index])
        result[index] = running
    return result


def _power10(value: float) -> float:
    return 10.0 ** value if value > -320 else 0.0


class DifferentialCodonUsage:
    """İki gen grubu arasında kodon kodon diferansiyel kullanım testi.

    Args:
        group_a, group_b: ``CUB``, ``Gene`` listesi ya da ``write_crn_list``
            ile yazılmış dosya yolu. Tek dosyadaki gruplar için ``from_crn``.
        name_a, name_b: Tablo ve yorumda kullanılacak grup adları.
        genetic_code: Grup ``CUB`` olarak verilmediğinde kullanılacak kod.

    Raises:
        TypeError: Grup bu tiplerden biri değilse.
        ValueError: Genetik kodlar farklıysa ya da bir gen iki grupta birden varsa.
    """

    def __init__(self, group_a: CUB | List[Gene] | str | Path, group_b: CUB | List[Gene] | str | Path,
                 name_a: str = "A", name_b: str = "B",
                 genetic_code: GeneticCode = GeneticCode.STANDARD):
        self.a = self._as_cub("group_a", group_a, genetic_code)
        self.b = self._as_cub("group_b", group_b, self.a.genetic_code)
        check_cub("group_b", self.b, self.a)
        shared = sorted(set(self.a.names) & set(self.b.names))
        if shared:
            raise ValueError(f"Aynı gen iki grupta birden olamaz: {', '.join(shared)}")
        self.name_a, self.name_b = name_a, name_b
        self._rows: Dict[str, List[Dict[str, object]]] = {}
        self._permuted: Dict[Tuple[str, int, int], List[Dict[str, float]]] = {}

    @staticmethod
    def _as_cub(name: str, group, genetic_code: GeneticCode) -> CUB:
        if isinstance(group, CUB):
            return group
        if isinstance(group, list):
            return CUB(group, genetic_code)
        if isinstance(group, (str, Path)):
            return CUB(fileio.read_crn_list(group), genetic_code)
        raise TypeError(f"{name} bir CUB, Gene listesi ya da CRN dosya yolu olmalı, "
                        f"{type(group).__name__} verildi")

    @classmethod
    def from_crn(cls, file_name: str | Path, group_a: str, group_b: str, key: str = "group",
                 genetic_code: GeneticCode = GeneticCode.STANDARD) -> DifferentialCodonUsage:
        """Tek CRN listesindeki genleri features'taki grup bilgisine göre ayırır.

        Args:
            file_name: ``write_crn_list`` ile yazılmış gen listesi.
            group_a, group_b: Karşılaştırılacak grupların adları.
            key: Grup adının tutulduğu features anahtarı.
            genetic_code: Kullanılacak genetik kod.

        Raises:
            TypeError: Dosyada Gene olmayan kayıt varsa.
            ValueError: Gruplardan biri boşsa (mesajda dosyadaki gruplar yazılır).
        """
        items = fileio.read_crn_list(file_name)
        others = [type(item).__name__ for item in items if not isinstance(item, Gene)]
        if others:
            raise TypeError(f"{file_name}: yalnızca Gene kayıtları olmalı; "
                            f"{len(others)} başka kayıt var ({others[0]} ...)")
        groups: Dict[object, List[Gene]] = {}
        for gene in items:
            groups.setdefault(gene.features.first(key), []).append(gene)
        for wanted in (group_a, group_b):
            if wanted not in groups:
                found = sorted(str(g) for g in groups if g is not None)
                raise ValueError(f"{file_name}: '{wanted}' grubunda gen yok "
                                 f"(features['{key}'] değerleri: {', '.join(found) or 'yok'})")
        skipped = len(items) - len(groups[group_a]) - len(groups[group_b])
        if skipped > 0:
            warnings.warn(f"{skipped} gen iki gruptan birinde olmadığı için analize girmedi",
                          CubWarning, stacklevel=2)
        return cls(groups[group_a], groups[group_b], name_a=group_a, name_b=group_b,
                   genetic_code=genetic_code)

    # --- Test -------------------------------------------------------------------

    def _pools(self, level: str) -> List[Tuple[str, Tuple[str, ...]]]:
        """(havuz adı, havuz kodonları); düzeyin test birimleri."""
        code = self.a.genetic_code
        if level == "global":
            return [("all", tuple(code.sense_codons))]
        families = [(aa, tuple(codons)) for aa, codons in code.families.items() if len(codons) > 1]
        if level == "family":
            return families
        pools = []
        for aa, codons in families:
            prefixes = {codon[:2]: [c for c in codons if c[:2] == codon[:2]] for codon in codons}
            pools += [(f"{aa}_{prefix}", tuple(members)) for prefix, members in prefixes.items()
                      if len(members) > 1]
        return pools

    def _compute(self, level: str) -> List[Dict[str, object]]:
        if level not in self._rows:
            table = self.a.genetic_code.table
            counts_a = self.a.codon_counts()
            counts_b = self.b.codon_counts()
            rows: List[Dict[str, object]] = []
            for pool, codons in self._pools(level):
                pool_a = sum(counts_a[c] for c in codons)
                pool_b = sum(counts_b[c] for c in codons)
                for codon in codons:
                    a, c = counts_a[codon], counts_b[codon]
                    b, d = pool_a - a, pool_b - c
                    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                    log10_p = fisher_log10_p(a, b, c, d)
                    rows.append({"level": level, "codon": codon, "amino_acid": table[codon],
                                 "pool": pool, "a": a, "b": b, "c": c, "d": d,
                                 "share_a": a / pool_a if pool_a else 0.0,
                                 "share_b": c / pool_b if pool_b else 0.0,
                                 "odds_ratio": odds, "log2_or": math.log2(odds),
                                 "p": _power10(log10_p), "log10_p": log10_p})
            for row, log10_q in zip(rows, benjamini_hochberg_log10([r["log10_p"] for r in rows])):
                row["q"] = _power10(log10_q)
                row["log10_q"] = log10_q
            self._rows[level] = rows
        return self._rows[level]

    def table(self, level: str = "family", q: float = 0.05,
              file_name: Optional[str | Path] = None) -> List[Dict[str, object]]:
        """Düzeydeki her kodon için bir satır.

        Satır anahtarları: ``level, codon, amino_acid, pool, a, b, c, d, share_a,
        share_b, odds_ratio, log2_or, p, log10_p, q, log10_q, significant,
        favored``. ``share_*`` kodonun havuzdaki payıdır; ``favored`` anlamlı
        satırda kodonu daha çok seçen grubun adı, değilse boştur.

        Args:
            level: ``"global"``, ``"family"`` ya da ``"subfamily"``.
            q: Anlamlılık eşiği.
            file_name: Verilirse tablo sekmeyle ayrılmış (TSV) olarak da yazılır.

        Raises:
            ValueError: level geçersizse ya da q 0-1 aralığında değilse.
        """
        check_choice("level", level, LEVELS)
        if not 0 < q <= 1:
            raise ValueError(f"q 0 ile 1 arasında olmalı, {q} verildi")
        rows = []
        for row in self._compute(level):
            significant = row["log10_q"] < math.log10(q)
            favored = "" if not significant else self.name_a if row["log2_or"] > 0 else self.name_b
            rows.append({**row, "significant": significant, "favored": favored})
        if file_name is not None:
            with open(file_name, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
        return rows

    def significant(self, level: str = "family", q: float = 0.05,
                    min_log2_or: float = 0.0) -> List[Dict[str, object]]:
        """q eşiğini geçen satırlar, etki büyüklüğüne (|log2 OR|) göre büyükten küçüğe.

        Args:
            level: Karşılaştırma düzeyi.
            q: Anlamlılık eşiği.
            min_log2_or: Satırın kalması için gereken en küçük |log2 OR|.
        """
        rows = [r for r in self.table(level, q)
                if r["significant"] and abs(r["log2_or"]) >= min_log2_or]
        return sorted(rows, key=lambda r: abs(r["log2_or"]), reverse=True)

    # --- Gen düzeyinde permütasyon ----------------------------------------------

    def _gene_counts(self, cub: CUB, codons: Sequence[str]) -> Any:
        import numpy as np
        return np.array([[cub.codon_counts(gene)[codon] for codon in codons] for gene in cub.names],
                        dtype=float)

    def _permute(self, level: str, n: int, seed: int) -> List[Dict[str, float]]:
        """Her satır için ``perm_p`` ve ``perm_q``; sonuç saklanır."""
        key = (level, n, seed)
        if key in self._permuted:
            return self._permuted[key]
        import numpy as np

        rows = self._compute(level)
        codons = [row["codon"] for row in rows]
        counts = np.vstack([self._gene_counts(self.a, codons), self._gene_counts(self.b, codons)])
        pools: Dict[str, List[int]] = {}
        for index, row in enumerate(rows):
            pools.setdefault(row["pool"], []).append(index)
        pool_counts = np.zeros_like(counts)
        for members in pools.values():
            pool_counts[:, members] = counts[:, members].sum(axis=1, keepdims=True)
        total, pool_total = counts.sum(axis=0), pool_counts.sum(axis=0)

        def log2_or(a: Any, pool_a: Any) -> Any:
            b, c = pool_a - a, total - a
            d = pool_total - pool_a - c
            return np.log2((a + 0.5) * (d + 0.5) / ((b + 0.5) * (c + 0.5)))

        size_a = len(self.a)
        observed = np.abs(log2_or(counts[:size_a].sum(axis=0), pool_counts[:size_a].sum(axis=0)))
        generator = np.random.default_rng(seed)
        extreme = np.zeros(len(rows))
        for _ in range(n):
            chosen = generator.permutation(len(counts))[:size_a]
            extreme += np.abs(log2_or(counts[chosen].sum(axis=0),
                                      pool_counts[chosen].sum(axis=0))) >= observed
        log10_p = [math.log10(value) for value in (extreme + 1) / (n + 1)]
        log10_q = benjamini_hochberg_log10(log10_p)
        result = [{"perm_p": 10.0 ** p, "perm_q": 10.0 ** adjusted}
                  for p, adjusted in zip(log10_p, log10_q)]
        self._permuted[key] = result
        return result

    def permutation(self, level: str = "family", q: float = 0.05, n: int = PERMUTATIONS,
                    seed: int = 0, file_name: Optional[str | Path] = None) -> List[Dict[str, object]]:
        """Gen düzeyinde permütasyon testi; ``table`` satırlarına üç anahtar ekler.

        Genler iki gruba, grup büyüklükleri korunarak rastgele yeniden dağıtılır
        ve her kodon için log2 OR yeniden hesaplanır. ``perm_p``, rastgele
        grupların en az gerçek fark kadar uç sonuç verme oranıdır,
        ``(uç sonuç + 1) / (n + 1)``. ``perm_q`` bunun Benjamini-Hochberg
        düzeltmesidir. ``robust`` True ise kodon hem Fisher'da hem permütasyonda
        q eşiğini geçmiştir.

        Sonuç seed için tekrar üretilebilir ve (level, n, seed) başına saklanır.

        Args:
            level: Karşılaştırma düzeyi.
            q: Anlamlılık eşiği.
            n: Karıştırma sayısı. En küçük p ``1 / (n + 1)`` olduğundan n, düzeydeki
                test sayısının ``1 / q`` katından az olmamalıdır.
            seed: Rastgele sayı üretecinin tohumu.
            file_name: Verilirse tablo TSV olarak da yazılır.

        Raises:
            ValueError: level, q ya da n geçersizse.
        """
        rows = self.table(level, q)
        if n < 1:
            raise ValueError(f"n en az 1 olmalı, {n} verildi")
        limit = math.log10(q)
        result = []
        for row, values in zip(rows, self._permute(level, n, seed)):
            result.append({**row, **values, "robust": row["significant"]
                           and math.log10(values["perm_q"]) < limit})
        if file_name is not None:
            with open(file_name, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=list(result[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(result)
        return result

    def robust(self, level: str = "family", q: float = 0.05, min_log2_or: float = 0.0,
               n: int = PERMUTATIONS, seed: int = 0) -> List[Dict[str, object]]:
        """Fisher ve permütasyonda birlikte anlamlı satırlar, |log2 OR|'a göre büyükten küçüğe.

        Args:
            level: Karşılaştırma düzeyi.
            q: Anlamlılık eşiği.
            min_log2_or: Satırın kalması için gereken en küçük |log2 OR|.
            n, seed: ``permutation`` ile aynı.
        """
        rows = [r for r in self.permutation(level, q, n, seed)
                if r["robust"] and abs(r["log2_or"]) >= min_log2_or]
        return sorted(rows, key=lambda r: abs(r["log2_or"]), reverse=True)

    def _smallest_p(self, n: int) -> float:
        """Ulaşılabilecek en küçük perm_p; grup küçükse ayrık bölünme sayısıyla sınırlıdır.

        BH'de perm_q en iyi durumda (bütün kodonlar bu p'deyse) bu değere iner; tek
        kodonun fark gösterdiği durumda ise ``test sayısı × p`` olur.
        """
        size_a, size_b = len(self.a), len(self.b)
        splits = math.comb(size_a + size_b, size_a)
        if size_a == size_b:        # iki grup yer değiştirince |log2 OR| aynı kalır
            splits //= 2
        return 1 / min(splits, n + 1)

    # --- Yorum ------------------------------------------------------------------

    def compare_levels(self, q: float = 0.05, permutations: int = 0,
                       seed: int = 0) -> List[Dict[str, object]]:
        """Eş anlamlı kodonların düzeyler arası okuması.

        Satır: ``codon, amino_acid, log2_or_global, log2_or_family,
        log2_or_subfamily`` (alt ailesi olmayan kodonda None), ``reading``.
        ``reading`` şu değerlerden biridir:

        * ``"sinonim"``: ailede anlamlı, küreselde değil; fark saf sinonim tercih.
        * ``"bilesim"``: küreselde anlamlı, ailede değil ya da ailede etkisi
          yarıdan aza iniyor; fark büyük ölçüde amino asit bileşiminden.
        * ``"karma"``: iki düzeyde de anlamlı ve etki korunuyor.
        * ``"wobble"``: yalnızca alt ailede anlamlı; tercih wobble konumunda.
        * ``""``: hiçbir düzeyde anlamlı değil.

        Args:
            q: Anlamlılık eşiği.
            permutations: 0 ise "anlamlı" Fisher'a göre okunur; büyükse o kadar
                karıştırmayla ``robust`` satırlar kullanılır. Fisher tek başına
                gen çeşitliliğini aşırı sayabildiği için ikincisi önerilir.
            seed: Permütasyon tohumu.
        """
        tables = {}
        for level in LEVELS:
            rows = (self.permutation(level, q, permutations, seed) if permutations
                    else self.table(level, q))
            if permutations:
                rows = [{**row, "significant": row["robust"]} for row in rows]
            tables[level] = {r["codon"]: r for r in rows}
        rows = []
        for codon, family in tables["family"].items():
            overall, sub = tables["global"][codon], tables["subfamily"].get(codon)
            in_family, in_global = family["significant"], overall["significant"]
            in_sub = sub is not None and sub["significant"]
            if in_family and not in_global:
                reading = "sinonim"
            elif in_global and not in_family:
                reading = "wobble" if in_sub else "bilesim"
            elif in_family and in_global:
                weakened = abs(family["log2_or"]) < abs(overall["log2_or"]) / 2
                reading = "bilesim" if weakened else "karma"
            else:
                reading = "wobble" if in_sub else ""
            rows.append({"codon": codon, "amino_acid": family["amino_acid"],
                         "log2_or_global": overall["log2_or"], "log2_or_family": family["log2_or"],
                         "log2_or_subfamily": sub["log2_or"] if sub else None, "reading": reading})
        return rows

    def _opposing_pairs(self, rows: Sequence[Dict[str, object]]) -> List[Tuple[str, str, str]]:
        """Aynı amino asitte zıt yönde anlamlı en güçlü kodon çiftleri (aa, A lehine, B lehine)."""
        pairs = []
        for aa in dict.fromkeys(r["amino_acid"] for r in rows):
            members = [r for r in rows if r["amino_acid"] == aa]
            up = [r for r in members if r["log2_or"] > 0]
            down = [r for r in members if r["log2_or"] < 0]
            if up and down:
                best_up = max(up, key=lambda r: r["log2_or"])
                best_down = min(down, key=lambda r: r["log2_or"])
                pairs.append((aa, best_up["codon"], best_down["codon"]))
        return pairs

    def _third_base(self, rows: Sequence[Dict[str, object]]) -> List[Tuple[int, int]]:
        """A ve B lehine anlamlı kodonlar için (3. konumu G/C olan, toplam) sayıları."""
        result = []
        for sign in (1, -1):
            codons = [r["codon"] for r in rows if r["log2_or"] * sign > 0]
            result.append((sum(c[2] in "GC" for c in codons), len(codons)))
        return result

    def _dominant_gene(self, cub: CUB) -> Tuple[str, float]:
        name = max(cub.names, key=cub.total_codons)
        return name, cub.total_codons(name) / cub.total_codons()

    def interpret(self, q: float = 0.05, top: int = 5, small_effect: float = 0.25,
                  min_count: int = 10, permutations: int = PERMUTATIONS, seed: int = 0) -> str:
        """Sonuçların Türkçe yorumu; yalnızca permütasyondan geçen kodonlara dayanır.

        Yön olasılık oranından, güven Fisher ile gen düzeyinde permütasyonun
        birlikte verdiği sonuçtan okunur. Fisher'ın bulup permütasyonun
        doğrulamadığı kodonlar yalnızca sayı olarak anılır, biyolojik yorum
        yapılmaz. Esas düzey ailedir; ardından düzeyler arası okuma, GC bağlamı
        ve güvenlik uyarıları gelir.

        Args:
            q: Anlamlılık eşiği.
            top: Her grup için anılacak en güçlü kodon sayısı.
            small_effect: Bu |log2 OR| değerinin altındaki anlamlı satırlar
                "etkisi küçük" sayılır (öntanımlı 0,25; OR yaklaşık 1,19).
            min_count: İki gruptaki toplam geçişi bundan az olan kodonlar
                düşük sayımlı sayılır.
            permutations: Gen etiketi karıştırma sayısı.
            seed: Permütasyon tohumu.
        """
        name_a, name_b = self.name_a, self.name_b
        rows = self.permutation("family", q, permutations, seed)
        hits = [r for r in rows if r["robust"]]
        fisher = sum(r["significant"] for r in rows)
        smallest = self._smallest_p(permutations)
        lines = [f"{name_a} ({len(self.a)} gen, {self.a.total_codons():,} kodon) ile {name_b} "
                 f"({len(self.b)} gen, {self.b.total_codons():,} kodon) karşılaştırıldı.",
                 f"Aile düzeyinde {len(rows)} eş anlamlı kodonun {fisher} tanesi Fisher testinde "
                 f"q < {q}. Gen etiketleri {permutations} kez karıştırıldığında {len(hits)} tanesi "
                 f"kaldı; yorum yalnızca bunlara dayanır."]
        if smallest >= q:
            lines.append(f"Permütasyon bu grup büyüklüğü ve karıştırma sayısıyla q < {q} eşiğine "
                         f"ulaşamaz (en küçük olası p {smallest:.3f}). Sonuç olumsuz kanıt "
                         f"sayılmamalıdır; daha çok gen gerekir.")
        elif not hits and smallest * len(rows) >= q:
            lines.append(f"Hiçbir kodon geçemedi, ama bu ayarla yalnızca çok sayıda kodonda birden "
                         f"görülen fark yakalanabilir (en küçük olası p {smallest:.3f}). Az sayıda "
                         f"kodona özgü fark kaçmış olabilir; sonuç olumsuz kanıt sayılmamalıdır.")
        elif fisher and not hits:
            lines.append("Fisher'ın bulduğu kodonların hiçbiri permütasyondan geçmedi. Bu fark "
                         "rastgele gen gruplarında da görülen doğal gen çeşitliliğiyle uyumludur; "
                         "gruba özgü sinonim tercih gösterilemedi.")
        elif not hits:
            lines.append("Bu düzeyde gruplar arasında sinonim tercih farkı gösterilemedi. Bu, "
                         "farkın olmadığını kanıtlamaz; gruplar küçükse güç yetersiz olabilir.")

        for name, sign in ((name_a, 1), (name_b, -1)):
            best = sorted((r for r in hits if r["log2_or"] * sign > 0),
                          key=lambda r: abs(r["log2_or"]), reverse=True)
            if best:
                listed = ", ".join(f"{r['codon']} ({r['amino_acid']}, OR {r['odds_ratio']:.2f})"
                                   for r in best[:top])
                lines.append(f"{name} grubunda göreli daha çok seçilen {len(best)} sağlam kodon var; "
                             f"en güçlüleri: {listed}.")

        pairs = self._opposing_pairs(hits)
        if pairs:
            listed = ", ".join(f"{aa}: {up} / {down}" for aa, up, down in pairs[:top + 3])
            lines.append(f"{len(pairs)} amino asitte tercih yön değiştiriyor. Her çiftte ilk kodon "
                         f"{name_a}, ikincisi {name_b} lehinedir ({listed}). Böyle çiftler sinonim "
                         f"tercihin yönünü doğrudan gösterir. Birinci grup ötekinden yüksek ifade "
                         f"ediliyorsa ilk kodonlar çeviri için optimal kodon adayıdır.")

        readings = self.compare_levels(q, permutations, seed)
        groups = {key: [r["codon"] for r in readings if r["reading"] == key]
                  for key in ("sinonim", "bilesim", "karma", "wobble")}
        if groups["sinonim"]:
            lines.append(f"Ailede anlamlı olup küresel düzeyde kaybolan {len(groups['sinonim'])} kodon "
                         f"var ({', '.join(groups['sinonim'][:8])}); bu, farkın saf sinonim tercihten "
                         f"geldiğini düşündürür.")
        if groups["bilesim"]:
            lines.append(f"Küreselde güçlü olup ailede zayıflayan ya da kaybolan "
                         f"{len(groups['bilesim'])} kodonun ({', '.join(groups['bilesim'][:8])}) farkının "
                         f"büyük ölçüde amino asit bileşiminden kaynaklandığı düşünülebilir.")
        if groups["karma"]:
            lines.append(f"{len(groups['karma'])} kodonda fark iki düzeyde de korunuyor; burada "
                         f"sinonim tercih ile amino asit bileşimi birlikte rol oynuyor olabilir.")
        if groups["wobble"]:
            lines.append(f"Yalnızca alt ailede beliren {len(groups['wobble'])} kodon var "
                         f"({', '.join(groups['wobble'][:8])}); bu, tercihin wobble konumuna indiğini "
                         f"gösterir.")

        gc3_a, gc3_b = Composition(self.a).gc3s(), Composition(self.b).gc3s()
        enc_a, enc_b = EffectiveCodons(self.a).enc(), EffectiveCodons(self.b).enc()
        lines.append(f"Bileşim bağlamı: GC3s {name_a} {gc3_a:.3f}, {name_b} {gc3_b:.3f}; "
                     f"ENC {name_a} {enc_a:.1f}, {name_b} {enc_b:.1f} (düşük ENC güçlü yanlılık).")
        (gc_a, total_a), (gc_b, total_b) = self._third_base(hits)
        if min(total_a, total_b) >= _MIN_THIRD_BASE:
            share_a, share_b = 100 * gc_a / total_a, 100 * gc_b / total_b
            lines.append(f"Üçüncü konumu G/C olan kodonların payı {name_a} lehine anlamlılarda "
                         f"%{share_a:.0f} ({gc_a}/{total_a}), {name_b} lehine anlamlılarda "
                         f"%{share_b:.0f} ({gc_b}/{total_b}).")
            if abs(share_a - share_b) >= 30 and abs(gc3_a - gc3_b) >= _GC3S_GAP:
                lines.append("Bu ayrışma GC3s farkıyla aynı doğrultudaysa kodon farklarının bir "
                             "bölümü GC bileşimiyle açıklanabilir; sonucu çeviri seçilimine "
                             "bağlamadan önce ENC ve tAI ile birlikte değerlendirin.")
        elif hits:
            lines.append("Anlamlı kodon sayısı az olduğu için üçüncü konum bileşimine göre bir "
                         "yorum yapılmadı.")

        warning_lines = []
        small = [r for r in hits if abs(r["log2_or"]) < small_effect]
        if small:
            warning_lines.append(f"{len(small)} anlamlı kodonun etkisi küçük (|log2 OR| < "
                                 f"{small_effect}); q anlamlılığı söyler, büyüklüğü söylemez ve "
                                 f"büyük gruplarda önemsiz farklar bile anlamlı çıkabilir.")
        noisy = [r["codon"] for r in rows if abs(r["log2_or"]) >= _LARGE_EFFECT and not r["robust"]]
        if noisy:
            warning_lines.append(f"Büyük görünen ama anlamlı olmayan oranlar küçük örneklem "
                                 f"yanılsaması olabilir: {', '.join(noisy[:8])}.")
        rare = [r["codon"] for r in rows if r["a"] + r["c"] < min_count]
        if rare:
            warning_lines.append(f"İki grupta toplam {min_count} kezden az geçen kodonlarda güven "
                                 f"düşüktür: {', '.join(rare[:8])}.")
        for cub, name in ((self.a, name_a), (self.b, name_b)):
            gene, share = self._dominant_gene(cub)
            if len(cub) > 1 and share >= _DOMINANT_GENE:
                warning_lines.append(f"{name} grubunun kodonlarının %{100 * share:.0f}'i tek genden "
                                     f"({gene}) geliyor; sonuç bu gene bağlı olabilir.")
        lines += [f"Uyarı: {line}" for line in warning_lines]

        lines.append("Test bir ilişkiyi gösterir, nedeni değil: zenginleşme çeviri seçilimi, "
                     "mutasyon yanlılığı ya da GC bileşimi farkından gelebilir. Permütasyon genler "
                     "arası çeşitliliği hesaba katar, ama gruplar uzunluk ya da bölgesel GC "
                     "bakımından eşleştirilmediyse bu farklar sonuca girer.")
        return "\n".join(lines)

    # --- Grafik -----------------------------------------------------------------

    def plot_volcano(self, level: str = "family", q: float = 0.05, labels: bool = True,
                     ax: Optional[Axes] = None, file_name: Optional[str | Path] = None) -> Axes:
        """Yanardağ grafiği: x log2 OR, y −log10 Fisher q; anlamlılar grup rengindedir.

        Grafik Fisher sonucunu gösterir; gen düzeyinde doğrulama için ``robust``.

        Args:
            level: Karşılaştırma düzeyi.
            q: Anlamlılık eşiği (yatay kılavuz çizgisi).
            labels: True ise anlamlı kodonların adları yazılır.
            ax: Çizilecek eksen; verilmezse yeni şekil açılır.
            file_name: Verilirse şekil kaydedilir.
        """
        rows = self.table(level, q)
        ax = _plotting.axes(ax)
        colors = dict(zip((self.name_a, self.name_b), _plotting.PLOT_GROUP_COLORS))
        style = {k: v for k, v in _plotting.PLOT_POINT_STYLE.items() if k != "color"}
        handles = []
        for label, color, members in (
                (self.name_a, colors[self.name_a], [r for r in rows if r["favored"] == self.name_a]),
                (self.name_b, colors[self.name_b], [r for r in rows if r["favored"] == self.name_b]),
                ("anlamlı değil", "#b0b0b0", [r for r in rows if not r["significant"]])):
            if members:
                handles.append(ax.scatter([r["log2_or"] for r in members],
                                          [-r["log10_q"] for r in members],
                                          color=color, label=label, **style))
        ax.axhline(-math.log10(q), **_plotting.PLOT_GUIDE_STYLE)
        ax.axvline(0, **_plotting.PLOT_GUIDE_STYLE)
        if labels:
            for r in rows:
                if r["significant"]:
                    ax.annotate(r["codon"], (r["log2_or"], -r["log10_q"]), xytext=(2, 2),
                                textcoords="offset points", fontsize=_plotting.PLOT_LABEL_SIZE,
                                color="#333333")
        ax.set_xlabel(f"log2 OR ({self.name_b} ← → {self.name_a})")
        ax.set_ylabel("−log10 q")
        ax.set_title(f"Diferansiyel kodon kullanımı ({level})")
        _plotting.group_legend(ax, handles)
        return _plotting.save(ax, file_name)

    def __repr__(self) -> str:
        return f"DifferentialCodonUsage: {self.name_a} ({len(self.a)} gen) / {self.name_b} ({len(self.b)} gen)"
