"""Kodon kullanım yanlılığı (codon usage bias) — C. Sinonim kodon tercihi.

Bir ``CUB`` nesnesinin kodon ve amino asit sayımlarından
(``CUB.codon_counts``, ``CUB.amino_acid_counts``) ve eş anlamlı kodon
ailelerinden (``CUB.codon_families``) eş anlamlı kodon tercihini hesaplar.

C. Sinonim kodon tercihi:

====  ================================================  ======================
 No   Ölçüm                                             Metot
====  ================================================  ======================
 14   Kodon oranı: N_kodon / N_amino asit               ``codon_ratio``
 15   RSCU: N_kodon / (N_amino asit / aile büyüklüğü)   ``rscu``
 16   Tercih edilen kodon: RSCU > 1,6 (< 0,6 az)        ``preferred_codons``,
                                                        ``underused_codons``
 17   Fop: optimal kodon / eş anlamlı kodon             ``fop``
====  ================================================  ======================

RSCU'da 1 tarafsız kullanımı, 1'in üstü tercihi gösterir; bir ailenin RSCU
değerlerinin toplamı aile büyüklüğüne eşittir. Hiç kullanılmayan amino
asidin kodonlarında kodon oranı ve RSCU 0'dır. Tek kodonlu amino asitlerin
(standart kodda Met, Trp) RSCU'su kullanıldıklarında her zaman 1'dir; bu
yüzden tercih edilen/az kullanılan kodon listelerine girmezler.

Fop için **optimal kodonlar** bir referans gen setinden gelir. Referans
olarak genellikle yüksek ifade edilen genlerden (ör. ribozomal proteinler)
kurulan ikinci bir ``CUB`` verilir. Referans verilmezse ``cub``'ın kendisi
kullanılır; bu durumda Fop genin kendi tercihini kendisiyle ölçer ve tek
gende anlamı sınırlıdır. Optimal kodonlar iki yöntemle seçilir (``method``):

- ``"rscu"`` (varsayılan): referansın (tüm genleri birlikte) RSCU'su
  eşiği (1,6) aşan eş anlamlı kodonlar.
- ``"chi2"`` (Ikemura; CodonW): referans ile karşılaştırma seti
  (``background``, ör. düşük ifadeli genler) arasında, kendi amino asidi
  içindeki payı referansta anlamlı biçimde (2×2 χ², ``alpha``) yüksek olan
  eş anlamlı kodonlar. ``background`` verilmezse ``cub`` kullanılır.

Optimal kodon listesi doğrudan da verilebilir. Fop'un paydası da seçilebilir
(``denominator``):

- ``"synonymous"`` (varsayılan): tüm eş anlamlı kodonlar.
- ``"optimal"`` (CodonW): yalnızca en az bir optimal kodonu olan amino
  asitlerin kodonları.

Kullanım::

    from crn.cub import CUB, Synonymous

    cub = CUB(isoforms)
    syn = Synonymous(cub, reference=CUB(ribosomal_genes),
                     background=CUB(low_expression_genes))
    syn.rscu(gene="MEFV_NM_000243.3")
    syn.preferred_codons(gene="MEFV_NM_000243.3")
    syn.fop(gene="MEFV_NM_000243.3")                       # rscu, tüm eş anlamlılar
    syn.fop(gene="MEFV_NM_000243.3", method="chi2", denominator="optimal")
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from crn.cub.core import CUB


# Tercih edilen ve az kullanılan kodon için yaygın RSCU eşikleri.
PREFERRED_RSCU = 1.6
UNDERUSED_RSCU = 0.6
# χ² yöntemiyle optimal kodon seçiminde anlamlılık düzeyi.
OPTIMAL_ALPHA = 0.01

# Optimal kodon seçim yöntemleri ve Fop paydaları.
OPTIMAL_METHODS = ("rscu", "chi2")
FOP_DENOMINATORS = ("synonymous", "optimal")


def check_choice(name: str, value: str, choices: Sequence[str]) -> None:
    if value not in choices:
        raise ValueError(f"{name} {', '.join(map(repr, choices))} olmalı, {value!r} verildi")


def check_cub(name: str, value, cub: Optional[CUB] = None) -> None:
    """value bir CUB değilse TypeError, genetik kodu cub'ınkinden farklıysa ValueError."""
    if not isinstance(value, CUB):
        raise TypeError(f"{name} bir CUB olmalı, {type(value).__name__} verildi")
    if cub is not None and value.genetic_code is not cub.genetic_code:
        raise ValueError(f"{name}'ın genetik kodu ({value.genetic_code.name}) "
                         f"cub'ınkinden ({cub.genetic_code.name}) farklı")


def _chi2_p(a: int, b: int, c: int, d: int) -> float:
    """2×2 tablonun [[a, b], [c, d]] Pearson χ² testinin p değeri (1 sd); boş satır/sütunda 1."""
    n = a + b + c + d
    margins = (a + b) * (c + d) * (a + c) * (b + d)
    if not margins:
        return 1.0
    chi2 = n * (a * d - b * c) ** 2 / margins
    return math.erfc(math.sqrt(chi2 / 2))


class Synonymous:
    """Bir ``CUB`` nesnesinin eş anlamlı kodon tercihi.

    Metotlar ``gene`` verilmezse tüm genleri birlikte, verilirse yalnızca o
    geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.
        reference: Optimal kodonların belirleneceği (yüksek ifadeli) gen seti;
            verilmezse ``cub``.
        background: ``method="chi2"`` iken referansla karşılaştırılan gen seti
            (ör. düşük ifadeli genler); verilmezse ``cub``.

    Raises:
        TypeError: cub, reference ya da background bir CUB değilse.
        ValueError: reference ya da background'un genetik kodu farklıysa.
    """

    def __init__(self, cub: CUB, reference: Optional[CUB] = None,
                 background: Optional[CUB] = None):
        check_cub("cub", cub)
        self.cub = cub
        self.reference = cub if reference is None else reference
        self.background = cub if background is None else background
        check_cub("reference", self.reference, cub)
        check_cub("background", self.background, cub)
        self._optimal: Dict[Tuple, List[str]] = {}

    # --- C. Sinonim kodon tercihi -----------------------------------------------

    def codon_ratio(self, gene: Optional[str] = None) -> Dict[str, float]:
        """14. Kodon oranı: kodonun kendi amino asidi içindeki payı (N_kodon / N_amino asit).

        Bir ailenin oranlarının toplamı 1'dir (amino asit hiç yoksa hepsi 0).
        """
        table = self.cub.genetic_code.table
        counts = self.cub.codon_counts(gene)
        aa_counts = self.cub.amino_acid_counts(gene)
        return {codon: (counts[codon] / aa_counts[table[codon]] if aa_counts[table[codon]] else 0.0)
                for codon in self.cub.sense_codons()}

    def rscu(self, gene: Optional[str] = None) -> Dict[str, float]:
        """15. RSCU: N_kodon / (N_amino asit / aile büyüklüğü) = kodon oranı × aile büyüklüğü."""
        table = self.cub.genetic_code.table
        families = self.cub.genetic_code.families
        return {codon: ratio * len(families[table[codon]])
                for codon, ratio in self.codon_ratio(gene).items()}

    def preferred_codons(self, gene: Optional[str] = None,
                         threshold: float = PREFERRED_RSCU) -> List[str]:
        """16. Tercih edilen kodonlar: RSCU'su ``threshold``'u aşan eş anlamlı kodonlar."""
        synonymous = self.cub.genetic_code.synonymous_codons
        return [codon for codon, value in self.rscu(gene).items()
                if codon in synonymous and value > threshold]

    def underused_codons(self, gene: Optional[str] = None,
                         threshold: float = UNDERUSED_RSCU) -> List[str]:
        """16. Az kullanılan kodonlar: RSCU'su ``threshold``'un altındaki eş anlamlı kodonlar.

        Hiç kullanılmayan amino asitlerin kodonları (RSCU 0) dahil edilmez.
        """
        table = self.cub.genetic_code.table
        synonymous = self.cub.genetic_code.synonymous_codons
        aa_counts = self.cub.amino_acid_counts(gene)
        return [codon for codon, value in self.rscu(gene).items()
                if codon in synonymous and aa_counts[table[codon]] and value < threshold]

    def optimal_codons(self, method: str = "rscu",
                       threshold: float = PREFERRED_RSCU,
                       alpha: float = OPTIMAL_ALPHA) -> List[str]:
        """Fop/CBI için optimal kodonlar, ``CODONS`` sırasıyla (sonuç önbelleğe alınır).

        Args:
            method: ``"rscu"``: referansın (tüm genleri birlikte) RSCU'su
                ``threshold``'u aşan eş anlamlı kodonlar. ``"chi2"`` (Ikemura;
                CodonW): her eş anlamlı kodon için referans ve background'da
                [kodon, ailenin diğer kodonları] sayılarından 2×2 χ²; payı
                referansta yüksek ve p < ``alpha`` olanlar.

        Raises:
            ValueError: method geçersizse ya da ``"chi2"`` için reference ile
                background aynı setse.
        """
        check_choice("method", method, OPTIMAL_METHODS)
        key = (method, threshold) if method == "rscu" else (method, alpha)
        if key not in self._optimal:
            self._optimal[key] = self._find_optimal(method, threshold, alpha)
        return list(self._optimal[key])

    def _find_optimal(self, method: str, threshold: float, alpha: float) -> List[str]:
        if method == "rscu":
            return Synonymous(self.reference).preferred_codons(None, threshold)
        if self.reference is self.background:
            raise ValueError("method='chi2' için reference ve background farklı gen setleri olmalı")
        table = self.cub.genetic_code.table
        high, low = self.reference.codon_counts(), self.background.codon_counts()
        high_aa, low_aa = self.reference.amino_acid_counts(), self.background.amino_acid_counts()
        optimal = []
        for codon in self.cub.synonymous_codons():
            aa = table[codon]
            a, b = high[codon], high_aa[aa] - high[codon]
            c, d = low[codon], low_aa[aa] - low[codon]
            if not (a + b and c + d) or a / (a + b) <= c / (c + d):
                continue
            if _chi2_p(a, b, c, d) < alpha:
                optimal.append(codon)
        return optimal

    def fop(self, gene: Optional[str] = None,
            optimal_codons: Optional[Sequence[str]] = None,
            method: str = "rscu",
            denominator: str = "synonymous",
            threshold: float = PREFERRED_RSCU,
            alpha: float = OPTIMAL_ALPHA) -> float:
        """17. Fop: optimal kodon sayısı / payda (0-1; payda 0 ise 0).

        Args:
            optimal_codons: Optimal kodonlar; verilmezse
                ``optimal_codons(method, threshold, alpha)``.
            denominator: ``"synonymous"``: tüm eş anlamlı kodonlar;
                ``"optimal"`` (CodonW): en az bir optimal kodonu olan amino
                asitlerin kodonları.

        Raises:
            ValueError: method ya da denominator geçersizse.
        """
        check_choice("denominator", denominator, FOP_DENOMINATORS)
        if optimal_codons is None:
            optimal_codons = self.optimal_codons(method, threshold, alpha)
        optimal_codons = set(optimal_codons)

        pool = self.cub.synonymous_codons()
        if denominator == "optimal":
            table = self.cub.genetic_code.table
            amino_acids = {table[codon] for codon in optimal_codons if codon in table}
            pool = [codon for codon in pool if table[codon] in amino_acids]
        counts = self.cub.codon_counts(gene, pool)
        total = sum(counts.values())
        optimal = sum(n for codon, n in counts.items() if codon in optimal_codons)
        return optimal / total if total else 0.0

    def __repr__(self) -> str:
        def label(other: CUB) -> str:
            return "kendisi" if other is self.cub else repr(other)
        return (f"Synonymous({self.cub!r}, referans: {label(self.reference)}, "
                f"karşılaştırma: {label(self.background)})")
