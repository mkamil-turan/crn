"""Kodon kullanım yanlılığı (codon usage bias) — D. Referans setine dayalı indeksler.

Bir ``CUB`` nesnesinin kodonlarını bir referans gen setine (genellikle
yüksek ifade edilen genler) göre değerlendirir. Referans, karşılaştırma seti
ve optimal kodon seçimi C bölümündeki ``Synonymous`` ile aynıdır.

D. Referans setine dayalı indeksler:

====  ==========================================================  ===========================
 No   Ölçüm                                                       Metot
====  ==========================================================  ===========================
 18   w (göreli uyum): referansta N_kodon / ailenin en büyük N'i  ``relative_adaptiveness``
 19   CAI: w'lerin geometrik ortalaması, exp(Σ ln w / L)          ``cai``
 20   CBI: (N_opt − N_rastgele) / (N_toplam − N_rastgele)         ``cbi``
====  ==========================================================  ===========================

**w ve CAI (Sharp ve Li, 1987).** w yalnızca eş anlamlı kodonlar için
tanımlıdır; tek kodonlu amino asitler (standart kodda Met, Trp) her zaman
w = 1 alacakları için CAI'ye girmez, stoplar da girmez. Bu yüzden CAI'deki L
genin eş anlamlı kodon sayısıdır. Referansta hiç geçmeyen bir kodonun w'si 0
olur ve ln 0 tanımsızdır; Sharp ve Li'deki gibi böyle kodonların sayısı
``absent_count`` (varsayılan 0,5) kabul edilir. Referansta hiç kullanılmayan
bir amino asidin tüm kodonları böylece w = 1 alır.

**CBI (Bennetzen ve Hall, 1982).** Yalnızca en az bir optimal kodonu olan
amino asitlerin kodonları sayılır:

- N_toplam: bu amino asitlerin genin içindeki kodon sayısı,
- N_opt: bunlardan optimal olanların sayısı,
- N_rastgele: kodonlar rastgele seçilseydi beklenecek optimal kodon sayısı,
  Σ N_aa × (amino asidin optimal kodon sayısı / aile büyüklüğü).

CBI 1 ise yalnızca optimal kodonlar, 0 ise rastgele kullanım vardır; optimal
kodonlar rastgeleden az kullanılırsa negatif olur.

Kullanım::

    from crn.cub import CUB, Adaptation

    ada = Adaptation(cub, reference=CUB(genes=ribosomal_genes))
    ada.relative_adaptiveness()             # {kodon: w}
    ada.cai(gene="MEFV_NM_000243.3")
    ada.cbi(gene="MEFV_NM_000243.3", method="chi2")
"""
from __future__ import annotations

import math
from typing import (Dict, Optional, Sequence)

from crn.io.crn_io_const import CODONS
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_synonymous import (Synonymous, PREFERRED_RSCU, OPTIMAL_ALPHA)


# Referansta hiç geçmeyen kodonun w hesabında kabul edilen sayısı (Sharp ve Li, 1987).
ABSENT_COUNT = 0.5


class Adaptation:
    """Bir ``CUB`` nesnesinin referans gen setine uyumu.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.
        reference: w'lerin ve optimal kodonların belirleneceği referans
            (yüksek ifadeli) gen seti; verilmezse ``cub`` kullanılır (bu
            durumda indeksler genin kendisine göre ölçülür).
        background: Optimal kodonların ``method="chi2"`` ile seçiminde
            referansla karşılaştırılan gen seti; verilmezse ``cub``.

    Raises:
        TypeError, ValueError: ``Synonymous`` ile aynı koşullarda.
    """

    def __init__(self, cub: CUB, reference: Optional[CUB] = None,
                 background: Optional[CUB] = None):
        self.__synonymous = Synonymous(cub, reference, background)

    @property
    def synonymous(self) -> Synonymous:
        """Optimal kodon seçiminde kullanılan ``Synonymous`` nesnesi."""
        return self.__synonymous

    @property
    def cub(self) -> CUB:
        return self.__synonymous.cub

    @property
    def reference(self) -> CUB:
        return self.__synonymous.reference

    @property
    def background(self) -> CUB:
        return self.__synonymous.background

    # --- D. Referans setine dayalı indeksler ------------------------------------

    def relative_adaptiveness(self, absent_count: float = ABSENT_COUNT,
                              codons: Sequence[str] = CODONS) -> Dict[str, float]:
        """18. w: referans setteki (tüm genleri birlikte) kodon sayısı / ailesinin en büyük sayısı.

        Yalnızca eş anlamlı kodonlar için döner. Kodon sayısı ile kodon oranı
        ya da RSCU aynı aile içinde orantılı olduğundan w üçünden de aynı
        çıkar.

        Args:
            absent_count: Referansta hiç geçmeyen kodonun sayısı yerine
                kullanılacak değer (0'dan büyük olmalı).

        Raises:
            ValueError: absent_count 0 ya da negatifse.
        """
        if absent_count <= 0:
            raise ValueError(f"absent_count 0'dan büyük olmalı, {absent_count!r} verildi")
        reference = self.reference
        counts = reference.codon_counts(None, codons)
        synonymous = set(reference.synonymous_codons(codons))

        w = {}
        for family in reference.codon_families().values():
            family = [codon for codon in family if codon in synonymous]
            if not family:
                continue
            values = {codon: (counts[codon] or absent_count) for codon in family}
            largest = max(values.values())
            w.update({codon: n / largest for codon, n in values.items()})
        return {codon: w[codon] for codon in codons if codon in w}

    def cai(self, gene: Optional[str] = None,
            w: Optional[Dict[str, float]] = None,
            absent_count: float = ABSENT_COUNT,
            codons: Sequence[str] = CODONS) -> float:
        """19. CAI: genin eş anlamlı kodonlarının w değerlerinin geometrik ortalaması.

        Args:
            w: Kodon -> w; verilmezse ``relative_adaptiveness(absent_count)``.

        Returns:
            0-1 arası değer; genin hiç eş anlamlı kodonu yoksa 0.

        Raises:
            ValueError: w içinde 0 ya da negatif bir değer varsa.
        """
        if w is None:
            w = self.relative_adaptiveness(absent_count, codons)
        counts = self.cub.codon_counts(gene, [codon for codon in codons if codon in w])
        total = sum(counts.values())
        if not total:
            return 0.0
        log_sum = 0.0
        for codon, n in counts.items():
            if n:
                if w[codon] <= 0:
                    raise ValueError(f"{codon} kodonunun w'si {w[codon]!r}; CAI için 0'dan büyük olmalı")
                log_sum += n * math.log(w[codon])
        return math.exp(log_sum / total)

    def cbi(self, gene: Optional[str] = None,
            optimal_codons: Optional[Sequence[str]] = None,
            method: str = "rscu",
            threshold: float = PREFERRED_RSCU,
            alpha: float = OPTIMAL_ALPHA,
            codons: Sequence[str] = CODONS) -> float:
        """20. CBI: (N_opt − N_rastgele) / (N_toplam − N_rastgele).

        Args:
            optimal_codons: Optimal kodonlar; verilmezse
                ``Synonymous.optimal_codons(method, threshold, alpha)``.
            method, threshold, alpha: Optimal kodon seçimi (bkz.
                ``Synonymous.optimal_codons``).

        Returns:
            1 yalnızca optimal kodon, 0 rastgele kullanım; payda 0 ise 0.

        Raises:
            ValueError: method geçersizse.
        """
        if optimal_codons is None:
            optimal_codons = self.__synonymous.optimal_codons(method, threshold, alpha, codons)
        optimal_codons = set(optimal_codons)

        cub = self.cub
        table = cub.genetic_code.table
        counts = cub.codon_counts(gene, codons)
        aa_counts = cub.amino_acid_counts(gene, codons=codons)
        synonymous = set(cub.synonymous_codons(codons))

        n_total = n_optimal = n_random = 0.0
        for aa, family in cub.codon_families().items():
            optimal = [codon for codon in family if codon in optimal_codons and codon in synonymous]
            if not optimal:
                continue
            n_total += aa_counts[aa]
            n_optimal += sum(counts[codon] for codon in optimal)
            n_random += aa_counts[aa] * len(optimal) / len(family)
        denominator = n_total - n_random
        return (n_optimal - n_random) / denominator if denominator else 0.0

    def __repr__(self) -> str:
        return f"Adaptation({self.__synonymous!r})"
