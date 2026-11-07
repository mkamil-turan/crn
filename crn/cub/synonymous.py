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
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple

from crn.cub.core import CUB

if TYPE_CHECKING:
    from matplotlib.axes import Axes


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

    # --- Çizimler --------------------------------------------------------------------

    def _synonymous_layout(self) -> Tuple[List[str], List[str]]:
        """Eş anlamlı kodonlar (amino asit sırasıyla gruplu) ve her kodonun amino asidi.

        Tek kodonlu amino asitler (Met, Trp, stop) RSCU'su hep 1 olduğundan dışta.
        """
        table = self.cub.genetic_code.table
        families = self.cub.genetic_code.families
        codons: List[str] = []
        amino_acids: List[str] = []
        for aa in sorted(families):
            family = families[aa]
            if len(family) < 2:
                continue
            for codon in family:
                codons.append(codon)
                amino_acids.append(table[codon])
        return codons, amino_acids

    def plot_rscu(self, gene: Optional[str] = None, ax: Optional[Axes] = None,
                  file_name: Optional[str | Path] = None) -> Axes:
        """15. RSCU çubuk grafiği: eş anlamlı kodonlar amino aside göre gruplanır.

        Çubuklar ardışık amino asit grupları için iki renkte boyanır; 1
        çizgisi tarafsız kullanımı (yanlılık yok), 1 üstü tercihi gösterir.
        ``gene`` verilmezse tüm genler birlikte (havuz) çizilir.
        """
        from crn.cub._plotting import PLOT_GUIDE_STYLE, PLOT_LABEL_SIZE, axes, save
        ax = axes(ax, figsize=(10, 4))
        codons, amino_acids = self._synonymous_layout()
        rscu = self.rscu(gene)
        values = [rscu[codon] for codon in codons]

        group_colors = ("#2a6fb0", "#8cb8dc")
        group_index: Dict[str, int] = {}
        colors = []
        for aa in amino_acids:
            group_index.setdefault(aa, len(group_index))
            colors.append(group_colors[group_index[aa] % 2])

        positions = range(len(codons))
        ax.bar(positions, values, color=colors, width=0.8, zorder=3)
        ax.axhline(1.0, **PLOT_GUIDE_STYLE)
        ax.set_xticks(list(positions))
        ax.set_xticklabels(codons, rotation=90, fontsize=PLOT_LABEL_SIZE,
                           family="monospace")
        # Her amino asit grubunun ortasına amino asit harfini yaz.
        seen = 0
        for aa in dict.fromkeys(amino_acids):
            idx = [i for i, a in enumerate(amino_acids) if a == aa]
            ax.text(sum(idx) / len(idx), ax.get_ylim()[1], aa, ha="center",
                    va="bottom", fontsize=8, fontweight="bold")
            seen += 1
        ax.set(xlim=(-0.6, len(codons) - 0.4), ylabel="RSCU",
               title=f"RSCU — {gene or 'tüm genler'}")
        ax.margins(y=0.08)
        return save(ax, file_name)

    def plot_rscu_heatmap(self, ax: Optional[Axes] = None, genes: Optional[Sequence[str]] = None,
                          cmap: str = "YlGnBu", groups=None,
                          file_name: Optional[str | Path] = None) -> Axes:
        """15. RSCU ısı haritası: satırlar genler, sütunlar eş anlamlı kodonlar.

        Türler/genler arası tercih örüntüsünü ve kümelenmeyi gösterir
        (literatürde sık sık kümelenmiş ısı haritası olarak). ``genes``
        verilmezse ``cub``'daki tüm genler. Satırlar gen sembolüyle yazılır.
        ``groups`` verilirse satırlar gruplara göre sıralanır, sembol
        renkleri grubu gösterir ve gruplar arasına ayraç çizilir.
        """
        import matplotlib.patches as mpatches
        from crn.cub._plotting import (PLOT_GENE_LABEL_SIZE, PLOT_LABEL_SIZE, axes, gene_groups,
                                       gene_symbols, group_colors, group_legend, save)
        names = list(genes) if genes is not None else list(self.cub.names)
        assigned = gene_groups(self.cub, names, groups)
        if assigned is not None:
            colors = group_colors(assigned)
            order = sorted(range(len(names)), key=lambda i: list(colors).index(assigned[i]))
            names, assigned = [names[i] for i in order], [assigned[i] for i in order]
        codons, amino_acids = self._synonymous_layout()
        matrix = [[self.rscu(name)[codon] for codon in codons] for name in names]

        many = len(names) > 40
        ax = axes(ax, figsize=(max(8, len(codons) * 0.18),
                               max(3, len(names) * (0.09 if many else 0.22))))
        image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0.0,
                          interpolation="nearest")
        ax.set_xticks(range(len(codons)))
        ax.set_xticklabels(codons, rotation=90, fontsize=PLOT_LABEL_SIZE, family="monospace")
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(gene_symbols(self.cub, names),
                           fontsize=PLOT_GENE_LABEL_SIZE if many else PLOT_LABEL_SIZE)
        if assigned is not None:
            for tick, group in zip(ax.get_yticklabels(), assigned):
                tick.set_color(colors[group])
            for i in range(1, len(names)):
                if assigned[i] != assigned[i - 1]:
                    ax.axhline(i - 0.5, color="white", linewidth=1.2)
        # Amino asit grupları arasına ince ayraç çizgileri.
        for i in range(1, len(codons)):
            if amino_acids[i] != amino_acids[i - 1]:
                ax.axvline(i - 0.5, color="white", linewidth=0.6)
        ax.figure.colorbar(image, ax=ax, label="RSCU", fraction=0.025, pad=0.02)
        if assigned is not None:
            group_legend(ax, [mpatches.Patch(color=c, label=g) for g, c in colors.items()])
            ax.get_legend().set_bbox_to_anchor((1.08, 1.0))   # renk çubuğunun sağına
        ax.set_title("RSCU ısı haritası")
        return save(ax, file_name)

    def __repr__(self) -> str:
        def label(other: CUB) -> str:
            return "kendisi" if other is self.cub else repr(other)
        return (f"Synonymous({self.cub!r}, referans: {label(self.reference)}, "
                f"karşılaştırma: {label(self.background)})")
