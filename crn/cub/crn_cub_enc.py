"""Kodon kullanım yanlılığı (codon usage bias) — F. Etkin kodon sayısı ve grafikler.

Bir ``CUB`` nesnesinin kodon sayımlarından (``CUB.codon_counts``) ve eş
anlamlı kodon ailelerinden (``CUB.codon_families``) etkin kodon sayısını;
B'deki ``Composition``'dan alınan GC3s, GC12, GC3 ve 3. konum baz
sayılarından ENC–GC3s, nötrallik ve PR2 grafiklerinin verilerini hesaplar.

F. Etkin kodon sayısı ve grafikler:

====  ==============================================================  ==========================
 No   Ölçüm                                                           Metot
====  ==============================================================  ==========================
 23   ENC (Nc, Wright 1990): F = (nΣp² − 1) / (n − 1);                ``enc``, ``homozygosity``
      Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
 24   Beklenen ENC: 2 + s + 29 / (s² + (1 − s)²), s = GC3s            ``expected_enc``
 25   ENC oranı: (ENC_beklenen − ENC_gözlenen) / ENC_beklenen         ``enc_ratio``
 26   Nötrallik grafiği: GC12 (y) ~ GC3 (x) regresyonu                ``neutrality``
 27   PR2 grafiği: x = G3 / (G3 + C3), y = A3 / (A3 + U3)             ``pr2``, ``pr2_points``
====  ==============================================================  ==========================

**ENC.** Her amino asit için F (homozigotluk) hesaplanır; aynı dejenerelik
sınıfındaki (aile büyüklüğü 2, 3, 4, 6) amino asitlerin F'leri ortalanır ve
Nc = (tek kodonlu amino asit sayısı) + Σ_k (sınıftaki amino asit sayısı) / F̄_k
bulunur; standart kodda bu Wright'ın formülüdür. En az iki kez kullanılmamış
(n < 2) ya da F'si 0 çıkan amino asitler hesaba girmez. Bir sınıfta hiç F
kalmazsa: 3'lü sınıf (standart kodda yalnızca Ile) için F3 = (F2 + F4) / 2
alınır (Wright), diğer sınıflarda ENC hesaplanamaz ve ``nan`` döner. Sonuç
amino asit kodlayan kodon sayısını (standart kodda 61) aşarsa o değere
indirilir.

**Grafikler.** ``*_points`` ve ``neutrality`` çizilecek verileri döner;
``plot_*`` metotları bunları literatürdeki biçimleriyle matplotlib ile
çizer. Nötrallik ve ENC–GC3s grafikleri genler arası karşılaştırmadır: her
gen bir noktadır. Nötrallik grafiğinde eğim 1'e yakınsa mutasyon baskısı,
0'a yakınsa seçilim baskın kabul edilir. PR2 varsayılan olarak dört katlı
dejenere kodonlarda (``CUB.fourfold_codons``; Sueoka) hesaplanır; (0,5; 0,5)
merkezi A=U ve G=C dengesidir.

==========================  ===================================================================
Metot                       Grafik
==========================  ===================================================================
``plot_enc``                ENC–GC3s (Wright 1990): noktalar ve beklenen ENC eğrisi; x 0-1,
                            y 20-61
``plot_enc_ratio``          ENC oranı frekans dağılımı (histogram); 0 çizgisi beklenen
``plot_neutrality``         Nötrallik (Sueoka 1988): GC12 ~ GC3, regresyon doğrusu,
                            denklem ve R², y = x köşegeni
``plot_pr2``                PR2 (Sueoka 1995): x 0-1, y 0-1, merkezde 0,5 çizgileri
==========================  ===================================================================

Her ``plot_*`` metodu ``ax`` verilirse onun üzerine, verilmezse yeni bir
şekle çizer ve ``Axes`` döner; ``file_name`` verilirse şekli kaydeder.

Kullanım::

    from crn.cub import CUB, EffectiveCodons

    eff = EffectiveCodons(cub)
    eff.enc(gene="MEFV_NM_000243.3")
    eff.enc_ratio(gene="MEFV_NM_000243.3")
    eff.neutrality()                         # {"points": [...], "slope": ..., ...}
    eff.pr2_points()
    eff.plot_enc(file_name="enc.png")
    eff.plot_neutrality(labels=True)
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import (Dict, List, Optional, Sequence, Tuple)

import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from crn.io.crn_io_const import CODONS
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_composition import Composition


# Grafik ayarları: eksen sınırları, eğri çözünürlüğü, nokta ve çizgi biçimleri.
ENC_AXIS_LIMITS = (20.0, 61.0)
ENC_CURVE_STEPS = 200
ENC_RATIO_BIN_WIDTH = 0.05
PLOT_POINT_STYLE = {"s": 18, "color": "#2a6fb0", "edgecolor": "none", "alpha": 0.8, "zorder": 3}
PLOT_CURVE_STYLE = {"color": "#c0392b", "linewidth": 1.4}
PLOT_GUIDE_STYLE = {"color": "#7f7f7f", "linestyle": "--", "linewidth": 0.8}
PLOT_LABEL_SIZE = 7


def _axes(ax: Optional[Axes]) -> Axes:
    """Verilen ``Axes``'i, yoksa yeni bir şeklin ``Axes``'ini döner."""
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4.5))
    return ax


def _finish(ax: Axes, points: Sequence[Tuple[float, float, str]], labels: bool,
            file_name: Optional[str | Path]) -> Axes:
    """Noktaları çizer, istenirse gen adlarını yazar ve şekli kaydeder."""
    ax.scatter([x for x, _, _ in points], [y for _, y, _ in points], **PLOT_POINT_STYLE)
    if labels:
        for x, y, name in points:
            ax.annotate(name, (x, y), xytext=(3, 3), textcoords="offset points",
                        fontsize=PLOT_LABEL_SIZE)
    if file_name is not None:
        ax.figure.savefig(file_name, dpi=300, bbox_inches="tight")
    return ax


def _linear_regression(points: Sequence[Tuple[float, float]]) -> Dict[str, float]:
    """En küçük kareler doğrusu y = eğim·x + kesişim ve r²; iki noktadan azsa ``nan``."""
    n = len(points)
    if n < 2:
        return {"slope": math.nan, "intercept": math.nan, "r2": math.nan}
    mean_x = sum(x for x, _ in points) / n
    mean_y = sum(y for _, y in points) / n
    sxx = sum((x - mean_x) ** 2 for x, _ in points)
    syy = sum((y - mean_y) ** 2 for _, y in points)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in points)
    if not sxx:
        return {"slope": math.nan, "intercept": math.nan, "r2": math.nan}
    slope = sxy / sxx
    return {"slope": slope, "intercept": mean_y - slope * mean_x,
            "r2": sxy * sxy / (sxx * syy) if syy else math.nan}


class EffectiveCodons:
    """Bir ``CUB`` nesnesinin etkin kodon sayısı ve CUB grafiklerinin verileri.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar; ``*_points`` ve
    ``neutrality`` her gen için ayrı bir nokta verir.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        self.__composition = Composition(cub)

    @property
    def cub(self) -> CUB:
        return self.__composition.cub

    def _names(self) -> List[str]:
        return [gene.name for gene in self.cub.genes]

    # --- 23-25. ENC ------------------------------------------------------------------

    def homozygosity(self, gene: Optional[str] = None,
                     codons: Sequence[str] = CODONS) -> Dict[str, float]:
        """23. Her eş anlamlı amino asidin F değeri: (nΣp² − 1) / (n − 1).

        n < 2 olan amino asitler dönmez.
        """
        counts = self.cub.codon_counts(gene, codons)
        result = {}
        for aa, family in self.cub.codon_families().items():
            if len(family) < 2:
                continue
            n = sum(counts.get(codon, 0) for codon in family)
            if n < 2:
                continue
            sum_p2 = sum((counts.get(codon, 0) / n) ** 2 for codon in family)
            result[aa] = (n * sum_p2 - 1) / (n - 1)
        return result

    def enc(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """23. Etkin kodon sayısı Nc (Wright 1990); 20 aşırı yanlılık, 61 yanlılık yok.

        Returns:
            Nc; bir dejenerelik sınıfı için F bulunamazsa ``nan``.
        """
        families = self.cub.codon_families()
        f_values = {aa: f for aa, f in self.homozygosity(gene, codons).items() if f > 0}

        classes: Dict[int, List[str]] = {}
        for aa, family in families.items():
            if family:
                classes.setdefault(len(family), []).append(aa)

        means: Dict[int, float] = {}
        for size, amino_acids in classes.items():
            if size == 1:
                continue
            values = [f_values[aa] for aa in amino_acids if aa in f_values]
            if values:
                means[size] = sum(values) / len(values)
        if 3 in classes and 3 not in means and 2 in means and 4 in means:
            means[3] = (means[2] + means[4]) / 2

        nc = float(len(classes.get(1, [])))
        for size, amino_acids in classes.items():
            if size == 1:
                continue
            if size not in means:
                return math.nan
            nc += len(amino_acids) / means[size]
        return min(nc, float(len(self.cub.sense_codons())))

    @staticmethod
    def expected_enc(gc3s: float) -> float:
        """24. Yanlılık yalnızca GC3s'ten gelseydi beklenecek Nc: 2 + s + 29 / (s² + (1 − s)²)."""
        return 2 + gc3s + 29 / (gc3s ** 2 + (1 - gc3s) ** 2)

    def enc_ratio(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """25. ENC oranı: (ENC_beklenen − ENC_gözlenen) / ENC_beklenen."""
        expected = self.expected_enc(self.__composition.gc3s(gene, codons))
        return (expected - self.enc(gene, codons)) / expected

    def enc_points(self, codons: Sequence[str] = CODONS) -> List[Dict[str, float | str]]:
        """24-25. ENC–GC3s grafiği: her gen için GC3s, gözlenen ve beklenen ENC, ENC oranı."""
        points = []
        for name in self._names():
            gc3s = self.__composition.gc3s(name, codons)
            observed = self.enc(name, codons)
            expected = self.expected_enc(gc3s)
            points.append({"gene": name, "gc3s": gc3s, "enc": observed,
                           "expected_enc": expected,
                           "enc_ratio": (expected - observed) / expected})
        return points

    # --- 26-27. Grafikler ------------------------------------------------------------

    def neutrality(self, codons: Sequence[str] = CODONS) -> Dict[str, object]:
        """26. Nötrallik grafiği: her gen için (GC3, GC12) ve GC12 = eğim·GC3 + kesişim.

        Returns:
            ``{"points": [{"gene", "gc3", "gc12"}, ...], "slope", "intercept",
            "r2"}``; iki genden azsa doğru ``nan``.
        """
        points = [{"gene": name, "gc3": self.__composition.gc3(name, codons),
                   "gc12": self.__composition.gc12(name, codons)} for name in self._names()]
        fit = _linear_regression([(p["gc3"], p["gc12"]) for p in points])
        return {"points": points, **fit}

    def pr2(self, gene: Optional[str] = None,
            codons: Optional[Sequence[str]] = None) -> Dict[str, float]:
        """27. PR2 konumu: x = G3 / (G3 + C3), y = A3 / (A3 + U3).

        Args:
            codons: 3. konumu sayılacak kodonlar; verilmezse dört katlı
                dejenere kodonlar (``CUB.fourfold_codons``).

        Returns:
            ``{"x": G3/(G3+C3), "y": A3/(A3+U3)}``; payda 0 ise 0,5.
        """
        if codons is None:
            codons = self.cub.fourfold_codons()
        third = self.__composition.nucleotide_counts(gene, 3, codons)
        gc, au = third["G"] + third["C"], third["A"] + third["U"]
        return {"x": third["G"] / gc if gc else 0.5, "y": third["A"] / au if au else 0.5}

    def pr2_points(self, codons: Optional[Sequence[str]] = None) -> List[Dict[str, float | str]]:
        """27. PR2 grafiği: her gen için (x, y)."""
        return [{"gene": name, **self.pr2(name, codons)} for name in self._names()]

    # --- Çizimler --------------------------------------------------------------------

    def plot_enc(self, codons: Sequence[str] = CODONS, ax: Optional[Axes] = None,
                 labels: bool = False, file_name: Optional[str | Path] = None) -> Axes:
        """24. ENC–GC3s grafiği: genler ve Wright'ın beklenen ENC eğrisi.

        Eğrinin altında kalan genlerde yanlılık yalnızca GC3s ile açıklanamaz
        (seçilim gibi başka etkenler). ENC'si ``nan`` olan genler çizilmez.
        """
        ax = _axes(ax)
        s = [i / ENC_CURVE_STEPS for i in range(ENC_CURVE_STEPS + 1)]
        ax.plot(s, [self.expected_enc(v) for v in s], label="Beklenen ENC", **PLOT_CURVE_STYLE)
        points = [(p["gc3s"], p["enc"], p["gene"]) for p in self.enc_points(codons)
                  if not math.isnan(p["enc"])]
        ax.set(xlim=(0, 1), ylim=ENC_AXIS_LIMITS, xlabel="GC3s", ylabel="ENC",
               title="ENC–GC3s")
        ax.legend(loc="lower center", frameon=False)
        return _finish(ax, points, labels, file_name)

    def plot_enc_ratio(self, codons: Sequence[str] = CODONS, ax: Optional[Axes] = None,
                       bin_width: float = ENC_RATIO_BIN_WIDTH,
                       file_name: Optional[str | Path] = None) -> Axes:
        """25. ENC oranının frekans dağılımı; 0'dan uzaklaşan genler GC3s dışı etki taşır."""
        ax = _axes(ax)
        ratios = [p["enc_ratio"] for p in self.enc_points(codons) if not math.isnan(p["enc_ratio"])]
        if ratios:
            low = math.floor(min(ratios) / bin_width) * bin_width
            high = math.ceil(max(ratios) / bin_width) * bin_width
            n_bins = max(1, round((high - low) / bin_width))
            ax.hist(ratios, bins=[low + i * bin_width for i in range(n_bins + 1)],
                    color=PLOT_POINT_STYLE["color"], edgecolor="white")
        ax.axvline(0, **PLOT_GUIDE_STYLE)
        ax.set(xlabel="(ENC beklenen − ENC gözlenen) / ENC beklenen", ylabel="Gen sayısı",
               title="ENC oranı dağılımı")
        if file_name is not None:
            ax.figure.savefig(file_name, dpi=300, bbox_inches="tight")
        return ax

    def plot_neutrality(self, codons: Sequence[str] = CODONS, ax: Optional[Axes] = None,
                        labels: bool = False, file_name: Optional[str | Path] = None) -> Axes:
        """26. Nötrallik grafiği: GC12 ~ GC3, regresyon doğrusu, denklem, R² ve y = x.

        Tek gende regresyon ``nan`` olduğundan yalnızca nokta çizilir.
        """
        ax = _axes(ax)
        data = self.neutrality(codons)
        points = [(p["gc3"], p["gc12"], p["gene"]) for p in data["points"]]
        values = [v for x, y, _ in points for v in (x, y)]
        low, high = (min(values), max(values)) if values else (0.0, 1.0)
        pad = (high - low) * 0.1 or 0.05
        low, high = max(0.0, low - pad), min(1.0, high + pad)
        ax.plot([low, high], [low, high], label="y = x", **PLOT_GUIDE_STYLE)
        if not math.isnan(data["slope"]):
            ax.plot([low, high], [data["slope"] * v + data["intercept"] for v in (low, high)],
                    label=f"y = {data['slope']:.4f}x + {data['intercept']:.4f}, "
                          f"R² = {data['r2']:.4f}", **PLOT_CURVE_STYLE)
        ax.set(xlim=(low, high), ylim=(low, high), xlabel="GC3", ylabel="GC12",
               title="Nötrallik grafiği")
        ax.legend(loc="upper left", frameon=False, fontsize=8)
        return _finish(ax, points, labels, file_name)

    def plot_pr2(self, codons: Optional[Sequence[str]] = None, ax: Optional[Axes] = None,
                 labels: bool = False, file_name: Optional[str | Path] = None) -> Axes:
        """27. PR2 grafiği: x = G3 / (G3 + C3), y = A3 / (A3 + U3); merkez (0,5; 0,5)."""
        ax = _axes(ax)
        ax.axhline(0.5, **PLOT_GUIDE_STYLE)
        ax.axvline(0.5, **PLOT_GUIDE_STYLE)
        points = [(p["x"], p["y"], p["gene"]) for p in self.pr2_points(codons)]
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="G3 / (G3 + C3)", ylabel="A3 / (A3 + U3)",
               title="PR2 grafiği")
        ax.set_aspect("equal")
        return _finish(ax, points, labels, file_name)

    def __repr__(self) -> str:
        return f"EffectiveCodons({self.cub!r})"
