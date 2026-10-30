"""Tüm genler için skaler CUB ve protein ölçümlerini tek tabloda toplar."""
from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Mapping, Optional, Sequence, Tuple

from crn.cub.core import CUB, CubWarning
from crn.cub.composition import Composition
from crn.cub.adaptation import Adaptation
from crn.cub.enc import EffectiveCodons
from crn.cub.divergence import Divergence
from crn.cub.context import CodonContext
from crn.protein import ProteinProperties

if TYPE_CHECKING:
    from matplotlib.axes import Axes

# Korelasyon ısı haritasında varsayılan olarak denenecek CUB indeksleri (varsa alınır).
DEFAULT_CORRELATION_COLUMNS = ("L", "GC%", "GC3", "GC3s", "ENC", "ENC_ratio",
                               "SCUO", "CDC", "RCBS", "Fop", "CAI", "CBI", "tAI", "CPB", "RCDI")


def summary_table(cub: CUB, reference: Optional[CUB] = None, background: Optional[CUB] = None,
                  host: Optional[CUB] = None, genome: Optional[CUB] = None,
                  trna_counts: Optional[Mapping[str, float]] = None, kingdom: str = "eukaryote",
                  protein: bool = True, file_name: Optional[str | Path] = None
                  ) -> List[Dict[str, object]]:
    """Her gen için bir satır: ``{"gene": ad, ölçüm: değer, ...}``.

    Her zaman: L, GC%, GC1, GC2, GC12, GC3, GC3s, ENC, ENC_expected, ENC_ratio,
    SCUO, CDC, RCBS. İsteğe bağlı sütunlar ilgili veri verilince eklenir:

    - ``reference`` (yüksek ifadeli genler): Fop, CAI, CBI; ``background``
      verilirse Fop_chi2 ve CBI_chi2 de (optimal kodonlar χ² ile).
    - ``host``: RCDI.  ``genome``: CPB.  ``trna_counts``: tAI.
    - ``protein`` True ise ``ProteinProperties.summary`` sütunları.

    Satır satır ``CubWarning`` uyarıları (RCDI/CPB'de referansta olmayan
    kodon/çiftler) bastırılır. ``file_name`` verilirse tablo sekmeyle ayrılmış
    (TSV) olarak da yazılır.
    """
    comp, eff, div = Composition(cub), EffectiveCodons(cub), Divergence(cub)
    ada = Adaptation(cub, reference, background) if reference is not None else None
    rcdi = Divergence(cub, host) if host is not None else None
    cpb = CodonContext(cub, genome) if genome is not None else None
    tai = CodonContext(cub) if trna_counts else None
    props = ProteinProperties(cub) if protein else None

    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CubWarning)
        for name in cub.names:
            gc3s = comp.gc3s(name)
            enc = eff.enc(name)
            expected = eff.expected_enc(gc3s)
            row: Dict[str, object] = {"gene": name, "L": cub.total_codons(name)}
            row.update({key: value for key, value in comp.summary(name).items()
                        if key in ("GC%", "GC1", "GC2", "GC12", "GC3", "GC3s")})
            row.update({"ENC": enc, "ENC_expected": expected,
                        "ENC_ratio": (expected - enc) / expected,
                        "SCUO": div.scuo(name), "CDC": div.cdc(name), "RCBS": div.rcbs(name)})
            if ada is not None:
                row.update({"Fop": ada.synonymous.fop(name), "CAI": ada.cai(name), "CBI": ada.cbi(name)})
                if background is not None:
                    row.update({"Fop_chi2": ada.synonymous.fop(name, method="chi2"),
                                "CBI_chi2": ada.cbi(name, method="chi2")})
            if rcdi is not None:
                row["RCDI"] = rcdi.rcdi(name)
            if cpb is not None:
                row["CPB"] = cpb.cpb(name)
            if tai is not None:
                row["tAI"] = tai.tai(name, trna_counts, kingdom=kingdom)
            if props is not None:
                row.update(props.summary(name))
            rows.append(row)

    if file_name is not None:
        with open(file_name, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    return rows


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    """İki dizinin Pearson korelasyonu; değişkenlik yoksa ``nan``."""
    n = len(xs)
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if not sxx or not syy:
        return math.nan
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy)


def correlation_matrix(rows: List[Dict[str, object]],
                       columns: Optional[Sequence[str]] = None
                       ) -> Tuple[List[str], List[List[float]]]:
    """``summary_table`` satırlarından seçili sütunların Pearson korelasyon matrisi.

    ``columns`` verilmezse ``DEFAULT_CORRELATION_COLUMNS``'tan satırlarda
    bulunan, sayısal ve değişkenliği olan sütunlar alınır. ``nan`` içeren
    değerler (ör. bazı genlerde ENC) o çiftte atlanır.

    Returns:
        ``(columns, matrix)``; ``matrix[i][j]`` i. ve j. sütunun korelasyonu.
    """
    candidates = list(columns) if columns is not None else [
        c for c in DEFAULT_CORRELATION_COLUMNS if rows and c in rows[0]]
    series: Dict[str, List[float]] = {}
    for col in candidates:
        values = [row.get(col) for row in rows]
        numeric = [v for v in values if isinstance(v, (int, float)) and not (
            isinstance(v, float) and math.isnan(v))]
        if len(set(numeric)) > 1:                      # değişkenliği olan sayısal sütun
            series[col] = values
    names = list(series)
    matrix = [[1.0] * len(names) for _ in names]
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            pairs = [(x, y) for x, y in zip(series[a], series[b])
                     if isinstance(x, (int, float)) and isinstance(y, (int, float))
                     and not (math.isnan(x) or math.isnan(y))]
            r = _pearson([p[0] for p in pairs], [p[1] for p in pairs]) if len(pairs) > 1 else math.nan
            matrix[i][j] = matrix[j][i] = r
    return names, matrix


def plot_correlations(cub: CUB, columns: Optional[Sequence[str]] = None,
                      ax: Optional["Axes"] = None, annotate: bool = True,
                      file_name: Optional[str | Path] = None, **summary_kwargs
                      ) -> "Axes":
    """CUB indeksleri arasındaki Pearson korelasyonlarının ısı haritası.

    ``summary_table`` ile gen başına indeksleri hesaplar, aralarındaki
    korelasyonu −1…1 arası renklerle çizer. Hangi sütunların hesaplanacağı
    ``summary_table``'ın isteğe bağlı argümanlarına bağlıdır (``reference``
    → Fop/CAI/CBI, ``trna_counts`` → tAI, ...); bunlar ``summary_kwargs`` ile
    geçirilir. ``columns`` ile sütunlar elle de seçilebilir.

    Örn. ``plot_correlations(cub, reference=high)`` CAI'yi de ekler.
    """
    from crn.cub._plotting import PLOT_LABEL_SIZE, axes, save
    rows = summary_table(cub, **summary_kwargs)
    names, matrix = correlation_matrix(rows, columns)
    if not names:
        raise ValueError("Korelasyon için değişkenliği olan sayısal sütun yok")

    ax = axes(ax, figsize=(max(5, len(names) * 0.6), max(4, len(names) * 0.6)))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-1.0, vmax=1.0, interpolation="nearest")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=90, fontsize=PLOT_LABEL_SIZE)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=PLOT_LABEL_SIZE)
    if annotate:
        for i in range(len(names)):
            for j in range(len(names)):
                r = matrix[i][j]
                if not math.isnan(r):
                    ax.text(j, i, f"{r:.2f}", ha="center", va="center", fontsize=6,
                            color="white" if abs(r) > 0.55 else "#222222")
    ax.figure.colorbar(image, ax=ax, label="Pearson r", fraction=0.045, pad=0.03)
    ax.set_title("CUB indeksleri korelasyonu")
    return save(ax, file_name)
