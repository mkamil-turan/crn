"""CUB grafiklerinde paylaşılan matplotlib yardımcıları ve biçimler.

Biçim sabitleri (renkler, nokta/çizgi stilleri) burada tanımlıdır; böylece
tüm ``plot_*`` metotları aynı görünümü paylaşır ve tek kaynaktan değişir.
matplotlib yalnızca bir çizim istendiğinde (``axes`` çağrılınca) yüklenir.

Gen başına nokta ya da değer çizen grafikler ``groups`` alır; her grup ayrı
renkte çizilir ve gösterge (legend) grafiğin sağına konur. ``groups``:

* metin: genin features'ındaki özellik adı, ör. ``"organism"``
* sözlük: gen adı -> grup
* fonksiyon: ``Gene`` alır, grup adını döner
"""
from __future__ import annotations

from pathlib import Path
from typing import (TYPE_CHECKING, Callable, Dict, List, Mapping, Optional, Sequence,
                    Tuple, Union)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from crn.cub.core import CUB
    from crn.gene import Gene

__all__ = ["PLOT_POINT_STYLE", "PLOT_CURVE_STYLE", "PLOT_GUIDE_STYLE", "PLOT_LABEL_SIZE",
           "PLOT_GENE_LABEL_SIZE", "PLOT_GROUP_COLORS", "Groups", "axes", "save",
           "gene_symbols", "gene_groups", "group_colors", "scatter", "histogram",
           "group_legend"]

PLOT_POINT_STYLE = {"s": 18, "color": "#2a6fb0", "edgecolor": "none", "alpha": 0.8, "zorder": 3}
PLOT_CURVE_STYLE = {"color": "#c0392b", "linewidth": 1.4}
PLOT_GUIDE_STYLE = {"color": "#7f7f7f", "linestyle": "--", "linewidth": 0.8}
PLOT_LABEL_SIZE = 7          # eksen yazıları (kodonlar, indeks adları)
PLOT_GENE_LABEL_SIZE = 4     # noktaların yanındaki gen sembolleri
PLOT_LEGEND_SIZE = 7
# Ayırt edilebilir renkler (Okabe-Ito + tab10); grup sayısı aşarsa başa döner.
PLOT_GROUP_COLORS = ["#0072b2", "#e69f00", "#009e73", "#d55e00", "#cc79a7", "#56b4e9",
                     "#882255", "#000000", "#8c564b", "#7f7f7f", "#bcbd22", "#17becf",
                     "#9467bd", "#ff9896"]

Groups = Union[None, str, Mapping[str, str], Callable[["Gene"], str]]


def axes(ax: Optional[Axes], figsize: Tuple[float, float] = (5, 4.5)) -> Axes:
    """Verilen ``Axes``'i, yoksa yeni bir şeklin ``Axes``'ini döner (matplotlib burada yüklenir)."""
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(figsize=figsize)
    return ax


def save(ax: Axes, file_name: Optional[str | Path]) -> Axes:
    """``file_name`` verilmişse şekli 300 dpi kaydeder ve ``Axes``'i döner."""
    if file_name is not None:
        ax.figure.savefig(file_name, dpi=300, bbox_inches="tight")
    return ax


# --- Gen sembolleri ve gruplar -------------------------------------------------------

def gene_symbols(cub: CUB, names: Sequence[str]) -> List[str]:
    """Genlerin sembolleri (features'taki ``gene_name``; yoksa genin adı)."""
    genes = {gene.name: gene for gene in cub.genes}
    return [genes[n].features.first("gene_name") or n for n in names]


def gene_groups(cub: CUB, names: Sequence[str], groups: Groups) -> Optional[List[str]]:
    """Her genin grubu; ``groups`` None ise None.

    Raises:
        ValueError: Bir genin grubu bulunamazsa.
    """
    if groups is None:
        return None
    genes = {gene.name: gene for gene in cub.genes}
    result = []
    for name in names:
        if isinstance(groups, str):
            value = genes[name].features.first(groups)
        elif isinstance(groups, Mapping):
            value = groups.get(name)
        else:
            value = groups(genes[name])
        if value in (None, ""):
            raise ValueError(f"{name}: grubu bulunamadı (groups={groups!r})")
        result.append(str(value))
    return result


def group_colors(groups: Sequence[str]) -> Dict[str, str]:
    """Gruplara ilk görünüş sırasıyla renk atar."""
    order = list(dict.fromkeys(groups))
    return {g: PLOT_GROUP_COLORS[i % len(PLOT_GROUP_COLORS)] for i, g in enumerate(order)}


def group_legend(ax: Axes, handles: Sequence, title: Optional[str] = None) -> None:
    """Grup göstergesini grafiğin sağına koyar; var olan göstergeyi korur."""
    existing = ax.get_legend()
    if existing is not None:
        ax.add_artist(existing)
    ax.legend(handles=handles, title=title, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, fontsize=PLOT_LEGEND_SIZE, title_fontsize=PLOT_LEGEND_SIZE,
              borderaxespad=0.0)


def scatter(ax: Axes, cub: CUB, points: Sequence[Tuple[float, float, str]],
            groups: Groups = None, labels: bool = False) -> None:
    """Gen noktalarını (x, y, gen adı) çizer; gruplar renkli, etiketler sembol."""
    names = [name for _, _, name in points]
    assigned = gene_groups(cub, names, groups)
    if assigned is None:
        ax.scatter([x for x, _, _ in points], [y for _, y, _ in points], **PLOT_POINT_STYLE)
    else:
        handles = []
        for group, color in group_colors(assigned).items():
            members = [p for p, g in zip(points, assigned) if g == group]
            handles.append(ax.scatter([x for x, _, _ in members], [y for _, y, _ in members],
                                      label=group, **{**PLOT_POINT_STYLE, "color": color}))
        group_legend(ax, handles)
    if labels:
        for (x, y, _), symbol in zip(points, gene_symbols(cub, names)):
            ax.annotate(symbol, (x, y), xytext=(2, 2), textcoords="offset points",
                        fontsize=PLOT_GENE_LABEL_SIZE, color="#333333")


def histogram(ax: Axes, cub: CUB, values: Sequence[Tuple[str, float]], bins,
              groups: Groups = None, **kwargs) -> None:
    """(gen adı, değer) çiftlerinin histogramı; gruplar varsa renkli ve yığılmış."""
    assigned = gene_groups(cub, [name for name, _ in values], groups)
    if assigned is None:
        ax.hist([v for _, v in values], bins=bins, color=PLOT_POINT_STYLE["color"],
                edgecolor="white", **kwargs)
        return
    colors = group_colors(assigned)
    data = [[v for (_, v), g in zip(values, assigned) if g == group] for group in colors]
    ax.hist(data, bins=bins, stacked=True, color=list(colors.values()), label=list(colors),
            edgecolor="white", linewidth=0.4, **kwargs)
    import matplotlib.patches as mpatches
    group_legend(ax, [mpatches.Patch(color=c, label=g) for g, c in colors.items()])
