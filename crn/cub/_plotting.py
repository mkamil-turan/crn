"""CUB grafiklerinde paylaşılan matplotlib yardımcıları ve biçimler.

Biçim sabitleri (renkler, nokta/çizgi stilleri) ``enc.py``'den alınır; böylece
tüm ``plot_*`` metotları aynı görünümü paylaşır ve tek kaynaktan değişir.
matplotlib yalnızca bir çizim istendiğinde (``axes`` çağrılınca) yüklenir.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional, Tuple

from crn.cub.enc import (PLOT_POINT_STYLE, PLOT_CURVE_STYLE, PLOT_GUIDE_STYLE,
                         PLOT_LABEL_SIZE)

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = ["PLOT_POINT_STYLE", "PLOT_CURVE_STYLE", "PLOT_GUIDE_STYLE",
           "PLOT_LABEL_SIZE", "axes", "save"]


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
