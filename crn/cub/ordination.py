"""Kodon kullanım yanlılığı (codon usage bias) — çok değişkenli sıralama (COA).

Genlerin RSCU değerlerinden **uyum çözümlemesi** (correspondence analysis,
COA) ile genleri birkaç ana eksen üzerine yerleştirir. Literatürde kodon
kullanımının başlıca değişim ekseninin (çoğu zaman GC3 / ifade düzeyi)
bulunması ve genlerin bu eksen boyunca saçılımının gösterilmesi için
kullanılır (ör. codonW'nin COA'sı).

Yöntem standart uyum çözümlemesidir. Gen (satır) × eş anlamlı kodon (sütun)
RSCU tablosu ``X`` için::

    P = X / ΣX,   satır kütlesi r_i = Σ_j P_ij,   sütun kütlesi c_j = Σ_i P_ij
    S_ij = (P_ij − r_i c_j) / √(r_i c_j)            (standartlaştırılmış artık)

``S = UΣVᵀ`` tekil değer ayrışımının öz-çözümü ``C = SᵀS``'nin (p×p simetrik)
özdeğerlerinden alınır; bu paket numpy'a bağlı olmadığından özçözüm küçük bir
saf-Python Jacobi algoritmasıyla yapılır (``enc.py``'deki elle yazılmış
regresyonla aynı çizgide). Eksen k'nin açıkladığı atalet oranı
λ_k / Σλ'dır; genlerin asıl koordinatları ``F = D_r^{-1/2} (S V)``.

Kullanım::

    from crn.cub import CUB, CorrespondenceAnalysis

    coa = CorrespondenceAnalysis(cub)
    coa.inertia()[:3]                 # ilk üç eksenin atalet oranı
    coa.coordinates(n_axes=2)         # [{"gene", "axis1", "axis2"}, ...]
    coa.plot_coa(labels=True, file_name="coa.png")
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Sequence, Tuple

from crn.cub.core import CUB
from crn.cub.synonymous import Synonymous

if TYPE_CHECKING:
    from matplotlib.axes import Axes

JACOBI_MAX_SWEEPS = 100
JACOBI_TOLERANCE = 1e-12


def _jacobi_eigen(matrix: List[List[float]]) -> Tuple[List[float], List[List[float]]]:
    """Simetrik matrisin özdeğer/özvektörleri (döngüsel Jacobi); özdeğerler azalan.

    Returns:
        ``(eigenvalues, eigenvectors)``; ``eigenvectors[k]`` k. özvektör (sütun).
    """
    n = len(matrix)
    a = [row[:] for row in matrix]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for _ in range(JACOBI_MAX_SWEEPS):
        off = math.sqrt(sum(a[p][q] ** 2 for p in range(n) for q in range(p + 1, n)))
        if off < JACOBI_TOLERANCE:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(a[p][q]) < JACOBI_TOLERANCE:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1))
                if theta == 0:
                    t = 1.0
                c = 1 / math.sqrt(t * t + 1)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq

    eigenvalues = [a[i][i] for i in range(n)]
    columns = [[v[row][col] for row in range(n)] for col in range(n)]
    order = sorted(range(n), key=lambda k: eigenvalues[k], reverse=True)
    return [eigenvalues[k] for k in order], [columns[k] for k in order]


class CorrespondenceAnalysis:
    """Genlerin RSCU tablosunun uyum çözümlemesi (COA).

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi (en az iki gen).

    Raises:
        TypeError: cub bir CUB değilse.
        ValueError: ikiden az gen varsa ya da tüm RSCU değerleri sabitse
            (çözümlenecek değişim yok).
    """

    def __init__(self, cub: CUB):
        if not isinstance(cub, CUB):
            raise TypeError(f"cub bir CUB olmalı, {type(cub).__name__} verildi")
        if len(cub.names) < 2:
            raise ValueError("COA için en az iki gen gerekir")
        self.cub = cub
        syn = Synonymous(cub)
        all_codons, _ = syn._synonymous_layout()
        names = list(cub.names)
        table = [[syn.rscu(name)[codon] for codon in all_codons] for name in names]

        # Hiçbir gende kullanılmayan kodonları (sütun kütlesi 0) çıkar.
        keep = [j for j in range(len(all_codons)) if any(row[j] for row in table)]
        if not keep:
            raise ValueError("RSCU tablosu boş; COA hesaplanamaz")
        self.genes = names
        self.codons = [all_codons[j] for j in keep]
        matrix = [[row[j] for j in keep] for row in table]
        self._fit(matrix)

    def _fit(self, x: List[List[float]]) -> None:
        n, p = len(x), len(x[0])
        total = sum(sum(row) for row in x)
        if total <= 0:
            raise ValueError("RSCU tablosu boş; COA hesaplanamaz")
        rows_mass = [sum(row) / total for row in x]
        cols_mass = [sum(x[i][j] for i in range(n)) / total for j in range(p)]

        # Standartlaştırılmış artık S.
        s = [[0.0] * p for _ in range(n)]
        for i in range(n):
            for j in range(p):
                denom = math.sqrt(rows_mass[i] * cols_mass[j])
                s[i][j] = ((x[i][j] / total) - rows_mass[i] * cols_mass[j]) / denom if denom else 0.0

        c = [[sum(s[k][i] * s[k][j] for k in range(n)) for j in range(p)] for i in range(p)]
        eigenvalues, eigenvectors = _jacobi_eigen(c)
        eigenvalues = [max(0.0, ev) for ev in eigenvalues]

        total_inertia = sum(eigenvalues) or 1.0
        # Öz-denetim: Σλ ≈ Σ S_ij² (toplam atalet).
        self._check = sum(s[i][j] ** 2 for i in range(n) for j in range(p))

        rank = min(n - 1, p - 1)
        self.eigenvalues = eigenvalues[:rank]
        self._inertia = [ev / total_inertia for ev in self.eigenvalues]

        # Satır (gen) asıl koordinatları: F = D_r^{-1/2} (S V).
        coords: List[List[float]] = []
        for i in range(n):
            row = []
            for k in range(rank):
                vk = eigenvectors[k]
                sv = sum(s[i][j] * vk[j] for j in range(p))
                row.append(sv / math.sqrt(rows_mass[i]) if rows_mass[i] else 0.0)
            coords.append(row)
        self._coords = coords

    def inertia(self) -> List[float]:
        """Her eksenin açıkladığı atalet oranı (azalan; toplamı 1)."""
        return list(self._inertia)

    def coordinates(self, n_axes: int = 2) -> List[Dict[str, float | str]]:
        """Her gen için ilk ``n_axes`` eksendeki asıl koordinat.

        Returns:
            ``[{"gene": ad, "axis1": .., "axis2": ..}, ...]``.
        """
        n_axes = max(1, min(n_axes, len(self.eigenvalues)))
        result = []
        for name, row in zip(self.genes, self._coords):
            point: Dict[str, float | str] = {"gene": name}
            for k in range(n_axes):
                point[f"axis{k + 1}"] = row[k]
            result.append(point)
        return result

    def plot_coa(self, x_axis: int = 1, y_axis: int = 2, ax: Optional[Axes] = None,
                 labels: bool = False, file_name: Optional[str | Path] = None) -> Axes:
        """COA saçılım grafiği: genler iki eksen üzerinde; atalet oranları eksende.

        Args:
            x_axis, y_axis: Çizilecek eksenler (1 tabanlı).
        """
        from crn.cub._plotting import PLOT_GUIDE_STYLE, PLOT_LABEL_SIZE, PLOT_POINT_STYLE, axes, save
        if max(x_axis, y_axis) > len(self.eigenvalues):
            raise ValueError(f"yalnızca {len(self.eigenvalues)} eksen var")
        ax = axes(ax)
        xi, yi = x_axis - 1, y_axis - 1
        points = [(row[xi], row[yi], name) for row, name in zip(self._coords, self.genes)]
        ax.axhline(0.0, **PLOT_GUIDE_STYLE)
        ax.axvline(0.0, **PLOT_GUIDE_STYLE)
        ax.scatter([x for x, _, _ in points], [y for _, y, _ in points], **PLOT_POINT_STYLE)
        if labels:
            for x, y, name in points:
                ax.annotate(name, (x, y), xytext=(3, 3), textcoords="offset points",
                            fontsize=PLOT_LABEL_SIZE)
        ax.set(xlabel=f"Eksen {x_axis} (%{self._inertia[xi] * 100:.1f})",
               ylabel=f"Eksen {y_axis} (%{self._inertia[yi] * 100:.1f})",
               title="Uyum çözümlemesi (COA)")
        return save(ax, file_name)

    def __repr__(self) -> str:
        return (f"CorrespondenceAnalysis({self.cub!r}): {len(self.genes)} gen, "
                f"{len(self.codons)} kodon, {len(self.eigenvalues)} eksen")
