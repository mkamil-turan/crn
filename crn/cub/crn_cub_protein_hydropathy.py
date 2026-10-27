"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, D. Hidrofobiklik ve polarite.

Bir ``CUB`` nesnesindeki proteinlerin hidrofobikliğini hesaplar. Ölçek
ortalamaları E-A'daki ``ProteinMass.scale_average``, oranlar
``ProteinMass.composition`` ve ``group_composition`` ile bulunur.

E-D. Hidrofobiklik ve polarite:

====  ==============================================================  ======================
 No   Ölçüm                                                           Metot
====  ==============================================================  ======================
 19   GRAVY: Σ(n_i × KD_i) / N, Kyte–Doolittle                        ``gravy``
 20   Diğer ölçeklerle ortalama: Σ(n_i × ölçek_i) / N                 ``hydropathy``
 21   Aromo (aromatiklik): (F + Y + W) / N                            ``aromaticity``
 22   Alifatik indeks (Ikai): X(A) + 2,9·X(V) + 3,9·[X(I) + X(L)]     ``aliphatic_index``
 23   Boman indeksi: Σ(n_i × ΔG_i) / N                                ``boman_index``
====  ==============================================================  ======================

Ölçekler ``crn.io.crn_io_const.HYDROPATHY_SCALES`` (``"Kyte-Doolittle"``,
``"Hopp-Woods"``, ``"Eisenberg"``, ``"Wimley-White"``, ``"Janin"``) ve
``BOMAN_SCALE``'den gelir; aynı biçimde bir sözlük de verilebilir.
Alifatik indekste X mol yüzdesidir (E-A ``composition``). Boman indeksi
kcal/mol'dür; yaklaşık 2,48'in üstü yüksek protein bağlama potansiyeli
gösterir.

Kullanım::

    from crn.cub import CUB, ProteinHydropathy

    hyd = ProteinHydropathy(cub)
    hyd.gravy(gene="MEFV_NM_000243.3")
    hyd.hydropathy(gene="MEFV_NM_000243.3", scale="Hopp-Woods")
"""
from __future__ import annotations

from typing import (Mapping, Optional)

from crn.io.crn_io_const import (HYDROPATHY_SCALES, BOMAN_SCALE, AMINO_ACID_GROUPS)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass


# Ikai (1980) alifatik indeksindeki Val ve Ile/Leu katsayıları (Ala'ya göre hacim oranları).
ALIPHATIC_V = 2.9
ALIPHATIC_IL = 3.9


class ProteinHydropathy:
    """Bir ``CUB`` nesnesindeki proteinlerin hidrofobikliği ve polaritesi.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        self.__mass = ProteinMass(cub)

    @property
    def cub(self) -> CUB:
        return self.__mass.cub

    @staticmethod
    def _scale(scale: str | Mapping[str, float]) -> Mapping[str, float]:
        if isinstance(scale, str):
            if scale not in HYDROPATHY_SCALES:
                raise ValueError(f"scale {', '.join(map(repr, HYDROPATHY_SCALES))} ya da bir "
                                 f"sözlük olmalı, {scale!r} verildi")
            return HYDROPATHY_SCALES[scale]
        return scale

    # --- E-D. Hidrofobiklik ve polarite ----------------------------------------------

    def gravy(self, gene: Optional[str] = None) -> float:
        """19. GRAVY: Kyte–Doolittle değerlerinin ortalaması."""
        return self.hydropathy(gene, "Kyte-Doolittle")

    def hydropathy(self, gene: Optional[str] = None,
                   scale: str | Mapping[str, float] = "Kyte-Doolittle") -> float:
        """20. Seçilen ölçeğin ortalaması: Σ(n_i × ölçek_i) / N.

        Args:
            scale: ``HYDROPATHY_SCALES`` içindeki ölçeğin adı ya da amino asit -> değer.

        Raises:
            ValueError: scale geçersiz bir ölçek adıysa.
        """
        return self.__mass.scale_average(self._scale(scale), gene)

    def aromaticity(self, gene: Optional[str] = None,
                    aromatic: str = AMINO_ACID_GROUPS["aromatic"]) -> float:
        """21. Aromo: (F + Y + W) / N."""
        return self.__mass.group_composition(gene, {"aromatic": aromatic})["aromatic"]

    def aliphatic_index(self, gene: Optional[str] = None,
                        v: float = ALIPHATIC_V, il: float = ALIPHATIC_IL) -> float:
        """22. Alifatik indeks (Ikai): X(A) + a·X(V) + b·[X(I) + X(L)], X mol yüzdesi."""
        x = self.__mass.composition(gene)
        return x["A"] + v * x["V"] + il * (x["I"] + x["L"])

    def boman_index(self, gene: Optional[str] = None,
                    scale: Mapping[str, float] = BOMAN_SCALE) -> float:
        """23. Boman indeksi (kcal/mol): Σ(n_i × ΔG_i) / N."""
        return self.__mass.scale_average(scale, gene)

    def __repr__(self) -> str:
        return f"ProteinHydropathy({self.cub!r})"
