"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, E. Optik ve hidrodinamik özellikler.

Bir ``CUB`` nesnesindeki proteinlerin absorbansını, kısmi spesifik hacmini ve
ortalama kalıntı hacmini hesaplar. Sayımlar, kütle, disülfit eşleri ve ölçek
ortalamaları E-A'daki ``ProteinMass``'tan gelir.

E-E. Optik ve hidrodinamik özellikler:

====  ==============================================================  ===============================
 No   Ölçüm                                                           Metot
====  ==============================================================  ===============================
 24   ε280: 5500·nTrp + 1490·nTyr + 125·nCystine (M⁻¹ cm⁻¹)           ``extinction_coefficient``
 25   %0,1 absorbans (1 mg/mL): ε / MW                                ``absorbance_01``
 26   Kısmi spesifik hacim (Cohn–Edsall): Σ n·M·v / Σ n·M (cm³/g)     ``partial_specific_volume``
 27   Ortalama kalıntı hacmi (Å³) / yüzey alanı (Å²): Σ n·x / N       ``mean_residue_volume``,
                                                                      ``mean_residue_surface``
====  ==============================================================  ===============================

Katsayılar ve ölçekler ``crn.io.crn_io_const`` içindeki
``EXTINCTION_COEFFICIENTS_280`` (Pace ve ark. 1995),
``RESIDUE_SPECIFIC_VOLUMES`` (Cohn ve Edsall), ``RESIDUE_VOLUMES``
(Zamyatnin) ve ``RESIDUE_MAX_ASA`` (Tien ve ark. 2013) sabitlerinden gelir.

ε'da cystine sayısı tüm Cys'lerin eşleştiği varsayımıyla ⌊nCys / 2⌋'dir
(``ProteinMass.disulfide_pairs``); ``reduced=True`` ise Cys'ler serbest
sayılır ve bu terim 0 olur. ExPASy ProtParam iki değeri de verir.

Kullanım::

    from crn.cub import CUB, ProteinOptics

    opt = ProteinOptics(cub)
    opt.extinction_coefficient(gene="MEFV_NM_000243.3")
    opt.absorbance_01(gene="MEFV_NM_000243.3", reduced=True)
"""
from __future__ import annotations

from typing import (Mapping, Optional)

from crn.io.crn_io_const import (EXTINCTION_COEFFICIENTS_280, RESIDUE_SPECIFIC_VOLUMES,
                                 RESIDUE_AVERAGE_MASS, RESIDUE_VOLUMES, RESIDUE_MAX_ASA)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass


class ProteinOptics:
    """Bir ``CUB`` nesnesindeki proteinlerin optik ve hidrodinamik özellikleri.

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

    # --- E-E. Optik ve hidrodinamik özellikler -----------------------------------------

    def extinction_coefficient(self, gene: Optional[str] = None, reduced: bool = False,
                               coefficients: Mapping[str, float] = EXTINCTION_COEFFICIENTS_280
                               ) -> float:
        """24. 280 nm molar ekstinksiyon katsayısı (M⁻¹ cm⁻¹).

        Args:
            reduced: True ise tüm Cys'ler serbest (cystine terimi 0); False ise
                tüm Cys'ler eşleşmiş varsayılır.
            coefficients: ``"W"``, ``"Y"`` ve ``"cystine"`` katkıları.
        """
        counts = self.__mass.residue_counts(gene, "WY")
        cystines = 0 if reduced else self.__mass.disulfide_pairs(gene)
        return (coefficients["W"] * counts["W"] + coefficients["Y"] * counts["Y"]
                + coefficients["cystine"] * cystines)

    def absorbance_01(self, gene: Optional[str] = None, reduced: bool = False,
                      coefficients: Mapping[str, float] = EXTINCTION_COEFFICIENTS_280) -> float:
        """25. %0,1 (1 mg/mL) çözeltinin 280 nm absorbansı: ε / MW."""
        mw = self.__mass.molecular_weight(gene)
        return self.extinction_coefficient(gene, reduced, coefficients) / mw if mw else 0.0

    def partial_specific_volume(self, gene: Optional[str] = None,
                                volumes: Mapping[str, float] = RESIDUE_SPECIFIC_VOLUMES,
                                masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS) -> float:
        """26. Kısmi spesifik hacim v̄ (cm³/g), Cohn–Edsall: Σ(n_i·M_i·v_i) / Σ(n_i·M_i)."""
        weights = {aa: masses[aa] for aa in volumes}
        weighted = {aa: masses[aa] * volumes[aa] for aa in volumes}
        total = self.__mass.scale_average(weights, gene)
        return self.__mass.scale_average(weighted, gene) / total if total else 0.0

    def mean_residue_volume(self, gene: Optional[str] = None,
                            volumes: Mapping[str, float] = RESIDUE_VOLUMES) -> float:
        """27. Ortalama kalıntı hacmi (Å³, Zamyatnin): Σ(n_i × hacim_i) / N."""
        return self.__mass.scale_average(volumes, gene)

    def mean_residue_surface(self, gene: Optional[str] = None,
                             surfaces: Mapping[str, float] = RESIDUE_MAX_ASA) -> float:
        """27. Ortalama kalıntı yüzey alanı (Å², en büyük ASA): Σ(n_i × ASA_i) / N."""
        return self.__mass.scale_average(surfaces, gene)

    def __repr__(self) -> str:
        return f"ProteinOptics({self.cub!r})"
