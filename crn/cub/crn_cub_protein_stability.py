"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, F. Kompozisyon temelli stabilite ve çözünürlük indeksleri.

Bir ``CUB`` nesnesindeki proteinlerin bileşiminden termostabilite ve
çözünürlük indekslerini hesaplar. Oranlar E-A'daki ``ProteinMass``'tan, NCPR
E-C'deki ``ProteinCharge``'tan gelir.

E-F. Kompozisyon temelli stabilite ve çözünürlük indeksleri:

====  ==============================================================  =================================
 No   Ölçüm                                                           Metot
====  ==============================================================  =================================
 28   IVYWREL oranı: (I + V + Y + W + R + E + L) / N                  ``ivywrel``
 29   CvP yanlılığı: [(D + E + K + R) − (N + Q + S + T)] / N × 100    ``cvp_bias``
 30   Wilkinson–Harrison: CV = λ1·(N + G + P + S)/N                   ``solubility_cv``,
      + λ2·|NCPR − 0,03|; λ1 = 15,43, λ2 = −29,56, CV' = 1,71         ``solubility_probability``
 31   Shannon entropisi: H = −Σ p_i log2 p_i (en çok log2 20)         ``composition_entropy``
====  ==============================================================  =================================

Kalıntı grupları ``crn.io.crn_io_const`` içindeki ``IVYWREL_RESIDUES``,
``CVP_CHARGED``, ``CVP_POLAR`` ve ``TURN_FORMING_RESIDUES`` sabitlerinden
gelir.

**Wilkinson–Harrison.** Sabitler, Davis ve ark. (1999, Biotechnol. Bioeng.
65:382) tarafından yeniden kalibre edilmiş iki parametreli modelinkilerdir.
Model E. coli'de aşırı ifade edilen proteinin çözünür kalıp kalmayacağını
tahmin eder: CV − CV' > 0 ise çözünmez (inklüzyon cismi), < 0 ise
çözünür tahmin edilir; olasılık
0,4934 + 0,276·|CV − CV'| − 0,0392·(CV − CV')² ile bulunur.

Kullanım::

    from crn.cub import CUB, ProteinStability

    stab = ProteinStability(cub)
    stab.ivywrel(gene="MEFV_NM_000243.3")
    stab.solubility_probability(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

import math
from typing import (Dict, Optional)

from crn.io.crn_io_const import (IVYWREL_RESIDUES, CVP_CHARGED, CVP_POLAR, TURN_FORMING_RESIDUES)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass
from crn.cub.crn_cub_protein_charge import ProteinCharge


# Wilkinson–Harrison modelinin sabitleri (Davis ve ark. 1999).
WH_LAMBDA_TURN = 15.43
WH_LAMBDA_CHARGE = -29.56
WH_CHARGE_OFFSET = 0.03
WH_DISCRIMINANT = 1.71
WH_PROBABILITY = (0.4934, 0.276, -0.0392)


class ProteinStability:
    """Bir ``CUB`` nesnesindeki proteinlerin bileşim temelli stabilite ve çözünürlük indeksleri.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        self.__mass = ProteinMass(cub)
        self.__charge = ProteinCharge(cub)

    @property
    def cub(self) -> CUB:
        return self.__mass.cub

    # --- E-F. Stabilite ve çözünürlük ------------------------------------------------

    def ivywrel(self, gene: Optional[str] = None, residues: str = IVYWREL_RESIDUES) -> float:
        """28. IVYWREL oranı: (I + V + Y + W + R + E + L) / N."""
        return self.__mass.group_composition(gene, {"ivywrel": residues})["ivywrel"]

    def cvp_bias(self, gene: Optional[str] = None,
                 charged: str = CVP_CHARGED, polar: str = CVP_POLAR) -> float:
        """29. CvP yanlılığı (%): [(D + E + K + R) − (N + Q + S + T)] / N × 100."""
        groups = self.__mass.group_composition(gene, {"charged": charged, "polar": polar})
        return (groups["charged"] - groups["polar"]) * 100

    def solubility_cv(self, gene: Optional[str] = None,
                      turn: str = TURN_FORMING_RESIDUES,
                      lambda_turn: float = WH_LAMBDA_TURN,
                      lambda_charge: float = WH_LAMBDA_CHARGE,
                      offset: float = WH_CHARGE_OFFSET) -> float:
        """30. Wilkinson–Harrison kanonik değişkeni: λ1·(N + G + P + S)/N + λ2·|NCPR − 0,03|."""
        fraction = self.__mass.group_composition(gene, {"turn": turn})["turn"]
        return lambda_turn * fraction + lambda_charge * abs(self.__charge.ncpr(gene) - offset)

    def solubility_probability(self, gene: Optional[str] = None,
                               discriminant: float = WH_DISCRIMINANT) -> Dict[str, float | bool]:
        """30. Wilkinson–Harrison tahmini: çözünür mü ve tahminin olasılığı.

        Returns:
            ``{"cv": CV, "soluble": CV − CV' <= 0, "probability": olasılık}``;
            olasılık tahmin edilen sınıfın (çözünür ya da çözünmez) olasılığıdır.
        """
        cv = self.solubility_cv(gene)
        d = cv - discriminant
        a, b, c = WH_PROBABILITY
        return {"cv": cv, "soluble": d <= 0, "probability": a + b * abs(d) + c * d * d}

    def composition_entropy(self, gene: Optional[str] = None) -> float:
        """31. Bileşimin Shannon entropisi (bit): H = −Σ p_i log2 p_i; en çok log2 20 ≈ 4,32."""
        return -sum(p / 100 * math.log2(p / 100)
                    for p in self.__mass.composition(gene).values() if p > 0)

    def __repr__(self) -> str:
        return f"ProteinStability({self.cub!r})"
