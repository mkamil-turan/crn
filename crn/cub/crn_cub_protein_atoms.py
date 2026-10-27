"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, B. Atom bileşimi.

Bir ``CUB`` nesnesinin amino asit sayımlarından (``CUB.amino_acid_counts``)
proteinin element formülünü hesaplar. Uzunluk ve zincir sayısı E-A'daki
``ProteinMass``'tan gelir.

E-B. Atom bileşimi:

====  ==============================================================  =====================
 No   Ölçüm                                                           Metot
====  ==============================================================  =====================
 8    Atom formülü: Σ n_i × amino asit formülü − (N − 1) H2O          ``formula``,
                                                                      ``formula_string``
 9    Toplam atom sayısı: C + H + N + O + S                           ``atom_count``
 10   Kalıntı başına atom sayısı: toplam atom / N                     ``atoms_per_residue``
 11   Ortalama karbon oksidasyon durumu (Dick): Z_C =                 ``carbon_oxidation_state``
      (Z − nH + 3nN + 2nO + 2nS) / nC
====  ==============================================================  =====================

Formüller ``crn.io.crn_io_const`` içindeki ``AMINO_ACID_FORMULAS`` (serbest
amino asitler, ``ELEMENTS`` sırasıyla) ve ``WATER_FORMULA`` sabitlerinden
gelir. ``gene`` verilmezse tüm genler birlikte hesaplanır; her protein ayrı
bir zincir sayılır, yani her biri için (N − 1) su çıkarılır.

Z_C'de Z proteinin net yüküdür; nötr protein için 0'dır. Değer −4 (CH4) ile
+4 (CO2) arasındadır; proteinlerde genellikle −0,4 ile 0 arasında çıkar.

Kullanım::

    from crn.cub import CUB, ProteinAtoms

    atoms = ProteinAtoms(cub)
    atoms.formula_string(gene="MEFV_NM_000243.3")      # "C3779H6038N1104O1167S26"
    atoms.carbon_oxidation_state(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

from typing import (Dict, Mapping, Optional, Sequence, Tuple)

from crn.io.crn_io_const import (ELEMENTS, AMINO_ACID_FORMULAS, WATER_FORMULA)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass


class ProteinAtoms:
    """Bir ``CUB`` nesnesindeki proteinlerin atom bileşimi.

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

    # --- E-B. Atom bileşimi -------------------------------------------------------

    def formula(self, gene: Optional[str] = None,
                formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                water: Tuple[int, ...] = WATER_FORMULA,
                elements: Sequence[str] = ELEMENTS) -> Dict[str, int]:
        """8. Atom formülü: Σ(n_i × amino asit formülü) − (N − zincir sayısı) × H2O.

        Args:
            formulas: Amino asit -> ``elements`` sırasıyla atom sayıları.
            water: Peptit bağında ayrılan suyun atom sayıları.
            elements: Formül tuple'larındaki element sırası.

        Returns:
            Element -> atom sayısı, ``elements`` sırasıyla; protein boşsa
            hepsi 0.

        Raises:
            KeyError: Sayılan bir amino asidin formülü ``formulas`` içinde yoksa.
        """
        counts = self.cub.amino_acid_counts(gene, tuple(formulas))
        n = sum(counts.values())
        bonds = n - self.__mass.chain_count(gene) if n else 0
        return {element: sum(count * formulas[aa][i] for aa, count in counts.items())
                         - bonds * water[i]
                for i, element in enumerate(elements)}

    def formula_string(self, gene: Optional[str] = None,
                       formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                       water: Tuple[int, ...] = WATER_FORMULA,
                       elements: Sequence[str] = ELEMENTS) -> str:
        """8. Atom formülü metin olarak, ör. ``"C3779H6038N1104O1167S26"`` (0 olan elementler yazılmaz)."""
        return "".join(f"{element}{n}" for element, n
                       in self.formula(gene, formulas, water, elements).items() if n)

    def atom_count(self, gene: Optional[str] = None,
                   formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                   water: Tuple[int, ...] = WATER_FORMULA,
                   elements: Sequence[str] = ELEMENTS) -> int:
        """9. Toplam atom sayısı: formüldeki tüm atomların toplamı."""
        return sum(self.formula(gene, formulas, water, elements).values())

    def atoms_per_residue(self, gene: Optional[str] = None,
                          formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                          water: Tuple[int, ...] = WATER_FORMULA,
                          elements: Sequence[str] = ELEMENTS) -> float:
        """10. Kalıntı başına atom sayısı: toplam atom / N."""
        n = self.__mass.length(gene, tuple(formulas))
        return self.atom_count(gene, formulas, water, elements) / n if n else 0.0

    def carbon_oxidation_state(self, gene: Optional[str] = None,
                               charge: float = 0.0,
                               formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                               water: Tuple[int, ...] = WATER_FORMULA,
                               elements: Sequence[str] = ELEMENTS) -> float:
        """11. Ortalama nominal karbon oksidasyon durumu (Dick): (Z − nH + 3nN + 2nO + 2nS) / nC.

        Args:
            charge: Proteinin net yükü Z; nötr protein için 0.

        Returns:
            Z_C; formülde karbon yoksa 0.
        """
        f = self.formula(gene, formulas, water, elements)
        if not f.get("C"):
            return 0.0
        return (charge - f.get("H", 0) + 3 * f.get("N", 0) + 2 * f.get("O", 0)
                + 2 * f.get("S", 0)) / f["C"]

    def __repr__(self) -> str:
        return f"ProteinAtoms({self.cub!r})"
