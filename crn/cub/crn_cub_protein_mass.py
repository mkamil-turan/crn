"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, A. Sayım ve kütle.

Bir ``CUB`` nesnesinin amino asit sayımlarından (``CUB.amino_acid_counts``)
proteinin uzunluğunu, bileşimini ve kütlesini hesaplar. ``CUB``'da her gen
proteinle birebir örtüşen ORF'ye indirgendiği için bu sayımlar proteinin
kendisidir (stop hariç).

E-A. Sayım ve kütle:

====  ================================================  ==========================
 No   Ölçüm                                             Metot
====  ================================================  ==========================
 1    Uzunluk (N): Σ n_i                                ``length``
 2    Amino asit bileşimi (%): n_i / N × 100            ``composition``
 3    Amino asit grubu bileşimi: Σ_grup n_i / N         ``group_composition``
 4    Ortalama moleküler ağırlık: Σ n_i·M_i + H2O       ``molecular_weight``
 5    Monoizotopik kütle: aynı, monoizotopik kütlelerle ``monoisotopic_mass``
 6    Ortalama kalıntı kütlesi: MW / N                  ``mean_residue_mass``
 7    Cys, Trp, Tyr sayısı; disülfit eşi ⌊nCys / 2⌋     ``residue_counts``,
                                                        ``disulfide_pairs``
====  ================================================  ==========================

Kütleler Da birimindedir ve ``crn.io.crn_io_const`` içindeki kalıntı kütlesi
tablolarından gelir (``RESIDUE_AVERAGE_MASS``, ``RESIDUE_MONOISOTOPIC_MASS``);
gruplar ``AMINO_ACID_GROUPS``'tan. Tablolar parametreyle değiştirilebilir.

``gene`` verilmezse tüm genler birlikte hesaplanır: bileşim tüm proteinlerin
toplamından, kütle de her protein ayrı bir zincir olarak (her biri için bir
H2O) toplanarak bulunur.

Kullanım::

    from crn.cub import CUB, ProteinMass

    mass = ProteinMass(cub)
    mass.molecular_weight(gene="MEFV_NM_000243.3")
    mass.group_composition(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

from typing import (Dict, Mapping, Optional, Sequence)

from crn.io.crn_io_const import (AMINO_ACIDS, AMINO_ACID_GROUPS,
                                 RESIDUE_AVERAGE_MASS, RESIDUE_MONOISOTOPIC_MASS,
                                 WATER_AVERAGE_MASS, WATER_MONOISOTOPIC_MASS)
from crn.cub.crn_cub_counts import CUB


class ProteinMass:
    """Bir ``CUB`` nesnesindeki proteinlerin uzunluğu, bileşimi ve kütlesi.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        if not isinstance(cub, CUB):
            raise TypeError(f"cub bir CUB olmalı, {type(cub).__name__} verildi")
        self.__cub = cub

    @property
    def cub(self) -> CUB:
        return self.__cub

    def chain_count(self, gene: Optional[str] = None) -> int:
        """Hesaba giren protein zinciri sayısı: gene verilmezse gen sayısı, verilirse 1.

        Raises:
            KeyError: Gen yoksa.
        """
        if gene is None:
            return len(self.__cub.genes)
        self.__cub.codons(gene)         # gen yoksa KeyError
        return 1

    # --- E-A. Sayım ve kütle ------------------------------------------------------

    def length(self, gene: Optional[str] = None,
               amino_acids: Sequence[str] = AMINO_ACIDS) -> int:
        """1. Uzunluk (N): amino asit sayılarının toplamı (stop hariç)."""
        return sum(self.__cub.amino_acid_counts(gene, amino_acids).values())

    def composition(self, gene: Optional[str] = None,
                    amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, float]:
        """2. Amino asit bileşimi: her amino asidin yüzdesi, n_i / N × 100."""
        counts = self.__cub.amino_acid_counts(gene, amino_acids)
        total = sum(counts.values())
        return {aa: (n / total * 100 if total else 0.0) for aa, n in counts.items()}

    def group_composition(self, gene: Optional[str] = None,
                          groups: Mapping[str, str] = AMINO_ACID_GROUPS,
                          amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, float]:
        """3. Amino asit grubu bileşimi: gruptaki amino asitlerin toplam oranı (0-1).

        Args:
            groups: Grup adı -> gruptaki amino asitler (tek harfli kodlar).
        """
        counts = self.__cub.amino_acid_counts(gene, amino_acids)
        total = sum(counts.values())
        return {name: (sum(counts.get(aa, 0) for aa in members) / total if total else 0.0)
                for name, members in groups.items()}

    def scale_average(self, scale: Mapping[str, float],
                      gene: Optional[str] = None) -> float:
        """Kalıntı ölçeğinin bileşimle ağırlıklı ortalaması: Σ(n_i × ölçek_i) / N.

        GRAVY, hidropati, Boman, kalıntı hacmi, TOP-IDP gibi ölçümlerin ortak
        formülüdür. Yalnızca ``scale`` içindeki amino asitler sayılır.

        Args:
            scale: Amino asit -> değer.
        """
        counts = self.__cub.amino_acid_counts(gene, tuple(scale))
        n = sum(counts.values())
        return sum(count * scale[aa] for aa, count in counts.items()) / n if n else 0.0

    def molecular_weight(self, gene: Optional[str] = None,
                         masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS,
                         water: float = WATER_AVERAGE_MASS) -> float:
        """4. Ortalama moleküler ağırlık (Da): Σ(n_i × M_i) + H2O.

        Args:
            masses: Amino asit -> kalıntı kütlesi.
            water: Zincir uçlarındaki H2O'nun kütlesi; her protein için bir kez
                eklenir.

        Raises:
            KeyError: Sayılan bir amino asidin kütlesi ``masses`` içinde yoksa.
        """
        counts = self.__cub.amino_acid_counts(gene, tuple(masses))
        if not sum(counts.values()):
            return 0.0
        return sum(n * masses[aa] for aa, n in counts.items()) + water * self.chain_count(gene)

    def monoisotopic_mass(self, gene: Optional[str] = None,
                          masses: Mapping[str, float] = RESIDUE_MONOISOTOPIC_MASS,
                          water: float = WATER_MONOISOTOPIC_MASS) -> float:
        """5. Monoizotopik kütle (Da): ``molecular_weight`` formülü, monoizotopik kütlelerle."""
        return self.molecular_weight(gene, masses, water)

    def mean_residue_mass(self, gene: Optional[str] = None,
                          masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS,
                          water: float = WATER_AVERAGE_MASS) -> float:
        """6. Ortalama kalıntı kütlesi (Da): MW / N."""
        n = self.length(gene, tuple(masses))
        return self.molecular_weight(gene, masses, water) / n if n else 0.0

    def residue_counts(self, gene: Optional[str] = None,
                       residues: Sequence[str] = "CWY") -> Dict[str, int]:
        """7. Seçilen kalıntıların sayısı; varsayılan Cys, Trp, Tyr."""
        return self.__cub.amino_acid_counts(gene, residues)

    def disulfide_pairs(self, gene: Optional[str] = None) -> int:
        """7. Olası disülfit bağı sayısı: ⌊nCys / 2⌋ (tüm Cys'lerin eşleştiği varsayımı).

        gene verilmezse her protein için ayrı hesaplanıp toplanır (farklı
        zincirlerdeki Cys'ler eşleştirilmez).
        """
        if gene is None:
            return sum(self.disulfide_pairs(g.name) for g in self.__cub.genes)
        return self.residue_counts(gene, "C")["C"] // 2

    def __repr__(self) -> str:
        return f"ProteinMass({self.__cub!r})"
