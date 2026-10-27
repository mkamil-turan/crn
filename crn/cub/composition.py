"""Kodon kullanım yanlılığı (codon usage bias) — B. Nükleotit bileşimi.

Bir ``CUB`` nesnesinin kodon sayımlarından (``CUB.codon_counts``) nükleotit
bileşimini hesaplar. Tüm ölçümler L ile aynı kodonlar üzerindendir: stop
kodonları ve belirsiz bazlı kodonlar sayılmaz.

B. Nükleotit bileşimi:

====  ==============================================  ======================
 No   Ölçüm                                           Metot
====  ==============================================  ======================
 6    GC%: (G + C) / toplam nükleotit × 100           ``gc_percent``
 7    GC1: 1. konumdaki G + C / L                     ``gc1``
 8    GC2: 2. konumdaki G + C / L                     ``gc2``
 9    GC12: (GC1 + GC2) / 2                           ``gc12``
 10   GC3: 3. konumdaki G + C / L                     ``gc3``
 11   GC3s: eş anlamlı kodonlarda 3. konum G + C      ``gc3s``
 12   A3s, U3s, G3s, C3s: 3. konum baz / eş anl. kod  ``base3s``
 13   AT3s: 1 − GC3s                                  ``at3s``
====  ==============================================  ======================

GC% yüzde, diğerleri 0-1 arası orandır. Eş anlamlı kodonlar
``CUB.synonymous_codons`` ile belirlenir: tek kodonlu amino asitler (standart
kodda Met, Trp) ve stoplar dışarıda kalır. Diziler RNA olduğu için T yerine
U yazılır (T3s = ``base3s()["U"]``). ``base3s`` değerlerinin toplamı 1'dir ve
``G3s + C3s = GC3s``.

Kullanım::

    from crn.cub import CUB, Composition

    cub = CUB(isoforms)
    comp = Composition(cub)
    comp.gc3s(gene="MEFV_NM_000243.3")
    comp.summary()                          # tüm genler birlikte, 6-13 hepsi
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

from crn.constants import CODONS, RNA_BASES
from crn.cub.core import CUB


def _ratio(part: float, total: float) -> float:
    return part / total if total else 0.0


class Composition:
    """Bir ``CUB`` nesnesinin nükleotit bileşimi.

    Metotlar ``gene`` verilmezse tüm genleri birlikte, verilirse yalnızca o
    geni hesaplar.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        if not isinstance(cub, CUB):
            raise TypeError(f"cub bir CUB olmalı, {type(cub).__name__} verildi")
        self.cub = cub

    # --- Temel sayım ----------------------------------------------------------

    def nucleotide_counts(self, gene: Optional[str] = None,
                          position: Optional[int] = None,
                          codons: Sequence[str] = CODONS,
                          bases: str = RNA_BASES) -> Dict[str, int]:
        """``codons`` içindeki amino asit kodlayan kodonlarda baz sayıları.

        Args:
            position: None ise kodonun üç konumu birden; 1, 2 ya da 3 ise o konum.
            codons: Sayıma girecek kodonlar (stoplar her durumda çıkarılır).
            bases: Sayılacak bazlar, sırasıyla.

        Raises:
            ValueError: position 1, 2, 3 ya da None değilse.
        """
        if position not in (None, 1, 2, 3):
            raise ValueError(f"position 1, 2, 3 ya da None olmalı, {position!r} verildi")
        counts = dict.fromkeys(bases, 0)
        for codon, n in self.cub.codon_counts(gene, self.cub.sense_codons(codons)).items():
            for base in (codon if position is None else codon[position - 1]):
                if base in counts:
                    counts[base] += n
        return counts

    def nucleotide_freqs(self, gene: Optional[str] = None,
                         position: Optional[int] = None,
                         codons: Sequence[str] = CODONS,
                         bases: str = RNA_BASES) -> Dict[str, float]:
        """Baz oranları; toplamları 1'dir (bkz. ``nucleotide_counts``)."""
        counts = self.nucleotide_counts(gene, position, codons, bases)
        total = sum(counts.values())
        return {base: _ratio(n, total) for base, n in counts.items()}

    def _gc(self, gene: Optional[str], position: Optional[int],
            codons: Sequence[str] = CODONS) -> float:
        freqs = self.nucleotide_freqs(gene, position, codons)
        return freqs["G"] + freqs["C"]

    # --- B. Nükleotit bileşimi --------------------------------------------------

    def gc_percent(self, gene: Optional[str] = None) -> float:
        """6. GC%: (G + C) / toplam nükleotit × 100."""
        return self._gc(gene, None) * 100

    def gc1(self, gene: Optional[str] = None) -> float:
        """7. GC1: kodonların 1. konumundaki G + C / L."""
        return self._gc(gene, 1)

    def gc2(self, gene: Optional[str] = None) -> float:
        """8. GC2: kodonların 2. konumundaki G + C / L."""
        return self._gc(gene, 2)

    def gc12(self, gene: Optional[str] = None) -> float:
        """9. GC12: (GC1 + GC2) / 2."""
        return (self.gc1(gene) + self.gc2(gene)) / 2

    def gc3(self, gene: Optional[str] = None) -> float:
        """10. GC3: kodonların 3. konumundaki G + C / L."""
        return self._gc(gene, 3)

    def gc3s(self, gene: Optional[str] = None) -> float:
        """11. GC3s: eş anlamlı kodonlarda (Met, Trp ve stop hariç) 3. konum G + C."""
        return self._gc(gene, 3, self.cub.synonymous_codons())

    def base3s(self, gene: Optional[str] = None) -> Dict[str, float]:
        """12. A3s, U3s, G3s, C3s: eş anlamlı kodonlarda 3. konumdaki baz / eş anlamlı kodon sayısı."""
        return self.nucleotide_freqs(gene, 3, self.cub.synonymous_codons())

    def at3s(self, gene: Optional[str] = None) -> float:
        """13. AT3s: 1 − GC3s."""
        return 1 - self.gc3s(gene)

    def summary(self, gene: Optional[str] = None) -> Dict[str, float]:
        """6-13 arasındaki tüm ölçümler tek sözlükte."""
        base3s = self.base3s(gene)
        return {
            "GC%": self.gc_percent(gene),
            "GC1": self.gc1(gene),
            "GC2": self.gc2(gene),
            "GC12": self.gc12(gene),
            "GC3": self.gc3(gene),
            "GC3s": self.gc3s(gene),
            **{f"{base}3s": base3s[base] for base in "AUGC"},
            "AT3s": self.at3s(gene),
        }

    def __repr__(self) -> str:
        return f"Composition({self.cub!r})"
