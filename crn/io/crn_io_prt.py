from __future__ import annotations

from typing import Self

from crn.io.crn_io_seqs import Sequence
from crn.io.crn_io_type import Molecule
from crn.io.crn_io_const import PRT_ALPHABET


class PRT(Sequence):
    """Protein dizisi; tek harfli amino asit kodlarıyla N-ucundan C-ucuna yazılır.

    Nükleik asitlerden farklı olarak strand ve tümleyici kavramı yoktur;
    bu yüzden root alanı yoktur.
    """

    _molecule = Molecule.PRT
    _alphabet = PRT_ALPHABET

    def __init__(self, seq, name = ""):
        super().__init__(seq, name)

    @classmethod
    def from_sequence(cls, sequence:Sequence) -> Self:
        """Bir Sequence'ten aynı ad, dizi ve features kopyasıyla protein oluşturur.

        Raises:
            ValueError: Dizide proteinde olmayan bir karakter varsa.
        """
        result = cls(seq=sequence.seq, name=sequence.name)
        result.features = sequence.features.copy()
        return result

    def __repr__(self):
        return f"{super().__repr__()} N->C"
