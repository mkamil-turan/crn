"""Dizi alfabeleri, kodonlar, tümleme tabloları ve genetik kodlar.

Protein ve tRNA ölçümlerinin sayısal tabloları (kütleler, pKa, ölçekler ...)
``crn.tables`` içindedir.
"""
from __future__ import annotations

from enum import Enum
from functools import cached_property
from typing import Dict, FrozenSet, List, Tuple

__all__ = ["DNA_ALPHABET", "RNA_ALPHABET", "PRT_ALPHABET", "DNA_BASES", "RNA_BASES",
           "CODONS", "AMINO_ACIDS", "STOP_SYMBOL", "GAP", "TRANSCRIPT_PRIORITY",
           "DNA_COMPLEMENT", "RNA_COMPLEMENT", "DNA_TO_RNA", "RNA_TO_DNA",
           "GeneticCode"]

# IUPAC nükleotit kodları ve boşluk (gap) karakteri.
DNA_ALPHABET = frozenset("ACGTRYKMSWBDHVN-")
# RNA'da T yerine U kullanılır; diğer IUPAC kodları aynıdır.
RNA_ALPHABET = frozenset("ACGURYKMSWBDHVN-")
# 20 standart amino asit, IUPAC belirsizlik kodları (B=Asx, Z=Glx, J=Xle,
# X=herhangi), U=selenosistein, O=pirolizin, '*' stop ve '-' boşluk.
PRT_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY" "BZJXUO" "*-")

GAP = "-"
# Kodon tablolarında ve proteinlerde stop kodonunun simgesi.
STOP_SYMBOL = "*"

# Belirsizliği olmayan bazlar, NCBI genetik kod tablolarının sırasıyla.
DNA_BASES = "TCAG"
RNA_BASES = "UCAG"
# 64 kodon (RNA), RNA_BASES sırasıyla: UUU, UUC, UUA, UUG, UCU, ... GGG.
CODONS = tuple(a + b + c for a in RNA_BASES for b in RNA_BASES for c in RNA_BASES)
# 20 standart amino asit, tek harfli kodlarıyla, alfabetik.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# IUPAC tümleyicileri: R<->Y, K<->M, B<->V, D<->H; S, W, N ve '-' kendisidir.
DNA_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVN-", "TGCAYRMKSWVHDBN-")
RNA_COMPLEMENT = str.maketrans("ACGURYKMSWBDHVN-", "UGCAYRMKSWVHDBN-")
# Kodlayan zincir <-> mRNA: dizi aynı kalır, yalnızca T <-> U.
DNA_TO_RNA = str.maketrans("T", "U")
RNA_TO_DNA = str.maketrans("U", "T")

# Aynı CDS'i veren izoformlardan hangisinin tutulacağı: RefSeq önekleri,
# öncelik sırasıyla (NM_ deneysel kanıtlı, XM_ öngörülen mRNA). Listede
# olmayan önekler en sona kalır.
TRANSCRIPT_PRIORITY = ("NM_", "XM_")


def _codon_table(amino_acids: str) -> Dict[str, str]:
    """NCBI biçimindeki 64 harflik amino asit dizisinden kodon tablosu (CODONS sırasıyla)."""
    return dict(zip(CODONS, amino_acids))


class GeneticCode(Enum):
    """Translasyon, ORF araması ve kodon sayımında kullanılan genetik kodlar.

    Değer: (NCBI tablo no, ad, 64 harflik amino asit dizisi, başlangıç kodonları).
    Kodonlar RNA alfabesiyle (U) yazılır; stop kodonları ``*`` ile gösterilir.
    Türetilen kümeler (``stops``, ``families`` ...) ilk kullanımda bir kez
    hesaplanır.
    """

    # NCBI tablo 1. AUG dışındaki başlangıç kodonları: UUG, CUG, GUG.
    STANDARD = (1, "Standard",
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
                ("AUG", "UUG", "CUG", "GUG"))
    # NCBI tablo 2. Farklar: UGA=Trp, AUA=Met, AGA ve AGG=stop.
    VERTEBRATE_MITOCHONDRIAL = (2, "Vertebrate Mitochondrial",
                                "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
                                ("AUG", "AUU", "AUC", "AUA", "GUG"))

    def __repr__(self) -> str:
        return f"<GeneticCode.{self.name}: NCBI {self.value[0]}, {self.value[1]}>"

    @property
    def ncbi_id(self) -> int:
        return self.value[0]

    @cached_property
    def table(self) -> Dict[str, str]:
        """Kodon -> amino asit sözlüğü (RNA kodonları)."""
        return _codon_table(self.value[2])

    @cached_property
    def starts(self) -> FrozenSet[str]:
        """AUG dahil tüm başlangıç kodonları (RNA kodonları)."""
        return frozenset(self.value[3])

    @cached_property
    def stops(self) -> FrozenSet[str]:
        """Stop kodonları (RNA kodonları)."""
        return frozenset(c for c, aa in self.table.items() if aa == STOP_SYMBOL)

    @cached_property
    def sense_codons(self) -> Tuple[str, ...]:
        """Amino asit kodlayan kodonlar (stop hariç), CODONS sırasıyla."""
        return tuple(c for c in CODONS if c not in self.stops)

    @cached_property
    def families(self) -> Dict[str, Tuple[str, ...]]:
        """Amino asit -> onu kodlayan kodonlar (eş anlamlı kodon ailesi), CODONS sırasıyla."""
        families: Dict[str, List[str]] = {}
        for codon in self.sense_codons:
            families.setdefault(self.table[codon], []).append(codon)
        return {aa: tuple(codons) for aa, codons in families.items()}

    @cached_property
    def synonymous_codons(self) -> FrozenSet[str]:
        """Amino asidi birden çok kodonla kodlanan kodonlar (standart kodda Met, Trp ve stop hariç)."""
        return frozenset(c for c in self.sense_codons if len(self.families[self.table[c]]) > 1)

    @cached_property
    def fourfold_codons(self) -> FrozenSet[str]:
        """3. konumu tümüyle eş anlamlı xyN kutularının kodonları (standart kodda 32 kodon)."""
        return frozenset(c for c in self.sense_codons
                         if len({self.table[c[:2] + b] for b in RNA_BASES}) == 1)
