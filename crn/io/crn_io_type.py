from enum import Enum

class Molecule(Enum):
    DNA = ("Dna","DNA")
    RNA = ("Rna", "RNA")
    PRT = ("Protein", "PRT")
    SEQ = ("Sequence", 'SEQ')

class Strand(Enum):
    """Dizinin yazılış yönü: hangi uçtan hangi uca okunduğu."""
    FivePrime = ("5'", "3'", "Five prime")
    ThreePrime = ("3'", "5'", "Three prime")

class GenomicStrand(Enum):
    """RNA'nın transkript edildiği genomik zincir (GFF/GTF/BED notasyonu).

    Yalnızca RNA özelliğidir; DNA dizileri için kullanılmaz.

    Yazılış yönünden (Strand) bağımsızdır: ``+`` referans genomda yazılan
    zinciri, ``-`` onun tümleyenini gösterir. Örneğin MEFV (16p13.3) ``-``
    zincirde yer alır; mRNA'sı ``-`` zincirin 5'→3' dizisiyle aynıdır.
    """
    Plus = ("+", "Plus")
    Minus = ("-", "Minus")
    Unstranded = (".", "Unstranded")   # zincir bilgisi anlamsız (ör. tekrar bölgesi)
    Unknown = ("?", "Unknown")         # zincir anlamlı ama bilinmiyor
