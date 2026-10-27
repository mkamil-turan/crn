from __future__ import annotations

from typing import (Any, Dict, Self)

from crn.io.crn_io_seqs import Sequence
from crn.io.crn_io_dna import DNA
from crn.io.crn_io_prt import PRT
from crn.io.crn_io_type import Molecule, Strand, GenomicStrand
from crn.io.crn_io_const import (Maps, GeneticCode, RNA_ALPHABET, STOP_SYMBOL)


class RNA(Sequence):

    _molecule = Molecule.RNA
    _alphabet = RNA_ALPHABET

    def __init__(self,
                 seq, name = "",
                 strand:Strand = Strand.FivePrime,
                 genomic_strand:GenomicStrand = GenomicStrand.Unknown,
                 maps:Maps = Maps.DEFAULT_COMPLEMENT_MAP_FOR_RNA):
        super().__init__(seq, name)
        self.strand = strand
        self.genomic_strand = genomic_strand
        self.__maps = maps

    def _root_fields(self) -> Dict[str, Any]:
        return {"strand": self.__strand,
                "genomic_strand": self.__genomic_strand,
                "maps": self.__maps}

    @classmethod
    def from_sequence(cls, sequence:Sequence) -> Self:
        """Bir Sequence'ten aynı ad, dizi ve features kopyasıyla RNA oluşturur.

        features içinde ``strand`` anahtarında bir ``Strand`` varsa RNA'nın
        strand'i o olur; yoksa ``Strand.FivePrime``. ``genomic_strand``
        anahtarında bir ``GenomicStrand`` varsa o kullanılır; yoksa
        ``GenomicStrand.Unknown``.

        Raises:
            ValueError: Dizide RNA'da olmayan bir karakter varsa (ör. T).
        """
        strand = sequence.features.first("strand")
        if not isinstance(strand, Strand):
            strand = Strand.FivePrime
        genomic_strand = sequence.features.first("genomic_strand")
        if not isinstance(genomic_strand, GenomicStrand):
            genomic_strand = GenomicStrand.Unknown
        result = cls(seq=sequence.seq, name=sequence.name,
                     strand=strand, genomic_strand=genomic_strand)
        result.features = sequence.features.copy()
        return result

    @property
    def strand(self) -> Strand:
        return self.__strand

    @strand.setter
    def strand(self, value:Strand) -> None:
        if not isinstance(value, Strand):
            raise TypeError(f"strand bir Strand olmalı, {type(value).__name__} verildi")
        self.__strand = value

    @property
    def genomic_strand(self) -> GenomicStrand:
        """RNA'nın transkript edildiği genomik zincir (+, -, ., ?)."""
        return self.__genomic_strand

    @genomic_strand.setter
    def genomic_strand(self, value:GenomicStrand) -> None:
        if not isinstance(value, GenomicStrand):
            raise TypeError(f"genomic_strand bir GenomicStrand olmalı, {type(value).__name__} verildi")
        self.__genomic_strand = value

    @property
    def maps(self) -> Maps:
        return self.__maps

    def complement(self) -> Self:
        """Bu RNA'yı yerinde tümleyiciye çevirir ve self'i döndürür (zincirleme için).

        Yalnızca yazılış yönü (strand) çevrilir; genomic_strand RNA'nın
        geldiği genomik zinciri gösterdiği için değişmez.
        """
        self.seq = self.seq.translate(self.__maps.value)
        self.__strand = (Strand.FivePrime
                         if self.__strand == Strand.ThreePrime
                         else Strand.ThreePrime)
        return self

    def complementary(self) -> Self:
        """Tümleyici olarak yeni bir RNA döndürür; bu nesneyi değiştirmez.

        Yeni RNA aynı adı, genomic_strand'i, maps'i ve features'ın bağımsız
        bir kopyasını taşır.
        """
        return self._new(self.seq).complement()

    def reverse_transcription(self, maps: Maps = Maps.REVERSE_TRANSCRIPTION_MAP) -> DNA:
        """Ters transkripsiyon: bu RNA'dan tümleyici cDNA (ilk zincir) üretir.

        Ters transkriptaz RNA'yı kalıp olarak kullanır ve ona tümleyici,
        antiparalel bir DNA sentezler: A->T, U->A, G->C, C->G. Sonuç RNA'nın
        altına hizalı yazılır; bu yüzden strand ters çevrilir
        (5'->3' RNA -> 3'->5' cDNA).

        mRNA ile aynı anlamlı (U yerine T) ikinci zincir için sonucun
        ``complementary()`` metodu kullanılır.

        Args:
            maps: RNA bazlarını DNA bazlarına çeviren tablo.

        Yeni DNA aynı adı ve features'ın kopyasını taşır; yalnızca RNA'ya ait
        ``genomic_strand`` anahtarı features'tan çıkarılır.
        """
        strand = (Strand.FivePrime
                  if self.__strand == Strand.ThreePrime
                  else Strand.ThreePrime)
        result = DNA(self.seq.translate(maps.value),
                     name=self.name, strand=strand)
        result.features = self.features.copy()
        result.features.remove("genomic_strand")
        return result

    def translation(self, genetic_code: GeneticCode = GeneticCode.STANDARD) -> PRT:
        """Translasyon: bu RNA'yı seçilen genetik koda göre proteine çevirir.

        Dizi ilk nükleotitten başlanarak son tam kodona kadar üçerli okunur.
        Başlangıç ya da stop kodonu aranmaz; stop kodonları ``*``
        olarak yazılır ve okuma devam eder.

        Ribozom RNA'yı 5'->3' yönünde okur; RNA 3'->5' yazılmışsa önce ters
        çevrilir. Boşluk (``-``) karakterleri okumadan önce atılır. Belirsiz
        baz (N, R, Y ...) içeren kodonlar ``X`` olur. Sondaki 3'ten kısa kodon
        parçası okunmaz. Başlangıç kodonları özel işlenmez; her kodon
        tablodaki amino asidine çevrilir (ör. GUG her yerde V olur).

        Args:
            genetic_code: Kullanılacak genetik kod; varsayılan standart kod.

        Yeni PRT aynı adı ve features'ın kopyasını taşır; protein için anlamı
        olmayan ``strand`` ve ``genomic_strand`` anahtarları çıkarılır.
        """
        seq = self.seq.replace("-", "")
        if self.__strand == Strand.ThreePrime:
            seq = seq[::-1]
        codons = [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]
        amino_acids = [genetic_code.table.get(codon, "X") for codon in codons]

        result = PRT("".join(amino_acids), name=self.name)
        result.features = self.features.copy()
        result.features.remove("strand")
        result.features.remove("genomic_strand")
        return result

    def match_orf(self, protein: PRT | str,
                  genetic_code: GeneticCode = GeneticCode.STANDARD) -> RNA | None:
        """Çevirisi ``protein`` ile birebir aynı olan ORF'yi bulur.

        RNA kodlayan DNA'ya çevrilip ``DNA.find_orfs`` ile taranır; yalnızca
        RNA'nın kendi (+) zincirindeki ORF'lere bakılır (mRNA'nın tümleyicisi
        çevrilmez). Karşılaştırmada proteinin ve ORF çevirisinin sonundaki
        stop simgesi atılır.

        Args:
            protein: Protein (``PRT``) ya da tek harfli amino asit dizisi.
            genetic_code: ORF arama ve translasyon için genetik kod.

        Returns:
            ORF'yi 5'->3' bir RNA olarak (stop kodonu dahil) ya da örtüşen ORF
            yoksa None. Yeni RNA bu RNA'nın adını, genomic_strand'ini ve
            features kopyasını taşır; ayrıca ``orf_start``, ``orf_end`` (ORF'nin
            bu RNA üzerindeki 1 tabanlı konumu) yazılır, ``cds`` ``(1, uzunluk)``
            olur, ``utr5`` ve ``utr3`` silinir.
        """
        target = (protein.seq if isinstance(protein, PRT) else protein.upper()).rstrip(STOP_SYMBOL)
        seq = self.seq.replace("-", "")
        if self.__strand == Strand.ThreePrime:
            seq = seq[::-1]

        dna = DNA(seq.translate(Maps.RNA_TO_CODING_MAP.value), name=self.name)
        for orf in dna.find_orfs(genetic_code=genetic_code):
            if orf.features.first("orf_strand") != "+":
                continue
            mrna = orf.transcription(template=False, genomic_strand=self.__genomic_strand)
            if mrna.translation(genetic_code).seq.rstrip(STOP_SYMBOL) == target:
                mrna.features = self.features.copy()
                mrna.features.remove("utr5")
                mrna.features.remove("utr3")
                mrna.features.set("cds", (1, len(mrna)))
                mrna.features.set("orf_start", orf.features.first("orf_start"))
                mrna.features.set("orf_end", orf.features.first("orf_end"))
                return mrna
        return None

    def __repr__(self):
        return (f"{super().__repr__()} {self.strand.value[0]}->{self.strand.value[1]}"
                f" ({self.genomic_strand.value[0]})")
