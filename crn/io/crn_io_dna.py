from __future__ import annotations

from typing import (Any, Dict, Self, TYPE_CHECKING)

from crn.io.crn_io_seqs import Sequence
from crn.io.crn_io_type import Molecule, Strand, GenomicStrand
from crn.io.crn_io_const import (Maps, GeneticCode, DNA_ALPHABET)

if TYPE_CHECKING:
    from crn.io.crn_io_rna import RNA


class DNA(Sequence):

    _molecule = Molecule.DNA
    _alphabet = DNA_ALPHABET

    def __init__(self,
                 seq, name = "",
                 strand:Strand = Strand.FivePrime,
                 maps:Maps = Maps.DEFAULT_COMPLEMENT_MAP_FOR_DNA):
        super().__init__(seq, name)
        self.strand = strand
        self.__maps = maps

    def _root_fields(self) -> Dict[str, Any]:
        return {"strand": self.__strand, "maps": self.__maps}

    @classmethod
    def from_sequence(cls, sequence:Sequence) -> Self:
        """Bir Sequence'ten aynı ad, dizi ve features kopyasıyla DNA oluşturur.

        features içinde ``strand`` anahtarında bir ``Strand`` varsa DNA'nın
        strand'i o olur; yoksa ``Strand.FivePrime``.

        Raises:
            ValueError: Dizide DNA'da olmayan bir karakter varsa.
        """
        strand = sequence.features.first("strand")
        if not isinstance(strand, Strand):
            strand = Strand.FivePrime
        result = cls(seq=sequence.seq, name=sequence.name, strand=strand)
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
    def maps(self) -> Maps:
        return self.__maps

    def complement(self) -> Self:
        """Bu DNA'yı yerinde tümleyiciye çevirir ve self'i döndürür (zincirleme için)."""
        self.seq = self.seq.translate(self.__maps.value)
        self.__strand = (Strand.FivePrime
                         if self.__strand == Strand.ThreePrime
                         else Strand.ThreePrime)
        return self

    def complementary(self) -> Self:
        """Tümleyici olarak yeni bir DNA döndürür; bu nesneyi değiştirmez.

        Yeni DNA aynı adı, maps'i ve features'ın bağımsız bir kopyasını taşır.
        """
        return self._new(self.seq).complement()

    def transcription(self, template: bool = True,
                      genomic_strand: GenomicStrand = GenomicStrand.Unknown,
                      maps: Maps | None = None) -> RNA:
        """Transkripsiyon: bu DNA'dan RNA üretir.

        Args:
            template: True ise bu DNA kalıp (template, antisense) zincirdir.
                RNA polimeraz gibi ona tümleyici, antiparalel RNA üretilir
                (A->U, T->A, G->C, C->G) ve strand ters çevrilir. Böylece
                ``rna.reverse_transcription().transcription()`` ilk RNA'yı
                geri verir.
                False ise bu DNA kodlayan (coding, sense) zincirdir; NCBI gibi
                kaynaklardan gelen gen dizileri çoğunlukla böyledir. RNA,
                kodlayan zincirle aynı diziye sahip olduğu için yalnızca T->U
                yapılır ve strand korunur.
            genomic_strand: Oluşan RNA'nın transkript edildiği genomik zincir
                (yalnızca RNA özelliğidir).
            maps: DNA bazlarını RNA bazlarına çeviren tablo. None ise
                ``template`` değerine göre ``Maps.TRANSCRIPTION_MAP`` ya da
                ``Maps.CODING_TO_RNA_MAP`` kullanılır.

        Yeni RNA aynı adı ve features'ın kopyasını taşır.
        """
        from crn.io.crn_io_rna import RNA

        if maps is None:
            maps = Maps.TRANSCRIPTION_MAP if template else Maps.CODING_TO_RNA_MAP
        seq = self.seq.translate(maps.value)
        if template:
            strand = (Strand.FivePrime
                      if self.__strand == Strand.ThreePrime
                      else Strand.ThreePrime)
        else:
            strand = self.__strand

        result = RNA(seq, name=self.name, strand=strand, genomic_strand=genomic_strand)
        result.features = self.features.copy()
        return result

    def find_orfs(self, minimum_length_of_orf: int = -1,
                  genetic_code: GeneticCode = GeneticCode.STANDARD,
                  alternative_starts: bool = False) -> list[Self]:
        """Bu DNA'nın iki zincirindeki tüm geçerli açık okuma çerçevelerini (ORF) bulur.

        Geçerli ORF: bir başlangıç kodonuyla başlayan, aynı okuma çerçevesindeki
        ilk stop kodonuyla biten dizidir. Bu zincir ve tümleyici zincir 5'->3'
        yönünde, üçer çerçeveden toplam altı çerçevede taranır. Stop kodonu
        olmadan dizinin sonuna ulaşan kısmi ORF'ler alınmaz. Aynı stop kodonunu
        paylaşan iç içe ORF'lerin hepsi ayrı ORF olarak döner. Dairesel diziler
        desteklenmez. Boşluk (``-``) karakterleri taramadan önce atlanır.

        Args:
            minimum_length_of_orf: Stop kodonu dahil en küçük ORF uzunluğu
                (nükleotit). -1 ise sınır uygulanmaz.
            genetic_code: Stop ve başlangıç kodonlarının alınacağı genetik kod.
            alternative_starts: False ise yalnızca ATG başlangıç kabul edilir;
                True ise genetik kodun tüm başlangıç kodonları (standart kodda
                TTG, CTG, GTG) da kabul edilir.

        Returns:
            Konuma göre sıralı ORF listesi. Her ORF kodlayan yönde, 5'->3'
            yazılmış bir DNA'dır; adı ``<ad>_ORF<n>`` olur, maps'i ve
            features'ın kopyasını taşır. Ayrıca features'a şunlar yazılır:

            - ``orf_start``, ``orf_end``: ORF'un bu DNA'nın ``seq`` dizisinde
              kapladığı bölgenin ilk ve son konumu (1 tabanlı, ilk nükleotitten
              uzaklık). Tümleyici zincirdeki ORF'lerde de ``orf_start <= orf_end``.
            - ``orf_strand``: ``"+"`` bu zincirde, ``"-"`` tümleyici zincirde
            - ``orf_frame``: o zincirin 5' ucuna göre okuma çerçevesi (0, 1, 2)

        Raises:
            ValueError: minimum_length_of_orf -1'den küçükse.
        """
        if minimum_length_of_orf < -1:
            raise ValueError(f"minimum_length_of_orf -1 ya da 0'dan büyük olmalı, "
                             f"{minimum_length_of_orf} verildi")

        # Bu zincirin 5'->3' okunuşu ve her bazın self.seq içindeki indeksi.
        seq = self.seq
        order = range(len(seq)) if self.__strand == Strand.FivePrime else range(len(seq) - 1, -1, -1)
        index = [i for i in order if seq[i] != "-"]
        forward = "".join(seq[i] for i in index)
        # Tümleyici zincirin 5'->3' okunuşu: ters tümleyici.
        reverse = forward[::-1].translate(self.__maps.value)

        starts = genetic_code.starts if alternative_starts else {"AUG"}
        starts = {c.replace("U", "T") for c in starts}
        stops = {c.replace("U", "T") for c in genetic_code.stops}

        found = []
        for sign, strand_seq, strand_index in (("+", forward, index), ("-", reverse, index[::-1])):
            for frame in range(3):
                pending = []
                for i in range(frame, len(strand_seq) - 2, 3):
                    codon = strand_seq[i:i + 3]
                    if codon in stops:
                        for start in pending:
                            first, last = strand_index[start], strand_index[i + 2]
                            found.append((min(first, last) + 1, max(first, last) + 1,
                                          sign, frame, strand_seq[start:i + 3]))
                        pending = []
                    elif codon in starts:
                        pending.append(i)
        found.sort()

        if minimum_length_of_orf != -1:
            found = [orf for orf in found if len(orf[4]) >= minimum_length_of_orf]

        orfs = []
        prefix = f"{self.name}_" if self.name else ""
        for n, (start, end, sign, frame, orf_seq) in enumerate(found, start=1):
            orf = self.__class__(orf_seq, name=f"{prefix}ORF{n}",
                                 strand=Strand.FivePrime, maps=self.__maps)
            orf.features = self.features.copy()
            orf.features.set("orf_start", start)
            orf.features.set("orf_end", end)
            orf.features.set("orf_strand", sign)
            orf.features.set("orf_frame", frame)
            orfs.append(orf)
        return orfs

    def __repr__(self):
        return f"{super().__repr__()} {self.strand.value[0]}->{self.strand.value[1]}"
