"""Dizi sınıfları: Sequence, DNA, RNA ve PRT.

::

    Sequence                 ad, dizi, features; alfabe denetimi; FASTA/CRN
     ├── NucleicAcid         strand (yazılış yönü), complement, reverse_complement
     │    ├── DNA            transcription, find_orfs
     │    └── RNA            genomic_strand, translation, reverse_transcription, match_orf
     └── PRT                 protein (N->C)

Sınıfın kendisinde yalnızca dizinin biyolojik anlamını değiştiren "root"
alanları vardır (``strand``, RNA'da ``genomic_strand``); diğer her bilgi
``features`` içindedir. Yeni nesne döndüren işlemler (dilimleme,
``complementary``, ``transcription`` ...) ``features``'ın bağımsız bir
kopyasını taşır.
"""
from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterator, List, Optional, Self

from crn.constants import (DNA_ALPHABET, RNA_ALPHABET, PRT_ALPHABET, GAP, STOP_SYMBOL,
                           DNA_COMPLEMENT, RNA_COMPLEMENT, DNA_TO_RNA, RNA_TO_DNA,
                           GeneticCode)
from crn.features import Features

__all__ = ["Molecule", "Strand", "GenomicStrand", "Sequence", "NucleicAcid",
           "DNA", "RNA", "PRT"]


class Molecule(Enum):
    DNA = "DNA"
    RNA = "RNA"
    PRT = "PRT"
    SEQ = "SEQ"


class Strand(Enum):
    """Dizinin yazılış yönü: hangi uçtan hangi uca okunduğu."""
    FivePrime = ("5'", "3'", "Five prime")
    ThreePrime = ("3'", "5'", "Three prime")

    @property
    def flipped(self) -> Strand:
        return Strand.ThreePrime if self is Strand.FivePrime else Strand.FivePrime


class GenomicStrand(Enum):
    """RNA'nın transkript edildiği genomik zincir (GFF/GTF/BED notasyonu).

    Yalnızca RNA özelliğidir ve yazılış yönünden (Strand) bağımsızdır: ``+``
    referans genomda yazılan zinciri, ``-`` onun tümleyenini gösterir. Örneğin
    MEFV (16p13.3) ``-`` zincirdedir; mRNA'sı ``-`` zincirin 5'->3' dizisidir.
    """
    Plus = ("+", "Plus")
    Minus = ("-", "Minus")
    Unstranded = (".", "Unstranded")   # zincir bilgisi anlamsız (ör. tekrar bölgesi)
    Unknown = ("?", "Unknown")         # zincir anlamlı ama bilinmiyor


# Molecule -> o molekülü temsil eden sınıf; alt sınıflar tanımlandıkça dolar.
SEQUENCE_CLASSES: Dict[Molecule, type] = {}


class Sequence:
    """Adı, dizisi ve ``features``'ı olan genel dizi; alfabe denetimi yapmaz."""

    molecule: Molecule = Molecule.SEQ
    # İzin verilen harfler; None ise denetim yapılmaz.
    alphabet: Optional[FrozenSet[str]] = None
    # Root alanları: __init__ parametresi -> tipi. CRN dosyasına yazılır,
    # kopyalarda korunur, from_sequence features'tan okur.
    root_types: Dict[str, type] = {}

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if "molecule" in cls.__dict__:
            SEQUENCE_CLASSES[cls.molecule] = cls

    def __init__(self, seq: str, name: str = ""):
        self._seq = self._check(seq)
        self.name = name
        self._features = Features()

    def _check(self, seq: str) -> str:
        """Diziyi büyük harfe çevirir; ``alphabet`` tanımlıysa harfleri doğrular."""
        seq = seq.upper()
        if self.alphabet is not None:
            invalid = sorted(set(seq) - self.alphabet)
            if invalid:
                raise ValueError(f"{self.molecule.value} dizisinde geçersiz "
                                 f"karakter(ler): {''.join(invalid)}")
        return seq

    # --- Alanlar --------------------------------------------------------------

    @property
    def seq(self) -> str:
        return self._seq

    @seq.setter
    def seq(self, value: str) -> None:
        self._seq = self._check(value)

    @property
    def features(self) -> Features:
        return self._features

    @features.setter
    def features(self, value: Features) -> None:
        if not isinstance(value, Features):
            raise TypeError(f"features bir Features olmalı, {type(value).__name__} verildi")
        self._features = value

    def root_fields(self) -> Dict[str, Any]:
        """Root alanları ve değerleri (ör. ``{"strand": Strand.FivePrime}``)."""
        return {key: getattr(self, key) for key in self.root_types}

    def _new(self, seq: str) -> Self:
        """Aynı sınıf, ad, root alanları ve features kopyasıyla yeni bir nesne."""
        result = type(self)(seq, name=self.name, **self.root_fields())
        result.features = self._features.copy()
        return result

    @classmethod
    def from_sequence(cls, sequence: Sequence) -> Self:
        """Başka bir diziden aynı ad, dizi ve features kopyasıyla ``cls`` oluşturur.

        Root alanları (ör. ``strand``) features'ta doğru tipte bir değer varsa
        oradan alınır, yoksa varsayılan kullanılır.

        Raises:
            ValueError: Dizide ``cls`` alfabesinde olmayan bir karakter varsa.
        """
        roots = {key: value for key, kind in cls.root_types.items()
                 if isinstance(value := sequence.features.first(key), kind)}
        result = cls(sequence.seq, name=sequence.name, **roots)
        result.features = sequence.features.copy()
        return result

    def sliding(self, window_size: int, step_size: int = 1) -> List[Self]:
        """Diziyi kayan pencerelere böler; her pencere ``s[a:b]`` ile aynı sınıftadır.

        Pencereler ``seq``'in yazıldığı sırayla, soldan sağa alınır. Yalnızca tam
        pencereler döner: sonda ``window_size``'dan kısa kalan parça alınmaz.

        Args:
            window_size: Pencere uzunluğu.
            step_size: İki pencerenin başlangıçları arasındaki uzaklık.

        Returns:
            Bu diziyle aynı sınıftan (DNA, RNA, PRT, Sequence), aynı ad ve root
            alanlarını, features'ın bağımsız kopyasını taşıyan diziler. Dizi
            ``window_size``'dan kısaysa liste boştur.

        Raises:
            ValueError: window_size ya da step_size 1'den küçükse.
        """
        if window_size < 1:
            raise ValueError(f"window_size 1 ya da daha büyük olmalı, {window_size} verildi")
        if step_size < 1:
            raise ValueError(f"step_size 1 ya da daha büyük olmalı, {step_size} verildi")
        return [self[start:start + window_size]
                for start in range(0, len(self) - window_size + 1, step_size)]

    # --- Dosyalar (biçimler crn.fileio içinde) ---------------------------------

    @classmethod
    def from_fasta_file(cls, file_name: str | Path) -> List[Self]:
        """FASTA dosyasındaki dizileri ``cls`` tipinde liste olarak okur (bkz. ``crn.fileio.read_fasta``)."""
        from crn.fileio import read_fasta
        return read_fasta(file_name, cls)

    def to_fasta(self, file_name: str | Path) -> None:
        """Diziyi FASTA dosyasına yazar (bkz. ``crn.fileio.write_fasta``)."""
        from crn.fileio import write_fasta
        write_fasta([self], file_name)

    def to_crn(self, file_name: str | Path) -> None:
        """Dizinin tüm içeriğini (tip, ad, dizi, root alanları, features) CRN dosyasına yazar."""
        from crn.fileio import write_crn
        write_crn(self, file_name)

    @classmethod
    def from_crn(cls, file_name: str | Path) -> Self:
        """``to_crn`` ile yazılmış dosyayı okur.

        ``Sequence.from_crn`` dosyadaki tipe uygun sınıfı döndürür; ``DNA.from_crn``
        gibi bir alt sınıftan çağrılırsa dosyadaki tip o sınıf olmalıdır.
        """
        from crn.fileio import read_crn
        return read_crn(file_name, cls)

    # --- İşleçler -------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._seq)

    def __eq__(self, value: object) -> bool:
        """Aynı tip ve aynı dizi ise eşittir; ``str`` ile büyük/küçük harf
        farkı gözetmeden karşılaştırılır. Ad ve features karşılaştırılmaz."""
        if isinstance(value, Sequence):
            return self.molecule is value.molecule and self._seq == value._seq
        if isinstance(value, str):
            return self._seq == value.upper()
        return NotImplemented

    # seq değiştirilebildiği için diziler hash'lenemez.
    __hash__ = None

    def __contains__(self, item: Sequence | str) -> bool:
        """``"ATG" in s`` ya da ``other in s`` biçiminde alt dizi arar."""
        if isinstance(item, Sequence):
            return item._seq in self._seq
        if isinstance(item, str):
            return item.upper() in self._seq
        raise TypeError(f"Sequence içinde aranamayan tip: {type(item).__name__}")

    def __iter__(self) -> Iterator[str]:
        return iter(self._seq)

    def __getitem__(self, key: int | slice) -> str | Self:
        """``s[i]`` tek harf, ``s[a:b]`` aynı sınıftan yeni dizi (features kopyasıyla)."""
        if isinstance(key, slice):
            return self._new(self._seq[key])
        return self._seq[key]

    def __repr__(self) -> str:
        seq = self._seq if len(self) < 80 else f"{self._seq[:40]}...{self._seq[-40:]}"
        return f"{self.name}: {seq} <{len(self)}> [{self.molecule.value}]"


class NucleicAcid(Sequence):
    """DNA ve RNA'nın ortak tabanı: yazılış yönü (strand) ve tümleme."""

    root_types = {"strand": Strand}
    # Tümleme tablosu (str.translate); alt sınıf belirler.
    complement_table: Dict[int, int] = {}

    def __init__(self, seq: str, name: str = "", strand: Strand = Strand.FivePrime):
        super().__init__(seq, name)
        self.strand = strand

    @property
    def strand(self) -> Strand:
        return self._strand

    @strand.setter
    def strand(self, value: Strand) -> None:
        if not isinstance(value, Strand):
            raise TypeError(f"strand bir Strand olmalı, {type(value).__name__} verildi")
        self._strand = value

    def seq_5to3(self) -> str:
        """Dizinin 5'->3' okunuşu, boşluklar (``-``) atılmış."""
        seq = self._seq.replace(GAP, "")
        return seq[::-1] if self._strand is Strand.ThreePrime else seq

    def complement(self) -> Self:
        """Yerinde tümleyiciye çevirir (strand da döner) ve self'i döndürür.

        Tümleyici, bu dizinin altına hizalı yazılır: 5'->3' ``ATGC`` -> 3'->5' ``TACG``.
        """
        self._seq = self._seq.translate(self.complement_table)
        self._strand = self._strand.flipped
        return self

    def complementary(self) -> Self:
        """Tümleyiciyi yeni bir nesne olarak döndürür (bkz. ``complement``)."""
        return self._new(self._seq).complement()

    def reverse_complement(self) -> Self:
        """Karşı zinciri bu diziyle aynı yönde okunacak biçimde döndürür.

        5'->3' ``ATGC`` -> 5'->3' ``GCAT``; yeni nesne, strand aynı kalır.
        """
        return self._new(self._seq.translate(self.complement_table)[::-1])

    def __repr__(self) -> str:
        return f"{super().__repr__()} {self._strand.value[0]}->{self._strand.value[1]}"


class DNA(NucleicAcid):

    molecule = Molecule.DNA
    alphabet = DNA_ALPHABET
    complement_table = DNA_COMPLEMENT

    def transcription(self, template: bool = True,
                      genomic_strand: GenomicStrand = GenomicStrand.Unknown) -> RNA:
        """Transkripsiyon: bu DNA'dan RNA üretir.

        Args:
            template: True ise bu DNA kalıp (antisense) zincirdir; ona tümleyici
                RNA üretilir ve strand döner, böylece
                ``rna.reverse_transcription().transcription()`` ilk RNA'yı verir.
                False ise kodlayan (sense) zincirdir (NCBI/GENCODE gen dizileri
                böyledir); yalnızca T->U yapılır, strand korunur.
            genomic_strand: Oluşan RNA'nın genomik zinciri.
        """
        if template:
            seq, strand = self._seq.translate(DNA_COMPLEMENT).translate(DNA_TO_RNA), self.strand.flipped
        else:
            seq, strand = self._seq.translate(DNA_TO_RNA), self.strand
        result = RNA(seq, name=self.name, strand=strand, genomic_strand=genomic_strand)
        result.features = self.features.copy()
        return result

    def find_orfs(self, minimum_length_of_orf: int = -1,
                  genetic_code: GeneticCode = GeneticCode.STANDARD,
                  alternative_starts: bool = False) -> List[Self]:
        """İki zincirdeki tüm açık okuma çerçevelerini (ORF) bulur.

        ORF: bir başlangıç kodonundan aynı çerçevedeki ilk stop kodonuna (dahil)
        kadar olan dizi. Bu zincir ve tümleyicisi 5'->3' yönünde, toplam altı
        çerçevede taranır. Stop'suz kısmi ORF'ler alınmaz; aynı stop'u paylaşan
        iç içe ORF'lerin hepsi döner. Boşluklar taramadan önce atlanır.

        Args:
            minimum_length_of_orf: Stop dahil en küçük uzunluk (nt); -1 sınırsız.
            genetic_code: Başlangıç ve stop kodonlarının alınacağı kod.
            alternative_starts: True ise kodun tüm başlangıç kodonları, False
                ise yalnızca ATG.

        Returns:
            Konuma göre sıralı, kodlayan yönde 5'->3' DNA listesi. Adı
            ``<ad>_ORF<n>``; features'a ``orf_start``/``orf_end`` (bu dizide
            1 tabanlı, start <= end), ``orf_strand`` (``+``/``-``) ve
            ``orf_frame`` (0, 1, 2) yazılır.

        Raises:
            ValueError: minimum_length_of_orf -1'den küçükse.
        """
        if minimum_length_of_orf < -1:
            raise ValueError(f"minimum_length_of_orf -1 ya da 0'dan büyük olmalı, "
                             f"{minimum_length_of_orf} verildi")

        # Bu zincirin 5'->3' okunuşu ve her bazın self.seq içindeki indeksi.
        seq = self._seq
        order = range(len(seq)) if self.strand is Strand.FivePrime else range(len(seq) - 1, -1, -1)
        index = [i for i in order if seq[i] != GAP]
        forward = "".join(seq[i] for i in index)
        reverse = forward[::-1].translate(DNA_COMPLEMENT)

        starts = genetic_code.starts if alternative_starts else {"AUG"}
        starts = {c.translate(RNA_TO_DNA) for c in starts}
        stops = {c.translate(RNA_TO_DNA) for c in genetic_code.stops}

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
            orf = type(self)(orf_seq, name=f"{prefix}ORF{n}")
            orf.features = self.features.copy()
            orf.features.set("orf_start", start)
            orf.features.set("orf_end", end)
            orf.features.set("orf_strand", sign)
            orf.features.set("orf_frame", frame)
            orfs.append(orf)
        return orfs


class RNA(NucleicAcid):

    molecule = Molecule.RNA
    alphabet = RNA_ALPHABET
    complement_table = RNA_COMPLEMENT
    root_types = {"strand": Strand, "genomic_strand": GenomicStrand}

    def __init__(self, seq: str, name: str = "", strand: Strand = Strand.FivePrime,
                 genomic_strand: GenomicStrand = GenomicStrand.Unknown):
        super().__init__(seq, name, strand)
        self.genomic_strand = genomic_strand

    @property
    def genomic_strand(self) -> GenomicStrand:
        """RNA'nın transkript edildiği genomik zincir (+, -, ., ?); complement ile değişmez."""
        return self._genomic_strand

    @genomic_strand.setter
    def genomic_strand(self, value: GenomicStrand) -> None:
        if not isinstance(value, GenomicStrand):
            raise TypeError(f"genomic_strand bir GenomicStrand olmalı, {type(value).__name__} verildi")
        self._genomic_strand = value

    def reverse_transcription(self) -> DNA:
        """Ters transkripsiyon: tümleyici, antiparalel cDNA (ilk zincir).

        A->T, U->A, G->C, C->G; strand döner (5'->3' RNA -> 3'->5' cDNA).
        mRNA ile aynı anlamlı ikinci zincir için sonucun ``complementary()``'si
        kullanılır. ``genomic_strand`` anahtarı features'tan çıkarılır.
        """
        result = DNA(self._seq.translate(RNA_COMPLEMENT).translate(RNA_TO_DNA),
                     name=self.name, strand=self.strand.flipped)
        result.features = self.features.copy()
        result.features.remove("genomic_strand")
        return result

    def translation(self, genetic_code: GeneticCode = GeneticCode.STANDARD) -> PRT:
        """Translasyon: ilk nükleotitten son tam kodona kadar üçerli okur.

        Başlangıç ya da stop aranmaz; stoplar ``*`` yazılır ve okuma sürer.
        3'->5' yazılmış RNA önce çevrilir, boşluklar atılır; belirsiz bazlı
        kodonlar ``X`` olur. ``strand`` ve ``genomic_strand`` anahtarları
        features'tan çıkarılır.
        """
        seq = self.seq_5to3()
        table = genetic_code.table
        result = PRT("".join(table.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3)),
                     name=self.name)
        result.features = self.features.copy()
        result.features.remove("strand")
        result.features.remove("genomic_strand")
        return result

    def match_orf(self, protein: PRT | str,
                  genetic_code: GeneticCode = GeneticCode.STANDARD) -> Optional[RNA]:
        """Çevirisi ``protein`` ile birebir aynı olan ORF'yi (AUG ... stop) bulur.

        Yalnızca RNA'nın kendi zincirine bakılır; sondaki stop simgesi iki
        tarafta da yok sayılır. Birden çok ORF uyarsa en baştaki seçilir.

        Returns:
            ORF'yi (stop dahil) 5'->3' bir RNA olarak ya da uyan ORF yoksa None.
            Yeni RNA bu RNA'nın adını, genomic_strand'ini ve features kopyasını
            taşır; ``orf_start``, ``orf_end`` (bu RNA üzerinde 1 tabanlı) ve
            ``cds`` = ``(1, uzunluk)`` yazılır, ``utr5`` ve ``utr3`` silinir.
        """
        target = (protein.seq if isinstance(protein, PRT) else protein.upper()).rstrip(STOP_SYMBOL)
        if not target or STOP_SYMBOL in target:
            return None
        seq = self.seq_5to3()
        table = genetic_code.table
        wanted = target + STOP_SYMBOL

        best = None
        for frame in range(3):
            translated = "".join(table.get(seq[i:i + 3], "X") for i in range(frame, len(seq) - 2, 3))
            at = translated.find(wanted)
            while at != -1:
                start = frame + 3 * at
                if seq[start:start + 3] == "AUG":
                    if best is None or start < best:
                        best = start
                    break
                at = translated.find(wanted, at + 1)
        if best is None:
            return None

        end = best + 3 * len(wanted)
        mrna = RNA(seq[best:end], name=self.name, genomic_strand=self.genomic_strand)
        mrna.features = self.features.copy()
        mrna.features.remove("utr5")
        mrna.features.remove("utr3")
        mrna.features.set("cds", (1, len(mrna)))
        mrna.features.set("orf_start", best + 1)
        mrna.features.set("orf_end", end)
        return mrna

    def __repr__(self) -> str:
        return f"{super().__repr__()} ({self.genomic_strand.value[0]})"


class PRT(Sequence):
    """Protein dizisi; tek harfli kodlarla N-ucundan C-ucuna yazılır."""

    molecule = Molecule.PRT
    alphabet = PRT_ALPHABET

    def __repr__(self) -> str:
        return f"{super().__repr__()} N->C"
