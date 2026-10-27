from __future__ import annotations

import struct
import warnings
from pathlib import Path
from typing import (Dict, List, Optional, Self, Sequence)

from crn.io.crn_io_const import (GeneticCode, TRANSCRIPT_PRIORITY)
from crn.io.crn_io_features import Features
from crn.io.crn_io_dna import DNA
from crn.io.crn_io_rna import RNA
from crn.io.crn_io_prt import PRT
from crn.io.crn_io_seqs import (_crn_pack_str, _crn_pack_features, _crn_write, _crn_open)


class Gene:
    """Bir gene ait dizileri bir arada tutan kapsayıcı.

    ::

        Gene (ör. MEFV)
         ├── dna          : DNA            genomik dizi (tek)
         ├── transcripts  : {id: RNA}      istenen sayıda RNA (ör. transkript varyantları)
         ├── proteins     : {id: PRT}      istenen sayıda protein
         └── features     : Features       gene ait diğer bilgiler (gene_id, locus ...)

    Transkriptler ve proteinler birbirinden bağımsız eklenir; hangi RNA'nın
    hangi proteini kodladığı gibi aralarındaki ilişkiler kullanıcının
    sorumluluğundadır (gerekirse features ile tutulabilir).
    """

    def __init__(self, name: str, dna: Optional[DNA] = None):
        self.__name = name
        self.__dna: Optional[DNA] = None
        self.__transcripts: Dict[str, RNA] = {}
        self.__proteins: Dict[str, PRT] = {}
        self.__features = Features()
        if dna is not None:
            self.dna = dna

    # --- Genomik DNA ----------------------------------------------------------

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, value: str) -> None:
        self.__name = value

    @property
    def dna(self) -> Optional[DNA]:
        """Genin genomik DNA'sı (tek)."""
        return self.__dna

    @dna.setter
    def dna(self, value: DNA) -> None:
        if not isinstance(value, DNA):
            raise TypeError(f"dna bir DNA olmalı, {type(value).__name__} verildi")
        self.__dna = value

    @property
    def features(self) -> Features:
        return self.__features

    @features.setter
    def features(self, value: Features) -> None:
        if not isinstance(value, Features):
            raise TypeError(f"features bir Features olmalı, {type(value).__name__} verildi")
        self.__features = value

    # --- Transkriptler ----------------------------------------------------------

    def add_transcript(self, rna: RNA, transcript_id: Optional[str] = None) -> RNA:
        """Bir RNA'yı transkript olarak ekler ve onu döndürür.

        ``transcript_id`` verilmezse RNA'nın adı kullanılır.

        Raises:
            TypeError: rna bir RNA değilse.
            ValueError: Kimlik boşsa ya da bu kimlikte bir transkript zaten varsa.
        """
        if not isinstance(rna, RNA):
            raise TypeError(f"transkript bir RNA olmalı, {type(rna).__name__} verildi")
        transcript_id = transcript_id or rna.name
        if not transcript_id:
            raise ValueError("Transkript kimliği boş olamaz")
        if transcript_id in self.__transcripts:
            raise ValueError(f"{self.__name}: '{transcript_id}' transkripti zaten var")
        self.__transcripts[transcript_id] = rna
        return rna

    def transcript(self, transcript_id: str) -> RNA:
        """Kimliği verilen transkripti döndürür.

        Raises:
            KeyError: Bu kimlikte transkript yoksa.
        """
        if transcript_id not in self.__transcripts:
            raise KeyError(f"{self.__name}: '{transcript_id}' transkripti yok")
        return self.__transcripts[transcript_id]

    def remove_transcript(self, transcript_id: str) -> bool:
        """Transkripti siler.

        Returns:
            Transkript varsa ve silindiyse ``True``, yoksa ``False``.
        """
        if transcript_id not in self.__transcripts:
            return False
        del self.__transcripts[transcript_id]
        return True

    @property
    def transcripts(self) -> Dict[str, RNA]:
        """{transkript kimliği: RNA}; ekleme sırasıyla. Sözlük bir kopyadır."""
        return dict(self.__transcripts)

    # --- Proteinler -----------------------------------------------------------

    def add_protein(self, prt: PRT, protein_id: Optional[str] = None) -> PRT:
        """Bir proteini ekler ve onu döndürür.

        ``protein_id`` verilmezse proteinin adı kullanılır.

        Raises:
            TypeError: prt bir PRT değilse.
            ValueError: Kimlik boşsa ya da bu kimlikte bir protein zaten varsa.
        """
        if not isinstance(prt, PRT):
            raise TypeError(f"protein bir PRT olmalı, {type(prt).__name__} verildi")
        protein_id = protein_id or prt.name
        if not protein_id:
            raise ValueError("Protein kimliği boş olamaz")
        if protein_id in self.__proteins:
            raise ValueError(f"{self.__name}: '{protein_id}' proteini zaten var")
        self.__proteins[protein_id] = prt
        return prt

    def protein(self, protein_id: str) -> PRT:
        """Kimliği verilen proteini döndürür.

        Raises:
            KeyError: Bu kimlikte protein yoksa.
        """
        if protein_id not in self.__proteins:
            raise KeyError(f"{self.__name}: '{protein_id}' proteini yok")
        return self.__proteins[protein_id]

    def remove_protein(self, protein_id: str) -> bool:
        """Proteini siler.

        Returns:
            Protein varsa ve silindiyse ``True``, yoksa ``False``.
        """
        if protein_id not in self.__proteins:
            return False
        del self.__proteins[protein_id]
        return True

    @property
    def proteins(self) -> Dict[str, PRT]:
        """{protein kimliği: PRT}; ekleme sırasıyla. Sözlük bir kopyadır."""
        return dict(self.__proteins)

    # --- İzoformlar -------------------------------------------------------------

    def match_transcripts_to_proteins(self, genetic_code: GeneticCode = GeneticCode.STANDARD,
                                      unique_cds: bool = True,
                                      priority: Sequence[str] = TRANSCRIPT_PRIORITY
                                      ) -> List[Gene]:
        """Geni izoformlarına ayırır: her biri bir DNA, bir transkript, bir protein.

        Her transkript için ORF'leri proteinlerle karşılaştırır
        (``RNA.match_orf``). Çevirisi bir proteinle birebir aynı olan ORF'si
        bulunan her transkript, o proteinle birlikte yeni bir gene konur.
        Bu gen değiştirilmez.

        Aynı proteini kodlayan birden çok transkriptin (ör. yalnızca UTR'leri
        farklı) CDS'i çoğu zaman birebir aynıdır. ``unique_cds`` True ise aynı
        CDS'i veren izoformlardan yalnızca biri tutulur: kimliği ``priority``
        içindeki önceliği en yüksek önekle başlayan (önce NM_, sonra XM_),
        eşitlikte ilk eklenen. Çıkarılanlar ``UserWarning`` ile bildirilir.
        False ise hepsi ayrı izoformlar olarak döner.

        Yeni genlerin adı ``<gen>_<transkript>`` olur; genomik DNA,
        transkript ve protein birer kopyadır. Features genin features
        kopyasıdır; ayrıca ``isoform_of`` (bu genin adı), ``transcript_id`` ve
        ``protein_id`` yazılır. Transkriptin kendisi kırpılmaz; ORF'nin konumu
        transkriptin ``orf_start`` ve ``orf_end`` özelliklerine yazılır.

        Hiçbir proteinle örtüşmeyen transkriptler (ör. kodlamayan RNA'lar) ve
        hiçbir transkripte eşlenmeyen proteinler ``UserWarning`` ile bildirilir.

        Args:
            genetic_code: ORF arama ve translasyon için genetik kod.
            unique_cds: Aynı CDS'i veren izoformlardan yalnızca biri tutulsun mu.
            priority: Tutulacak izoformu seçerken transkript kimliği önekleri,
                öncelik sırasıyla.

        Returns:
            Transkriptlerin ekleme sırasıyla izoform genleri.

        Raises:
            ValueError: Genin genomik DNA'sı yoksa.
        """
        if self.__dna is None:
            raise ValueError(f"{self.__name}: genomik DNA yok; izoformlar ayrılamaz")

        isoforms, cds, unmatched, used = [], [], [], set()
        for tid, rna in self.__transcripts.items():
            for pid, prt in self.__proteins.items():
                orf = rna.match_orf(prt, genetic_code)
                if orf is None:
                    continue
                transcript = rna[:]
                transcript.features.set("orf_start", orf.features.first("orf_start"))
                transcript.features.set("orf_end", orf.features.first("orf_end"))

                isoform = Gene(f"{self.__name}_{tid}", self.__dna[:])
                isoform.features = self.__features.copy()
                isoform.features.set("isoform_of", self.__name)
                isoform.features.set("transcript_id", tid)
                isoform.features.set("protein_id", pid)
                isoform.add_transcript(transcript, tid)
                isoform.add_protein(prt[:], pid)
                isoforms.append(isoform)
                cds.append(orf.seq)
                used.add(pid)
                break
            else:
                unmatched.append(tid)

        unused = [pid for pid in self.__proteins if pid not in used]
        if unmatched or unused:
            details = []
            if unmatched:
                details.append(f"proteini olmayan transkriptler: {', '.join(unmatched)}")
            if unused:
                details.append(f"transkripti olmayan proteinler: {', '.join(unused)}")
            warnings.warn(f"{self.__name}: {'; '.join(details)}", UserWarning, stacklevel=2)
        if unique_cds:
            isoforms = self._unique_cds(isoforms, cds, priority)
        return isoforms

    def _unique_cds(self, isoforms: List[Gene], cds: List[str],
                    priority: Sequence[str]) -> List[Gene]:
        """Aynı CDS'i veren izoformlardan önceliği en yüksek olanı tutar; sırayı korur."""
        def rank(index: int) -> tuple:
            tid = isoforms[index].features.first("transcript_id")
            prefix = next((i for i, p in enumerate(priority) if tid.startswith(p)), len(priority))
            return prefix, index

        groups: Dict[str, List[int]] = {}
        for index, seq in enumerate(cds):
            groups.setdefault(seq, []).append(index)
        kept, dropped = set(), []
        for indices in groups.values():
            best = min(indices, key=rank)
            kept.add(best)
            others = [isoforms[i].features.first("transcript_id") for i in indices if i != best]
            if others:
                dropped.append(f"{', '.join(others)} -> "
                               f"{isoforms[best].features.first('transcript_id')}")
        if dropped:
            warnings.warn(f"{self.__name}: aynı CDS'i veren izoformlar çıkarıldı "
                          f"({'; '.join(dropped)})", UserWarning, stacklevel=3)
        return [isoform for i, isoform in enumerate(isoforms) if i in kept]

    # --- Kaydetme --------------------------------------------------------------

    def to_crn(self, file_name: str | Path) -> None:
        """Geni tüm içeriğiyle CRN ikili dosyasına yazar.

        Yazılanlar: ad, features, genomik DNA, transkriptler ve proteinler
        (kimlikleriyle). Her dizi kendi tipi, adı, root alanları ve
        features'ıyla birlikte yazılır. pickle kullanılmaz; biçim
        ``crn_io_seqs`` modülünün başındaki açıklamada tanımlıdır. Dosya varsa
        üzerine yazılır.

        Raises:
            TypeError: Desteklenmeyen tipte bir feature değeri varsa. Bu durumda
                dosyaya hiçbir şey yazılmaz.
        """
        parts = [
            _crn_pack_str("GENE"),
            _crn_pack_str(self.__name),
            _crn_pack_features(self.__features),
        ]
        if self.__dna is None:
            parts.append(struct.pack("<B", 0))
        else:
            parts.append(struct.pack("<B", 1))
            parts.append(self.__dna._crn_record())

        for items in (self.__transcripts, self.__proteins):
            parts.append(struct.pack("<I", len(items)))
            for item_id, seq in items.items():
                parts.append(_crn_pack_str(item_id))
                parts.append(seq._crn_record())

        _crn_write(file_name, b"".join(parts))

    @classmethod
    def from_crn(cls, file_name: str | Path) -> Self:
        """``to_crn`` ile yazılmış gen dosyasını okur.

        Ad, features, genomik DNA, transkriptler ve proteinler kendi
        sınıflarıyla (DNA, RNA, PRT), kimlikleri, root alanları ve
        features'larıyla geri yüklenir; ekleme sırası korunur.

        Raises:
            FileNotFoundError: Dosya yoksa.
            ValueError: Dosya CRN değilse, bozuksa, sürümü desteklenmiyorsa ya
                da bir gen dosyası değilse (ör. tek bir dizinin CRN dosyası).
        """
        reader = _crn_open(file_name)
        kind = reader.text()
        if kind != "GENE":
            raise ValueError(f"{reader.source}: gen dosyası değil (dosya tipi: {kind})")

        gene = cls(reader.text())
        gene.features = reader.features()
        if reader.uint("<B"):
            gene.dna = DNA._crn_read(reader)

        for _ in range(reader.uint("<I")):
            transcript_id = reader.text()
            gene.add_transcript(RNA._crn_read(reader), transcript_id)
        for _ in range(reader.uint("<I")):
            protein_id = reader.text()
            gene.add_protein(PRT._crn_read(reader), protein_id)

        reader.finish()
        return gene

    def __repr__(self) -> str:
        dna = f"DNA <{len(self.__dna)}>" if self.__dna is not None else "DNA yok"
        return (f"{self.__name}: {dna}, {len(self.__transcripts)} transkript, "
                f"{len(self.__proteins)} protein")
