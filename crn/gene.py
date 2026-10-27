from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Self, Sequence

from crn.constants import GeneticCode, TRANSCRIPT_PRIORITY
from crn.features import Features
from crn.seq import DNA, RNA, PRT


class Gene:
    """Bir gene ait dizileri bir arada tutan kapsayıcı.

    ::

        Gene (ör. MEFV)
         ├── dna          : DNA            genomik dizi (tek)
         ├── transcripts  : {id: RNA}      transkript varyantları
         ├── proteins     : {id: PRT}      proteinler
         └── features     : Features       gene_id, locus ...

    Transkriptler ve proteinler bağımsız eklenir; aralarındaki eşleşme
    ``isoforms()`` ile kurulur. Veri kaynakları (Ncbi, Gencode) her proteinin
    ``transcript_id`` özelliğine kaynağın bildirdiği transkripti yazar.
    """

    def __init__(self, name: str, dna: Optional[DNA] = None):
        self.name = name
        self._dna: Optional[DNA] = None
        self._transcripts: Dict[str, RNA] = {}
        self._proteins: Dict[str, PRT] = {}
        self._features = Features()
        if dna is not None:
            self.dna = dna

    @property
    def dna(self) -> Optional[DNA]:
        """Genin genomik DNA'sı."""
        return self._dna

    @dna.setter
    def dna(self, value: DNA) -> None:
        if not isinstance(value, DNA):
            raise TypeError(f"dna bir DNA olmalı, {type(value).__name__} verildi")
        self._dna = value

    @property
    def features(self) -> Features:
        return self._features

    @features.setter
    def features(self, value: Features) -> None:
        if not isinstance(value, Features):
            raise TypeError(f"features bir Features olmalı, {type(value).__name__} verildi")
        self._features = value

    # --- Transkriptler ve proteinler -------------------------------------------

    @staticmethod
    def _add(items: Dict, item, item_id: Optional[str], kind: type, label: str, gene: str):
        if not isinstance(item, kind):
            raise TypeError(f"{label} bir {kind.__name__} olmalı, {type(item).__name__} verildi")
        item_id = item_id or item.name
        if not item_id:
            raise ValueError(f"{label.capitalize()} kimliği boş olamaz")
        if item_id in items:
            raise ValueError(f"{gene}: '{item_id}' {label}i zaten var")
        items[item_id] = item
        return item

    def add_transcript(self, rna: RNA, transcript_id: Optional[str] = None) -> RNA:
        """RNA'yı transkript olarak ekler (kimlik verilmezse RNA'nın adı) ve döndürür.

        Raises:
            TypeError: rna bir RNA değilse.
            ValueError: Kimlik boşsa ya da zaten varsa.
        """
        return self._add(self._transcripts, rna, transcript_id, RNA, "transkript", self.name)

    def add_protein(self, prt: PRT, protein_id: Optional[str] = None) -> PRT:
        """Proteini ekler (kimlik verilmezse proteinin adı) ve döndürür.

        Raises:
            TypeError: prt bir PRT değilse.
            ValueError: Kimlik boşsa ya da zaten varsa.
        """
        return self._add(self._proteins, prt, protein_id, PRT, "protein", self.name)

    def transcript(self, transcript_id: str) -> RNA:
        """Raises: KeyError: Bu kimlikte transkript yoksa."""
        if transcript_id not in self._transcripts:
            raise KeyError(f"{self.name}: '{transcript_id}' transkripti yok")
        return self._transcripts[transcript_id]

    def protein(self, protein_id: str) -> PRT:
        """Raises: KeyError: Bu kimlikte protein yoksa."""
        if protein_id not in self._proteins:
            raise KeyError(f"{self.name}: '{protein_id}' proteini yok")
        return self._proteins[protein_id]

    def remove_transcript(self, transcript_id: str) -> bool:
        return self._transcripts.pop(transcript_id, None) is not None

    def remove_protein(self, protein_id: str) -> bool:
        return self._proteins.pop(protein_id, None) is not None

    @property
    def transcripts(self) -> Dict[str, RNA]:
        """{kimlik: RNA}, ekleme sırasıyla; sözlük bir kopyadır."""
        return dict(self._transcripts)

    @property
    def proteins(self) -> Dict[str, PRT]:
        """{kimlik: PRT}, ekleme sırasıyla; sözlük bir kopyadır."""
        return dict(self._proteins)

    # --- İzoformlar -------------------------------------------------------------

    def isoforms(self, genetic_code: GeneticCode = GeneticCode.STANDARD,
                 unique_cds: bool = True,
                 priority: Sequence[str] = TRANSCRIPT_PRIORITY) -> List[Gene]:
        """Geni izoformlarına ayırır: her biri bir DNA, bir transkript, bir protein.

        Her transkript için çevirisi bir proteinle birebir aynı olan ORF aranır
        (``RNA.match_orf``); önce kaynağın o transkripte bağladığı protein
        (``transcript_id`` özelliği), sonra diğerleri denenir. Eşleşen her
        transkript proteiniyle birlikte yeni bir gene konur; bu gen değişmez.

        Yeni genin adı ``<gen>_<transkript>``; DNA, transkript ve protein
        kopyadır. features genin kopyasıdır, ayrıca ``isoform_of``,
        ``transcript_id`` ve ``protein_id`` yazılır. Transkript kırpılmaz; ORF
        konumu transkriptin ``orf_start``/``orf_end`` özelliklerine yazılır.

        ``unique_cds`` True ise aynı CDS'i veren izoformlardan yalnızca biri
        tutulur: kimliği ``priority`` içinde önce gelen önekle başlayan (önce
        NM_, sonra XM_), eşitlikte ilk eklenen. Çıkarılanlar ve eşleşmeyen
        transkript/proteinler ``UserWarning`` ile bildirilir.

        Raises:
            ValueError: Genin genomik DNA'sı yoksa.
        """
        if self._dna is None:
            raise ValueError(f"{self.name}: genomik DNA yok; izoformlar ayrılamaz")

        isoforms, cds, unmatched, used = [], [], [], set()
        for tid, rna in self._transcripts.items():
            linked = [pid for pid, p in self._proteins.items() if p.features.first("transcript_id") == tid]
            candidates = linked + [pid for pid in self._proteins if pid not in linked]
            for pid in candidates:
                prt = self._proteins[pid]
                orf = rna.match_orf(prt, genetic_code)
                if orf is None:
                    continue
                transcript = rna[:]
                transcript.features.set("orf_start", orf.features.first("orf_start"))
                transcript.features.set("orf_end", orf.features.first("orf_end"))
                isoform = Gene(f"{self.name}_{tid}", self._dna[:])
                isoform.features = self._features.copy()
                isoform.features.set("isoform_of", self.name)
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

        unused = [pid for pid in self._proteins if pid not in used]
        details = []
        if unmatched:
            details.append(f"proteini olmayan transkriptler: {', '.join(unmatched)}")
        if unused:
            details.append(f"transkripti olmayan proteinler: {', '.join(unused)}")
        if details:
            warnings.warn(f"{self.name}: {'; '.join(details)}", UserWarning, stacklevel=2)
        if unique_cds:
            isoforms = self._unique_cds(isoforms, cds, priority)
        return isoforms

    def _unique_cds(self, isoforms: List[Gene], cds: List[str],
                    priority: Sequence[str]) -> List[Gene]:
        """Aynı CDS'i veren izoformlardan önceliği en yüksek olanı tutar; sırayı korur."""
        def tid(index: int) -> str:
            return isoforms[index].features.first("transcript_id")

        def rank(index: int) -> tuple:
            prefix = next((i for i, p in enumerate(priority) if tid(index).startswith(p)), len(priority))
            return prefix, index

        groups: Dict[str, List[int]] = {}
        for index, seq in enumerate(cds):
            groups.setdefault(seq, []).append(index)
        kept, dropped = set(), []
        for indices in groups.values():
            best = min(indices, key=rank)
            kept.add(best)
            others = [tid(i) for i in indices if i != best]
            if others:
                dropped.append(f"{', '.join(others)} -> {tid(best)}")
        if dropped:
            warnings.warn(f"{self.name}: aynı CDS'i veren izoformlar çıkarıldı "
                          f"({'; '.join(dropped)})", UserWarning, stacklevel=3)
        return [isoform for i, isoform in enumerate(isoforms) if i in kept]

    # --- Dosya ------------------------------------------------------------------

    def to_crn(self, file_name: str | Path) -> None:
        """Geni tüm içeriğiyle CRN ikili dosyasına yazar (biçim: ``crn.fileio``)."""
        from crn.fileio import write_gene_crn
        write_gene_crn(self, file_name)

    @classmethod
    def from_crn(cls, file_name: str | Path) -> Self:
        """``to_crn`` ile yazılmış gen dosyasını okur."""
        from crn.fileio import read_gene_crn
        return read_gene_crn(file_name)

    def __repr__(self) -> str:
        dna = f"DNA <{len(self._dna)}>" if self._dna is not None else "DNA yok"
        return (f"{self.name}: {dna}, {len(self._transcripts)} transkript, "
                f"{len(self._proteins)} protein")
