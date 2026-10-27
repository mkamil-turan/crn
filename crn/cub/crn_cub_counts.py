"""Kodon kullanım yanlılığı (codon usage bias) — A. Ham sayımlar.

``crn.cub`` paketinin temel modülüdür. ``CUB`` sınıfı genleri hesaba hazırlar
(kodonlara böler) ve A bölümündeki ham sayımları verir. Diğer bölümlerin
modülleri (nükleotit bileşimi, sinonim kodon tercihi ...) bir ``CUB``
nesnesini girdi alır ve buradaki sayımları kullanır.

Girdi her zaman bir ``Gene`` listesidir; tek gen için ``genes=[gene]`` yazılır.
Liste elemanı bir ``Gene`` ya da ``Gene.to_crn`` ile yazılmış bir ``.crn``
dosyasının yolu olabilir.

Her gen **tam olarak bir genomik DNA, bir transkript ve bir protein**
içermelidir. Bunu kullanıcı sağlar; uymayan bir gen varsa ``CUB`` hesaba
başlamadan ``ValueError`` verir. Çok izoformlu bir gen
``Gene.match_transcripts_to_proteins`` ile bu yapıdaki izoform genlerine
ayrılır::

    isoforms = ncbi.gene("MEFV").match_transcripts_to_proteins()
    cub = CUB(genes=isoforms)       # MEFV_NM_000243.3, MEFV_NM_001198536.2

Kodonlar ancak bir okuma çerçevesi içinde ve translasyona uğrayan dizide
anlamlıdır; bu yüzden yalnızca protein kodlayan genler hesaba girer.
Transkriptteki ORF'ler bulunur ve çevirisi proteinle birebir aynı olan ORF
transkriptin yerine konur (transkript mRNA'nın kodlayan bölgesine, yani CDS'e
indirgenir). Gen **yerinde** düzenlenir ve ``CubWarning`` ile bildirilir.
Örtüşen ORF yoksa ``ValueError`` verilir. Kodlamayan RNA'lar (NR_, XR_)
için kodon kullanımı tanımsızdır.

A. Ham sayımlar:

====  ===================================  ==========================
 No   Ölçüm                                Metot
====  ===================================  ==========================
 1    Kodon sayısı (N_kodon)               ``codon_counts``
 2    Toplam kodon sayısı (L, stop hariç)  ``total_codons``
 3    Amino asit sayısı                    ``amino_acid_counts``
 4    Kodon frekansı (N_kodon / L)         ``codon_freqs``
 5    Binde frekans (frekans × 1000)       ``codon_per_thousand``
====  ===================================  ==========================

Metotlar ``gene`` verilmezse tüm genleri birlikte, verilirse yalnızca o geni
hesaplar. Kodon ve amino asit listeleri ``crn.io.crn_io_const`` içindeki
``CODONS`` ve ``AMINO_ACIDS`` sabitlerinden gelir ve parametreyle
değiştirilebilir. Stop kodonları ``codon_counts`` içinde sayılır, ancak L'ye,
amino asit sayılarına ve frekanslara girmez. Belirsiz baz (N, R ...) içeren
kodonlar hiçbir sayıma girmez.

Kullanım::

    from crn.cub import CUB

    cub = CUB(genes=isoforms)
    cub.total_codons(gene="MEFV_NM_000243.3")       # L
    cub.codon_counts()                              # tüm genler birlikte
    cub.codon_per_thousand(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import (Dict, List, Optional, Sequence)

from crn.io import (RNA, Gene, GeneticCode, Strand)
from crn.io.crn_io_const import (CODONS, AMINO_ACIDS, STOP_SYMBOL, RNA_BASES)


class CubWarning(UserWarning):
    """CUB hazırlığında genin değiştirildiğini ya da dizilerin nasıl okunduğunu bildirir."""


def _load_gene(item: Gene | str | Path) -> Gene:
    if isinstance(item, Gene):
        return item
    if isinstance(item, (str, Path)):
        return Gene.from_crn(item)
    raise TypeError(f"genes elemanı bir Gene ya da .crn dosya yolu olmalı, "
                    f"{type(item).__name__} verildi")


def _gene_problems(gene: Gene) -> List[str]:
    """Genin "bir DNA, bir transkript, bir protein" kuralına uymayan yanları."""
    problems = []
    if gene.dna is None:
        problems.append("genomik DNA yok")
    for kind, items in (("transkript", gene.transcripts), ("protein", gene.proteins)):
        if len(items) != 1:
            ids = f" ({', '.join(items)})" if items else ""
            problems.append(f"{len(items)} {kind}{ids}")
    return problems


def _reading_seq(rna: RNA) -> str:
    """RNA'nın ribozomun okuyacağı 5'->3' dizisi (boşluklar atılmış)."""
    seq = rna.seq.replace("-", "")
    return seq[::-1] if rna.strand is Strand.ThreePrime else seq


class CUB:
    """Genlerin kodon kullanımını hesaplar.

    Args:
        genes: ``Gene`` nesneleri ya da ``.crn`` gen dosyalarının yolları;
            **liste olmalı**, tek gen için ``[gene]``. Her gende bir genomik
            DNA, bir transkript ve bir protein bulunmalı; çok izoformlu genler
            için ``Gene.match_transcripts_to_proteins``. Transkript proteinle
            örtüşen ORF'ye indirgenir ve gen yerinde düzenlenir.
        genetic_code: ORF arama ve kodon -> amino asit çevirisi için genetik kod.

    Raises:
        TypeError: genes liste değilse ya da elemanı Gene/.crn yolu değilse.
        ValueError: genes boşsa, iki genin adı aynıysa, bir gen "bir DNA, bir
            transkript, bir protein" kuralına uymuyorsa ya da transkriptte
            proteinle örtüşen ORF yoksa.
    """

    def __init__(self, genes: List[Gene | str | Path],
                 genetic_code: GeneticCode = GeneticCode.STANDARD):
        if not isinstance(genes, list):
            raise TypeError(f"genes bir liste olmalı, {type(genes).__name__} verildi; "
                            f"tek gen için genes=[gene] yazın")
        if not genes:
            raise ValueError("genes boş olamaz")

        self.__genes = [_load_gene(item) for item in genes]
        names = [gene.name for gene in self.__genes]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"Aynı adlı genler var: {', '.join(duplicates)}")

        problems = {gene.name: _gene_problems(gene) for gene in self.__genes}
        problems = {name: p for name, p in problems.items() if p}
        if problems:
            lines = "\n".join(f"  {name}: {'; '.join(p)}" for name, p in problems.items())
            raise ValueError("Her gende bir genomik DNA, bir transkript ve bir protein olmalı "
                             f"(çok izoformlu genler için Gene.match_transcripts_to_proteins); "
                             f"hesap durduruldu:\n{lines}")

        self.__genetic_code = genetic_code

        # {gen adı: [kodon, ...]}
        self.__codons: Dict[str, List[str]] = {}
        for gene in self.__genes:
            self._fit_to_protein(gene)
            self.__codons[gene.name] = self._split_codons(self._transcript(gene))

    # --- Özellikler -----------------------------------------------------------

    @property
    def genes(self) -> List[Gene]:
        """Hesaba giren genler (transkriptleri CDS'e indirgenmiş hâlleri)."""
        return list(self.__genes)

    @property
    def genetic_code(self) -> GeneticCode:
        return self.__genetic_code

    def codons(self, gene: Optional[str] = None) -> List[str]:
        """Kodon listesi (belirsiz bazlı kodonlar dahil), okuma sırasıyla.

        Raises:
            KeyError: Gen yoksa.
        """
        if gene is None:
            return [codon for codons in self.__codons.values() for codon in codons]
        if gene not in self.__codons:
            raise KeyError(f"CUB: '{gene}' geni yok")
        return list(self.__codons[gene])

    # --- Hazırlık ---------------------------------------------------------------

    @staticmethod
    def _transcript(gene: Gene) -> RNA:
        return next(iter(gene.transcripts.values()))

    def _fit_to_protein(self, gene: Gene) -> None:
        """Transkripti proteinle örtüşen ORF'ye indirger; geni yerinde düzenler.

        Raises:
            ValueError: Transkriptte proteinle örtüşen ORF yoksa.
        """
        tid, rna = next(iter(gene.transcripts.items()))
        pid, prt = next(iter(gene.proteins.items()))
        orf = rna.match_orf(prt, self.__genetic_code)
        if orf is None:
            raise ValueError(f"{gene.name}: {tid} transkriptinde {pid} proteiniyle "
                             f"örtüşen ORF yok; hesap durduruldu")
        if orf.seq == rna.seq:
            return
        gene.remove_transcript(tid)
        gene.add_transcript(orf, tid)
        warnings.warn(f"{gene.name}: gen yeniden düzenlendi ({tid}, {pid} proteiniyle "
                      f"örtüşen ORF'ye indirgendi: {orf.features.first('orf_start')}-"
                      f"{orf.features.first('orf_end')}).", CubWarning, stacklevel=3)

    @staticmethod
    def _split_codons(rna: RNA) -> List[str]:
        """CDS'e indirgenmiş RNA'yı ilk nükleotitten başlayarak üçlülere böler."""
        seq = _reading_seq(rna)
        return [seq[i:i + 3] for i in range(0, len(seq) - 2, 3)]

    # --- A. Ham sayımlar ------------------------------------------------------

    def sense_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden bu genetik kodda amino asit kodlayanlar (stop hariç)."""
        stops = self.__genetic_code.stops
        return [codon for codon in codons if codon not in stops]

    def codon_families(self, amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, List[str]]:
        """Her amino asidi kodlayan kodonlar (eş anlamlı kodon ailesi), seçilen genetik koda göre.

        Aileler 64 kodonun tamamından kurulur; kodonlar ``CODONS`` sırasıyla
        dizilir. Yalnızca ``amino_acids`` içindeki amino asitler döner.
        """
        table = self.__genetic_code.table
        families: Dict[str, List[str]] = {aa: [] for aa in amino_acids}
        for codon in self.sense_codons(CODONS):
            if table[codon] in families:
                families[table[codon]].append(codon)
        return families

    def synonymous_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden eş anlamlı kodonlar: amino asidi birden çok kodonla kodlananlar.

        Dejenerelik seçilen genetik kodun 64 kodonu üzerinden belirlenir
        (``codon_families``); tek kodonlu amino asitler (standart kodda Met ve
        Trp) ve stop kodonları dışarıda kalır.
        """
        table = self.__genetic_code.table
        families = self.codon_families()
        return [codon for codon in self.sense_codons(codons)
                if len(families.get(table[codon], ())) > 1]

    def fourfold_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden dört katlı dejenere kodonlar: 3. konumu tümüyle eş anlamlı olanlar.

        xyN kutusundaki dört kodonun hepsi aynı amino asidi kodluyorsa o
        kutunun kodonları döner (standart kodda Ala, Gly, Pro, Thr, Val ile
        Leu CUN, Arg CGN, Ser UCN; 32 kodon).
        """
        table = self.__genetic_code.table
        return [codon for codon in self.sense_codons(codons)
                if len({table[codon[:2] + base] for base in RNA_BASES}) == 1]

    def codon_counts(self, gene: Optional[str] = None,
                     codons: Sequence[str] = CODONS) -> Dict[str, int]:
        """1. Kodon sayısı (N_kodon): ``codons`` içindeki her kodonun kaç kez geçtiği.

        Stop kodonları da sayılır. ``codons`` dışındaki (ör. belirsiz bazlı)
        kodonlar sayılmaz.
        """
        counts = dict.fromkeys(codons, 0)
        for codon in self.codons(gene):
            if codon in counts:
                counts[codon] += 1
        return counts

    def total_codons(self, gene: Optional[str] = None,
                     codons: Sequence[str] = CODONS) -> int:
        """2. Toplam kodon sayısı (L): stop kodonları hariç kodon sayısı."""
        counts = self.codon_counts(gene, codons)
        return sum(counts[codon] for codon in self.sense_codons(codons))

    def amino_acid_counts(self, gene: Optional[str] = None,
                          amino_acids: Sequence[str] = AMINO_ACIDS,
                          codons: Sequence[str] = CODONS) -> Dict[str, int]:
        """3. Amino asit sayısı: aynı amino asidi kodlayan kodonların toplamı (stop hariç).

        Yalnızca ``amino_acids`` içindeki amino asitler döner.
        """
        counts = dict.fromkeys(amino_acids, 0)
        table = self.__genetic_code.table
        for codon, n in self.codon_counts(gene, codons).items():
            if table[codon] in counts:
                counts[table[codon]] += n
        return counts

    def codon_freqs(self, gene: Optional[str] = None,
                    codons: Sequence[str] = CODONS) -> Dict[str, float]:
        """4. Kodon frekansı: N_kodon / L; yalnızca amino asit kodlayan kodonlar.

        Değerlerin toplamı 1'dir (hiç kodon yoksa hepsi 0).
        """
        counts = self.codon_counts(gene, codons)
        total = self.total_codons(gene, codons)
        return {codon: (counts[codon] / total if total else 0.0)
                for codon in self.sense_codons(codons)}

    def codon_per_thousand(self, gene: Optional[str] = None,
                           codons: Sequence[str] = CODONS) -> Dict[str, float]:
        """5. Binde frekans: kodon frekansı × 1000."""
        return {codon: freq * 1000
                for codon, freq in self.codon_freqs(gene, codons).items()}

    def __repr__(self) -> str:
        return f"CUB: {len(self.__genes)} gen, L={self.total_codons()}"
