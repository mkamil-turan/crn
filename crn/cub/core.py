"""Kodon kullanım yanlılığı (codon usage bias) — A. Ham sayımlar.

``CUB`` genleri hesaba hazırlar ve A bölümündeki ham sayımları verir; diğer
bölümler (``Composition``, ``Synonymous`` ...) bir ``CUB`` alır ve buradaki
sayımları kullanır.

Girdi bir ``Gene`` listesidir (ya da ``Gene.to_crn`` ile yazılmış ``.crn``
yolları); tek gen için ``[gene]``. Her gende **tam olarak bir genomik DNA,
bir transkript ve bir protein** bulunmalıdır; çok izoformlu genler
``Gene.isoforms()`` ile ayrılır::

    isoforms = ncbi.gene("MEFV").isoforms()
    cub = CUB(isoforms)             # MEFV_NM_000243.3, MEFV_NM_001198536.2

Kodonlar yalnızca okuma çerçevesinde anlamlı olduğundan, transkriptteki
çevirisi proteinle birebir aynı ORF bulunur ve hesap bu CDS üzerinden yapılır.
Verilen genler değiştirilmez; ``cub.genes`` transkripti CDS'e indirgenmiş
kopyaları verir (genomik DNA ve protein nesneleri paylaşılır). Örtüşen ORF
yoksa ``ValueError`` verilir.

Kodon sayımları kurulumda gen başına bir kez yapılır; her ölçüm bu sayımlardan
(64 kodon) hesaplandığı için genom ölçeğinde (~20 000 gen) de hızlıdır.

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
hesaplar. Stop kodonları ``codon_counts`` içinde sayılır, ancak L'ye, amino
asit sayılarına ve frekanslara girmez. Belirsiz baz (N, R ...) içeren kodonlar
hiçbir sayıma girmez.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from crn.constants import CODONS, AMINO_ACIDS, GeneticCode
from crn.gene import Gene


class CubWarning(UserWarning):
    """CUB hesabında bir verinin dışarıda bırakıldığını bildirir."""


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


class CUB:
    """Genlerin kodon kullanımı.

    Args:
        genes: ``Gene`` nesneleri ya da ``.crn`` gen dosyalarının yolları
            (liste). Her gende bir genomik DNA, bir transkript, bir protein.
        genetic_code: ORF arama ve kodon -> amino asit çevirisi için kod.

    Raises:
        TypeError: genes liste değilse ya da elemanı Gene/.crn yolu değilse.
        ValueError: genes boşsa, iki genin adı aynıysa, bir gen kurala
            uymuyorsa ya da transkriptte proteinle örtüşen ORF yoksa.

    İki CUB ``+`` ile birleştirilebilir (ör. ``high + low``); sayımlar
    yeniden yapılmaz.
    """

    def __init__(self, genes: List[Gene | str | Path],
                 genetic_code: GeneticCode = GeneticCode.STANDARD):
        if not isinstance(genes, list):
            raise TypeError(f"genes bir liste olmalı, {type(genes).__name__} verildi; "
                            f"tek gen için genes=[gene] yazın")
        if not genes:
            raise ValueError("genes boş olamaz")

        loaded = [_load_gene(item) for item in genes]
        names = [gene.name for gene in loaded]
        duplicates = sorted({n for n in names if names.count(n) > 1})
        if duplicates:
            raise ValueError(f"Aynı adlı genler var: {', '.join(duplicates)}")
        problems = {gene.name: p for gene in loaded if (p := _gene_problems(gene))}
        if problems:
            lines = "\n".join(f"  {name}: {'; '.join(p)}" for name, p in problems.items())
            raise ValueError("Her gende bir genomik DNA, bir transkript ve bir protein olmalı "
                             f"(çok izoformlu genler için Gene.isoforms()); "
                             f"hesap durduruldu:\n{lines}")

        self._code = genetic_code
        self._genes: Dict[str, Gene] = {}
        self._seqs: Dict[str, str] = {}
        self._counts: Dict[str, Counter] = {}
        for gene in loaded:
            fitted = self._fit_to_protein(gene)
            seq = next(iter(fitted.transcripts.values())).seq_5to3()
            self._genes[gene.name] = fitted
            self._seqs[gene.name] = seq
            self._counts[gene.name] = Counter(self._split(seq))
        self._finish()

    def _finish(self) -> None:
        self._pooled = sum(self._counts.values(), Counter())
        self._pairs: Dict[Optional[str], Counter] = {}

    def _fit_to_protein(self, gene: Gene) -> Gene:
        """Transkripti proteinle örtüşen ORF'ye indirgenmiş bir gen kopyası.

        Raises:
            ValueError: Transkriptte proteinle örtüşen ORF yoksa.
        """
        (tid, rna), = gene.transcripts.items()
        (pid, prt), = gene.proteins.items()
        orf = rna.match_orf(prt, self._code)
        if orf is None:
            raise ValueError(f"{gene.name}: {tid} transkriptinde {pid} proteiniyle "
                             f"örtüşen ORF yok; hesap durduruldu")
        fitted = Gene(gene.name, gene.dna)
        fitted.features = gene.features.copy()
        fitted.add_transcript(rna if orf.seq == rna.seq else orf, tid)
        fitted.add_protein(prt, pid)
        return fitted

    @staticmethod
    def _split(seq: str) -> Iterator[str]:
        return (seq[i:i + 3] for i in range(0, len(seq) - 2, 3))

    def __add__(self, other: CUB) -> CUB:
        if not isinstance(other, CUB):
            return NotImplemented
        if other._code is not self._code:
            raise ValueError("Farklı genetik kodlu CUB'lar birleştirilemez")
        duplicates = sorted(set(self._genes) & set(other._genes))
        if duplicates:
            raise ValueError(f"Aynı adlı genler var: {', '.join(duplicates)}")
        result = CUB.__new__(CUB)
        result._code = self._code
        result._genes = {**self._genes, **other._genes}
        result._seqs = {**self._seqs, **other._seqs}
        result._counts = {**self._counts, **other._counts}
        result._finish()
        return result

    # --- Özellikler -----------------------------------------------------------

    @property
    def genes(self) -> List[Gene]:
        """Hesaba giren genler (transkriptleri CDS'e indirgenmiş kopyalar)."""
        return list(self._genes.values())

    @property
    def names(self) -> List[str]:
        """Gen adları, ekleme sırasıyla."""
        return list(self._genes)

    @property
    def genetic_code(self) -> GeneticCode:
        return self._code

    def _check(self, gene: str) -> str:
        if gene not in self._genes:
            raise KeyError(f"CUB: '{gene}' geni yok")
        return gene

    def _selected(self, gene: Optional[str]) -> List[str]:
        """``gene`` None ise tüm gen adları, değilse yalnızca o (yoksa KeyError)."""
        return self.names if gene is None else [self._check(gene)]

    def _counter(self, gene: Optional[str]) -> Counter:
        return self._pooled if gene is None else self._counts[self._check(gene)]

    def cds(self, gene: str) -> str:
        """Genin hesaba giren CDS dizisi (RNA, 5'->3', stop dahil)."""
        return self._seqs[self._check(gene)]

    def codons(self, gene: Optional[str] = None) -> List[str]:
        """Kodon listesi (belirsiz bazlılar dahil), okuma sırasıyla; None ise genler art arda."""
        return [codon for name in self._selected(gene) for codon in self._split(self._seqs[name])]

    # --- Kodon kümeleri (genetik koda göre) ---------------------------------------

    def sense_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden amino asit kodlayanlar (stop hariç)."""
        stops = self._code.stops
        return [codon for codon in codons if codon not in stops]

    def codon_families(self, amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, List[str]]:
        """Amino asit -> onu kodlayan kodonlar (eş anlamlı kodon ailesi), ``CODONS`` sırasıyla."""
        return {aa: list(self._code.families.get(aa, ())) for aa in amino_acids}

    def synonymous_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden eş anlamlı kodonlar (tek kodonlu amino asitler ve stop hariç)."""
        return [codon for codon in codons if codon in self._code.synonymous_codons]

    def fourfold_codons(self, codons: Sequence[str] = CODONS) -> List[str]:
        """``codons`` içinden dört katlı dejenere kodonlar (3. konumu tümüyle eş anlamlı)."""
        return [codon for codon in codons if codon in self._code.fourfold_codons]

    # --- A. Ham sayımlar ------------------------------------------------------

    def codon_counts(self, gene: Optional[str] = None,
                     codons: Sequence[str] = CODONS) -> Dict[str, int]:
        """1. Kodon sayısı: ``codons`` içindeki her kodonun kaç kez geçtiği (stoplar dahil)."""
        counter = self._counter(gene)
        return {codon: counter.get(codon, 0) for codon in codons}

    def total_codons(self, gene: Optional[str] = None) -> int:
        """2. Toplam kodon sayısı (L): stop kodonları hariç kodon sayısı."""
        counter = self._counter(gene)
        return sum(counter.get(codon, 0) for codon in self._code.sense_codons)

    def amino_acid_counts(self, gene: Optional[str] = None,
                          amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, int]:
        """3. Amino asit sayısı: aynı amino asidi kodlayan kodonların toplamı (stop hariç)."""
        counts = dict.fromkeys(amino_acids, 0)
        table = self._code.table
        for codon, n in self._counter(gene).items():
            aa = table.get(codon)
            if aa in counts:
                counts[aa] += n
        return counts

    def codon_freqs(self, gene: Optional[str] = None) -> Dict[str, float]:
        """4. Kodon frekansı: N_kodon / L; yalnızca amino asit kodlayan kodonlar (toplam 1)."""
        counts = self.codon_counts(gene, self._code.sense_codons)
        total = sum(counts.values())
        return {codon: (n / total if total else 0.0) for codon, n in counts.items()}

    def codon_per_thousand(self, gene: Optional[str] = None) -> Dict[str, float]:
        """5. Binde frekans: kodon frekansı × 1000."""
        return {codon: freq * 1000 for codon, freq in self.codon_freqs(gene).items()}

    # --- Diğer bölümlerin kullandığı sayımlar ---------------------------------------

    def termini(self, gene: Optional[str] = None) -> List[Tuple[str, str]]:
        """Her proteinin (ilk, son) kalıntısı; son kalıntı stop'tan öncekidir."""
        table = self._code.table
        result = []
        for name in self._selected(gene):
            residues = [aa for codon in self._split(self._seqs[name])
                        if (aa := table.get(codon)) is not None and codon not in self._code.stops]
            if residues:
                result.append((residues[0], residues[-1]))
        return result

    def chain_count(self, gene: Optional[str] = None) -> int:
        """Hesaba giren protein zinciri sayısı: gene None ise gen sayısı, değilse 1."""
        return len(self._selected(gene))

    def codon_pair_counts(self, gene: Optional[str] = None) -> Dict[str, int]:
        """Komşu kodon çiftlerinin sayıları, 6 harflik anahtarla (ör. ``"GCCGAA"``).

        Stop ya da belirsiz baz içeren çiftler sayılmaz. Sonuç önbelleğe alınır.
        """
        key = None if gene is None else self._check(gene)
        if key not in self._pairs:
            if key is None:
                self._pairs[None] = sum((Counter(self.codon_pair_counts(n)) for n in self.names), Counter())
            else:
                sense = set(self._code.sense_codons)
                codons = self.codons(key)
                self._pairs[key] = Counter(a + b for a, b in zip(codons, codons[1:])
                                           if a in sense and b in sense)
        return dict(self._pairs[key])

    # --- Tablo ----------------------------------------------------------------------

    def table(self, **kwargs) -> List[Dict[str, object]]:
        """Her gen için tüm skaler ölçümler tek tabloda (bkz. ``crn.cub.report.summary_table``)."""
        from crn.cub.report import summary_table
        return summary_table(self, **kwargs)

    def __len__(self) -> int:
        return len(self._genes)

    def __repr__(self) -> str:
        return f"CUB: {len(self._genes)} gen, L={self.total_codons()}"
