"""Kodon kullanım yanlılığı (codon usage bias) — H. tRNA ve kodon bağlamı.

Bir ``CUB`` nesnesinin kodonlarını tRNA havuzuna ve komşu kodonlarına göre
değerlendirir. Kodon ve amino asit sayıları ``CUB``'dan, kodon konumlarındaki
baz oranları B'deki ``Composition``'dan gelir.

H. tRNA ve kodon bağlamı:

====  ==============================================================  ===========================
 No   Ölçüm                                                           Metot
====  ==============================================================  ===========================
 32   tAI: W = Σ (1 − s)·tGCN; w = W / W_en büyük;                    ``trna_weights``, ``tai``
      tAI = exp(Σ ln w / L)
 33   Dinükleotit ρ: f(XY) / (f(X)·f(Y))                              ``dinucleotide_counts``,
                                                                      ``rho``
 34   CPS: ln [N(AB) / (N(A)·N(B) / (N(X)·N(Y)) × N(XY))];            ``codon_pair_scores``,
      CPB: genin kodon çiftlerinin CPS ortalaması                     ``cpb``
====  ==============================================================  ===========================

**tAI (dos Reis ve ark. 2004).** tGCN, kodonu okuyabilen antikodonun tRNA
gen kopya sayısıdır; ``trna_counts`` (antikodon 5'→3' -> kopya sayısı) ile
verilir ve GtRNAdb FASTA dosyasından ``read_trna_counts`` ile okunabilir.
Kodonun 3. bazını okuyan antikodonlar ve wobble cezaları s,
``crn.tables`` içindeki ``TAI_WOBBLE_PAIRS`` ve ``TAI_S_VALUES``
sabitleridir (genomdaki A34 inozin kabul edilir). ``kingdom="prokaryote"``
ise AUA kodonu lizidinli tRNA-Ile CAU'dan da (L:A) okunur. Met kodonları
(başlatıcı tRNA ile karışmasın diye) ve stop kodonları hesaba girmez. W'si 0
olan kodonların w'si diğer w'lerin geometrik ortalaması alınır.

**Dinükleotit ρ (Karlin ve Burge 1995).** ``frame`` None ise kodlayan dizinin
tamamındaki dinükleotitler ve baz oranları kullanılır. ``"12"`` ve ``"23"``
kodon içindeki konum çiftleri, ``"31"`` kodon sınırını aşan çifttir (bir
kodonun 3. ve sonrakinin 1. bazı); bu durumlarda f(X) ve f(Y) o konumların
baz oranlarıdır. ρ < 0,78 az, > 1,23 çok temsil edilmiş kabul edilir; CpG
için ``rho()["CG"]``.

**CPS ve CPB (Coleman ve ark. 2008).** Sayımlar referans gen setinden
(genellikle genomun tüm ORF'leri) alınır: A ve B komşu kodonlar, X ve Y
onların amino asitleridir. Pozitif CPS çiftin beklenenden çok, negatif
beklenenden az kullanıldığını gösterir. CPB genin her kodon çiftinin CPS'sinin
ortalamasıdır; referansta hiç geçmeyen çiftler hesaptan çıkarılır ve
``CubWarning`` ile bildirilir. Stop içeren çiftler sayılmaz.

Kullanım::

    from crn.cub import CUB, CodonContext, read_trna_counts

    ctx = CodonContext(cub, reference=CUB(host_genes))
    ctx.tai(gene="MEFV_NM_000243.3", trna_counts=read_trna_counts("hg38-tRNAs.fa"))
    ctx.rho(gene="MEFV_NM_000243.3", frame="31")["CG"]
    ctx.cpb(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

import math
import re
import warnings
from functools import cached_property
from pathlib import Path
from typing import Dict, Mapping, Optional

from crn.constants import RNA_BASES, DNA_TO_RNA, RNA_COMPLEMENT
from crn.tables import TAI_S_VALUES, TAI_WOBBLE_PAIRS
from crn.cub.core import CUB, CubWarning
from crn.cub.composition import Composition
from crn.cub.synonymous import check_cub


# Dinükleotit sayımında kullanılabilecek konum çiftleri (None: dizinin tamamı)
# ve her birinde f(X), f(Y) için kullanılan kodon konumları.
DINUCLEOTIDE_FRAMES = {None: (None, None), "12": (1, 2), "23": (2, 3), "31": (3, 1)}
# tAI'de lizidin (L:A) kuralının uygulandığı canlı grubu.
TAI_KINGDOMS = ("eukaryote", "prokaryote")
# GtRNAdb FASTA başlıkları: ">Homo_sapiens_tRNA-Ala-AGC-1-1 ..." ya da ">... Ala (AGC) 73 bp ...".
_TRNA_HEADER_PATTERNS = (re.compile(r"tRNA-([A-Za-z]+)-([ACGTU]{3})-"),
                         re.compile(r"\s([A-Za-z]+) \(([ACGTU]{3})\)"))
# Kodon okumayan ya da standart dışı tRNA'lar.
_TRNA_SKIPPED = frozenset({"Undet", "Sup", "SeC", "Pseudo"})


def read_trna_counts(file_name: str | Path) -> Dict[str, int]:
    """GtRNAdb FASTA dosyasından antikodon başına tRNA gen kopya sayısı (tGCN).

    Antikodonlar 5'->3' RNA alfabesiyle döner (``"TTC"`` -> ``"UUC"``).
    Belirsiz, baskılayıcı ve selenosistein tRNA'ları sayılmaz.
    """
    counts: Dict[str, int] = {}
    with open(file_name, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            for pattern in _TRNA_HEADER_PATTERNS:
                match = pattern.search(line)
                if match:
                    amino_acid, anticodon = match.groups()
                    if amino_acid not in _TRNA_SKIPPED:
                        anticodon = anticodon.translate(DNA_TO_RNA)
                        counts[anticodon] = counts.get(anticodon, 0) + 1
                    break
    return counts


def _anticodon(codon: str, first: str) -> str:
    """Kodonun 1. ve 2. bazını eşleyen, 34. bazı ``first`` olan antikodon (5'->3')."""
    return first + codon[1::-1].translate(RNA_COMPLEMENT)


class CodonContext:
    """Bir ``CUB`` nesnesinin tRNA havuzuna uyumu ve kodon bağlamı.

    Metotlar ``gene`` verilmezse tüm genleri birlikte, verilirse yalnızca o
    geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.
        reference: CPS'nin hesaplandığı referans gen seti; verilmezse ``cub``.

    Raises:
        TypeError: cub ya da reference bir CUB değilse.
    """

    def __init__(self, cub: CUB, reference: Optional[CUB] = None):
        self.composition = Composition(cub)
        self.cub = cub
        self.reference = cub if reference is None else reference
        check_cub("reference", self.reference)

    # --- 32. tAI ---------------------------------------------------------------------

    def trna_weights(self, trna_counts: Mapping[str, float],
                     s_values: Mapping[str, float] = TAI_S_VALUES,
                     kingdom: str = "eukaryote") -> Dict[str, float]:
        """32. Her kodonun tAI ağırlığı w = W / W_en büyük, W = Σ (1 − s)·tGCN.

        Args:
            trna_counts: Antikodon (5'->3', T ya da U) -> tRNA gen kopya sayısı.
            s_values: Wobble cezaları; anahtarlar ``TAI_S_VALUES`` ile aynı.
            kingdom: ``"eukaryote"`` ya da ``"prokaryote"`` (AUA için L:A).

        Returns:
            Met ve stop dışındaki kodonların w'si; W'si 0 olanlara diğer w'lerin
            geometrik ortalaması verilir.

        Raises:
            ValueError: kingdom geçersizse ya da hiçbir kodon için tRNA yoksa.
        """
        if kingdom not in TAI_KINGDOMS:
            raise ValueError(f"kingdom {', '.join(TAI_KINGDOMS)} olmalı, {kingdom!r} verildi")
        trna = {anticodon.upper().translate(DNA_TO_RNA): n for anticodon, n in trna_counts.items()}
        table = self.cub.genetic_code.table
        raw: Dict[str, float] = {}
        for codon in self.cub.genetic_code.sense_codons:
            if table[codon] == "M":
                continue
            raw[codon] = sum((1 - s_values[key]) * trna.get(_anticodon(codon, first), 0)
                             for first, key in TAI_WOBBLE_PAIRS[codon[2]])
            if kingdom == "prokaryote" and codon == "AUA" and table[codon] == "I":
                raw[codon] += (1 - s_values["L:A"]) * trna.get("CAU", 0)

        top = max(raw.values(), default=0)
        if not top:
            raise ValueError("tAI: trna_counts hiçbir kodonu okuyamıyor")
        weights = {codon: w / top for codon, w in raw.items()}
        nonzero = [w for w in weights.values() if w > 0]
        mean = math.exp(sum(math.log(w) for w in nonzero) / len(nonzero))
        return {codon: (w if w > 0 else mean) for codon, w in weights.items()}

    def tai(self, gene: Optional[str] = None, trna_counts: Optional[Mapping[str, float]] = None,
            s_values: Mapping[str, float] = TAI_S_VALUES, kingdom: str = "eukaryote") -> float:
        """32. tAI: genin kodonlarının tAI ağırlıklarının geometrik ortalaması.

        Raises:
            ValueError: trna_counts verilmezse (bkz. ``trna_weights``).
        """
        if not trna_counts:
            raise ValueError("trna_counts verilmeli: antikodon -> tRNA gen kopya sayısı "
                             "(read_trna_counts ile GtRNAdb dosyasından okunabilir)")
        weights = self.trna_weights(trna_counts, s_values, kingdom)
        counts = self.cub.codon_counts(gene, list(weights))
        total = sum(counts.values())
        if not total:
            return 0.0
        return math.exp(sum(n * math.log(weights[c]) for c, n in counts.items() if n) / total)

    # --- 33. Dinükleotit ρ -----------------------------------------------------------

    def dinucleotide_counts(self, gene: Optional[str] = None, frame: Optional[str] = None,
                            bases: str = RNA_BASES) -> Dict[str, int]:
        """33. Kodlayan dizideki dinükleotit sayıları (stop kodonları hariç).

        Args:
            frame: None (tüm dizi), ``"12"``, ``"23"`` ya da ``"31"`` (kodon sınırı).

        Raises:
            ValueError: frame geçersizse.
        """
        if frame not in DINUCLEOTIDE_FRAMES:
            raise ValueError(f"frame None, '12', '23' ya da '31' olmalı, {frame!r} verildi")
        sense = set(self.cub.genetic_code.sense_codons)
        counts = {x + y: 0 for x in bases for y in bases}
        names = self.cub.names if gene is None else [gene]
        for name in names:
            reading = self.cub.codons(name)
            for i, codon in enumerate(reading):
                if codon not in sense:
                    continue
                pairs = []
                if frame in (None, "12"):
                    pairs.append(codon[:2])
                if frame in (None, "23"):
                    pairs.append(codon[1:])
                if frame in (None, "31") and i + 1 < len(reading) and reading[i + 1] in sense:
                    pairs.append(codon[2] + reading[i + 1][0])
                for pair in pairs:
                    if pair in counts:
                        counts[pair] += 1
        return counts

    def rho(self, gene: Optional[str] = None, frame: Optional[str] = None,
            bases: str = RNA_BASES) -> Dict[str, float]:
        """33. Dinükleotit göreli bolluğu ρ(XY) = f(XY) / (f(X)·f(Y)); 1 beklenen düzey.

        Returns:
            Dinükleotit -> ρ; bazlardan biri hiç geçmiyorsa ``nan``.
        """
        counts = self.dinucleotide_counts(gene, frame, bases)
        total = sum(counts.values())
        first, second = DINUCLEOTIDE_FRAMES[frame]
        fx = self.composition.nucleotide_freqs(gene, first, bases=bases)
        fy = self.composition.nucleotide_freqs(gene, second, bases=bases)
        result = {}
        for pair, n in counts.items():
            expected = fx[pair[0]] * fy[pair[1]]
            result[pair] = (n / total) / expected if total and expected else math.nan
        return result

    # --- 34. CPS ve CPB --------------------------------------------------------------

    @cached_property
    def _pair_scores(self) -> Dict[str, float]:
        table = self.cub.genetic_code.table
        reference = self.reference
        pairs = reference.codon_pair_counts()
        codon_n = reference.codon_counts()
        aa_n = reference.amino_acid_counts()
        aa_pairs: Dict[str, int] = {}
        for pair, n in pairs.items():
            key = table[pair[:3]] + table[pair[3:]]
            aa_pairs[key] = aa_pairs.get(key, 0) + n
        scores = {}
        for pair, n in pairs.items():
            a, b = pair[:3], pair[3:]
            x, y = table[a], table[b]
            expected = codon_n[a] * codon_n[b] / (aa_n[x] * aa_n[y]) * aa_pairs[x + y]
            scores[pair] = math.log(n / expected)
        return scores

    def codon_pair_scores(self) -> Dict[str, float]:
        """34. Referansta geçen her kodon çiftinin CPS'si (sonuç önbelleğe alınır).

        CPS = ln [N(AB) / (N(A)·N(B) / (N(X)·N(Y)) × N(XY))], tüm sayımlar
        referanstan.
        """
        return dict(self._pair_scores)

    def cpb(self, gene: Optional[str] = None) -> float:
        """34. CPB: genin kodon çiftlerinin CPS ortalaması.

        Referansta hiç geçmeyen çiftler çıkarılır ve ``CubWarning`` ile
        bildirilir. Hesaba girecek çift yoksa ``nan``.
        """
        scores = self._pair_scores
        pairs = self.cub.codon_pair_counts(gene)
        missing = [pair for pair in pairs if pair not in scores]
        if missing:
            warnings.warn(f"CPB: referansta hiç geçmeyen {len(missing)} kodon çifti "
                          f"hesaptan çıkarıldı", CubWarning, stacklevel=2)
        total = sum(n for pair, n in pairs.items() if pair in scores)
        if not total:
            return math.nan
        return sum(scores[pair] * n for pair, n in pairs.items() if pair in scores) / total

    def __repr__(self) -> str:
        return f"CodonContext({self.cub!r}, reference={self.reference!r})"
