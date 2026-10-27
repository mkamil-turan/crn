"""Kodon kullanım yanlılığı (codon usage bias) — G. Entropi ve sapma tabanlı ölçümler.

Bir ``CUB`` nesnesinin kodon kullanımının düzensizliğini ve beklenen
kullanımdan sapmasını ölçer. Kodon sayıları ve aileler ``CUB``'dan, kodon
konumlarındaki baz oranları B'deki ``Composition``'dan, kodon oranları C'deki
``Synonymous``'tan gelir.

G. Entropi ve sapma tabanlı ölçümler:

====  ==============================================================  ===========
 No   Ölçüm                                                           Metot
====  ==============================================================  ===========
 28   SCUO: Σ_aa (n_aa / n) × (Hmax − H) / Hmax                       ``scuo``
 29   CDC: 1 − (A·B) / (|A||B|); A gözlenen, B beklenen kodon oranı   ``cdc``
 30   RCDI: Σ (oran_gen / oran_ref) × N_kodon / N                     ``rcdi``
 31   RCBS: [Π (1 + d)]^(1/L) − 1, d = (f_xyz − f1x·f2y·f3z) /        ``rcbs``
      (f1x·f2y·f3z)
====  ==============================================================  ===========

**SCUO** (Wan ve ark. 2004): eş anlamlı her amino asit için kodonlarının
Shannon entropisi H = −Σ p log2 p ve en büyük entropi Hmax = log2(aile
büyüklüğü) hesaplanır; (Hmax − H) / Hmax amino asidin genteki payıyla
ağırlıklandırılıp toplanır. 0 yanlılık yok, 1 her amino asit için tek kodon.

**CDC ve RCBS**: beklenen kodon oranı kodon konumlarındaki baz
oranlarının çarpımıdır, f1(x)·f2(y)·f3(z); stop kodonları iki tarafta da
yoktur ve beklenen oranlar amino asit kodlayan kodonlar üzerinden toplamı 1
olacak biçimde normalize edilir. CDC'de Zhang ve ark. (2012) arka plan bileşimini her konumun GC ve
pürin içeriğiyle modeller ve anlamlılığı bootstrap ile sınar; burada
beklenen bileşim doğrudan konum baz oranlarından kurulur ve bootstrap
yapılmaz. RCBS (Roymondal ve ark. 2009) genin her kodonu için (1 + d)'nin
geometrik ortalamasıdır; 0 yanlılık yok demektir.

**RCDI** (Mueller ve ark. 2006): genin kodon oranının (kodonun kendi amino
asidi içindeki payı) referanstakine (konak) oranı, genin kodonlarıyla
ağırlıklandırılır. Kullanıcının tanımındaki RSCU_gen / RSCU_ref aynı aile
içinde kodon oranlarının oranına eşittir. Yalnızca eş anlamlı kodonlar
sayılır; 1 konakla aynı kullanım, büyük değer konağa uyumsuzluk demektir.
Referansta hiç geçmeyen kodonlar hesaptan çıkarılır ve ``CubWarning`` ile
bildirilir.

Kullanım::

    from crn.cub import CUB, Divergence

    div = Divergence(cub, reference=CUB(genes=host_genes))
    div.scuo(gene="MEFV_NM_000243.3")
    div.rcdi(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

import math
import warnings
from typing import (Dict, Optional, Sequence)

from crn.io.crn_io_const import CODONS
from crn.cub.crn_cub_counts import (CUB, CubWarning)
from crn.cub.crn_cub_composition import Composition
from crn.cub.crn_cub_synonymous import Synonymous


class Divergence:
    """Bir ``CUB`` nesnesinin kodon kullanımındaki düzensizlik ve sapma.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.
        reference: RCDI'deki konak (referans) gen seti; verilmezse ``cub``.

    Raises:
        TypeError, ValueError: ``Synonymous`` ile aynı koşullarda.
    """

    def __init__(self, cub: CUB, reference: Optional[CUB] = None):
        self.__synonymous = Synonymous(cub, reference)
        self.__composition = Composition(cub)

    @property
    def cub(self) -> CUB:
        return self.__synonymous.cub

    @property
    def reference(self) -> CUB:
        return self.__synonymous.reference

    def expected_codon_freqs(self, gene: Optional[str] = None,
                             codons: Sequence[str] = CODONS) -> Dict[str, float]:
        """Konum baz oranlarından beklenen kodon oranları: f1(x)·f2(y)·f3(z).

        Stop kodonları gözlenen kodon oranlarında olmadığı için beklenen
        oranlar da amino asit kodlayan kodonlar üzerinden toplamı 1 olacak
        biçimde normalize edilir.
        """
        f = [self.__composition.nucleotide_freqs(gene, position, codons) for position in (1, 2, 3)]
        products = {codon: f[0][codon[0]] * f[1][codon[1]] * f[2][codon[2]]
                    for codon in self.cub.sense_codons(codons)}
        total = sum(products.values())
        return {codon: (v / total if total else 0.0) for codon, v in products.items()}

    # --- G. Entropi ve sapma -----------------------------------------------------------

    def scuo(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """28. SCUO: amino asit payıyla ağırlıklı normalize entropi farkı, (Hmax − H) / Hmax."""
        counts = self.cub.codon_counts(gene, codons)
        orders: Dict[str, float] = {}
        totals: Dict[str, int] = {}
        for aa, family in self.cub.codon_families().items():
            n = sum(counts.get(codon, 0) for codon in family)
            if len(family) < 2 or not n:
                continue
            entropy = -sum(counts[codon] / n * math.log2(counts[codon] / n)
                           for codon in family if counts.get(codon))
            h_max = math.log2(len(family))
            orders[aa] = (h_max - entropy) / h_max
            totals[aa] = n
        total = sum(totals.values())
        return sum(totals[aa] / total * orders[aa] for aa in orders) if total else 0.0

    def cdc(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """29. CDC: gözlenen ve beklenen kodon oranları arasındaki kosinüs uzaklığı (0-1)."""
        observed = self.cub.codon_freqs(gene, codons)
        expected = self.expected_codon_freqs(gene, codons)
        dot = sum(observed[c] * expected[c] for c in observed)
        norm = math.sqrt(sum(v * v for v in observed.values())) \
            * math.sqrt(sum(v * v for v in expected.values()))
        return 1 - dot / norm if norm else 0.0

    def rcdi(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """30. RCDI: Σ (oran_gen / oran_ref) × N_kodon / N, eş anlamlı kodonlarda.

        Returns:
            RCDI; gende sayılabilecek kodon yoksa 0.
        """
        gene_ratio = self.__synonymous.codon_ratio(gene, codons)
        reference_ratio = Synonymous(self.reference).codon_ratio(None, codons)
        synonymous = self.cub.synonymous_codons(codons)
        counts = self.cub.codon_counts(gene, synonymous)

        missing = [c for c in synonymous if counts[c] and not reference_ratio[c]]
        if missing:
            warnings.warn(f"RCDI: referansta hiç geçmeyen kodonlar hesaptan çıkarıldı: "
                          f"{', '.join(missing)}", CubWarning, stacklevel=2)
        used = [c for c in synonymous if c not in missing]
        total = sum(counts[c] for c in used)
        if not total:
            return 0.0
        return sum(gene_ratio[c] / reference_ratio[c] * counts[c] for c in used if counts[c]) / total

    def rcbs(self, gene: Optional[str] = None, codons: Sequence[str] = CODONS) -> float:
        """31. RCBS: genin kodonları üzerinden (1 + d)'nin geometrik ortalaması − 1."""
        observed = self.cub.codon_freqs(gene, codons)
        expected = self.expected_codon_freqs(gene, codons)
        counts = self.cub.codon_counts(gene, list(observed))
        total = sum(counts.values())
        if not total:
            return 0.0
        log_sum = sum(n * math.log(observed[c] / expected[c]) for c, n in counts.items() if n)
        return math.exp(log_sum / total) - 1

    def __repr__(self) -> str:
        return f"Divergence({self.__synonymous!r})"
