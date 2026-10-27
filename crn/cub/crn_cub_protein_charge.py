"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, C. Yük ve iyonlaşma.

Bir ``CUB`` nesnesinin amino asit sayımlarından (``CUB.amino_acid_counts``)
proteinin yükünü hesaplar. Uzunluk ve zincir sayısı E-A'daki
``ProteinMass``'tan, zincir uçlarındaki kalıntılar ``CUB.codons``'tan gelir.

E-C. Yük ve iyonlaşma:

====  ==============================================================  ======================
 No   Ölçüm                                                           Metot
====  ==============================================================  ======================
 12   Pozitif (R + K) ve negatif (D + E) yüklü kalıntı sayısı         ``charged_counts``
 13   Kaba net yük (pH 7): (R + K) − (D + E), His ihmal edilir        ``simple_net_charge``
 14   FCR: (R + K + D + E) / N                                        ``fcr``
 15   NCPR: [(R + K) − (D + E)] / N                                   ``ncpr``
 16   Q(pH): Henderson–Hasselbalch                                    ``net_charge``
 17   Teorik pI: Q(pH) = 0 olan pH (ikiye bölme)                      ``isoelectric_point``
 18   Titrasyon eğrisi: pH 0-14 boyunca Q(pH)                         ``titration_curve``
====  ==============================================================  ======================

Q(pH) = Σ_pozitif n / (1 + 10^(pH − pKa)) − Σ_negatif n / (1 + 10^(pKa − pH)).
Pozitif gruplar N-ucu, K, R, H; negatif gruplar C-ucu, D, E, C, Y'dir. Uç
gruplar serbest varsayılır; tüm Cys'ler serbest (indirgenmiş) sayılır.

pKa kümeleri ``crn.io.crn_io_const.PKA_SETS``'ten gelir: ``"EMBOSS"``,
``"Lehninger"``, ``"Bjellqvist"`` (varsayılan; ExPASy Compute pI/Mw).
Bjellqvist'te uç grupların pKa'sı uçtaki kalıntıya göre değişir. Aynı
anahtarlarla kurulmuş bir sözlük de verilebilir.

``gene`` verilmezse tüm genler birlikte hesaplanır: yan zincirler tüm
proteinlerden, uç gruplar her zincir için ayrı ayrı sayılır.

Kullanım::

    from crn.cub import CUB, ProteinCharge

    charge = ProteinCharge(cub)
    charge.isoelectric_point(gene="MEFV_NM_000243.3")
    charge.net_charge(gene="MEFV_NM_000243.3", ph=7.4, pka="EMBOSS")
    charge.titration_curve(gene="MEFV_NM_000243.3", step=0.5)
"""
from __future__ import annotations

from typing import (Dict, List, Mapping, Optional, Tuple)

from crn.io.crn_io_const import (POSITIVE_RESIDUES, NEGATIVE_RESIDUES,
                                 IONIZABLE_POSITIVE, IONIZABLE_NEGATIVE,
                                 PKA_SETS, STOP_SYMBOL)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass


# Varsayılan pKa kümesi ve pI aramasının sınırları, duyarlılığı.
DEFAULT_PKA_SET = "Bjellqvist"
PH_RANGE = (0.0, 14.0)
PI_TOLERANCE = 1e-4


class ProteinCharge:
    """Bir ``CUB`` nesnesindeki proteinlerin yükü ve iyonlaşması.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        self.__mass = ProteinMass(cub)

    @property
    def cub(self) -> CUB:
        return self.__mass.cub

    # --- Yardımcılar -------------------------------------------------------------

    @staticmethod
    def _pka_set(pka: str | Mapping) -> Mapping:
        if isinstance(pka, str):
            if pka not in PKA_SETS:
                raise ValueError(f"pka {', '.join(map(repr, PKA_SETS))} ya da bir sözlük "
                                 f"olmalı, {pka!r} verildi")
            return PKA_SETS[pka]
        return pka

    def termini(self, gene: Optional[str] = None) -> List[Tuple[str, str]]:
        """Her zincirin (ilk, son) kalıntısı; son kalıntı stop'tan önceki kalıntıdır.

        Raises:
            KeyError: Gen yoksa.
        """
        cub = self.cub
        table = cub.genetic_code.table
        names = [g.name for g in cub.genes] if gene is None else [gene]
        result = []
        for name in names:
            residues = [table[c] for c in cub.codons(name) if table.get(c, STOP_SYMBOL) != STOP_SYMBOL]
            if residues:
                result.append((residues[0], residues[-1]))
        return result

    # --- E-C. Yük ve iyonlaşma -----------------------------------------------------

    def charged_counts(self, gene: Optional[str] = None,
                       positive: str = POSITIVE_RESIDUES,
                       negative: str = NEGATIVE_RESIDUES) -> Dict[str, int]:
        """12. Yüklü kalıntı sayıları: ``{"positive": R + K, "negative": D + E}``."""
        counts = self.cub.amino_acid_counts(gene, positive + negative)
        return {"positive": sum(counts[aa] for aa in positive),
                "negative": sum(counts[aa] for aa in negative)}

    def simple_net_charge(self, gene: Optional[str] = None,
                          positive: str = POSITIVE_RESIDUES,
                          negative: str = NEGATIVE_RESIDUES) -> int:
        """13. Kaba net yük (pH 7): (R + K) − (D + E); His ve uçlar ihmal edilir."""
        counts = self.charged_counts(gene, positive, negative)
        return counts["positive"] - counts["negative"]

    def fcr(self, gene: Optional[str] = None,
            positive: str = POSITIVE_RESIDUES,
            negative: str = NEGATIVE_RESIDUES) -> float:
        """14. FCR (yüklü kalıntı oranı): (R + K + D + E) / N."""
        n = self.__mass.length(gene)
        counts = self.charged_counts(gene, positive, negative)
        return (counts["positive"] + counts["negative"]) / n if n else 0.0

    def ncpr(self, gene: Optional[str] = None,
             positive: str = POSITIVE_RESIDUES,
             negative: str = NEGATIVE_RESIDUES) -> float:
        """15. NCPR (kalıntı başına net yük): [(R + K) − (D + E)] / N."""
        n = self.__mass.length(gene)
        return self.simple_net_charge(gene, positive, negative) / n if n else 0.0

    def net_charge(self, gene: Optional[str] = None, ph: float = 7.0,
                   pka: str | Mapping = DEFAULT_PKA_SET,
                   positive: str = IONIZABLE_POSITIVE,
                   negative: str = IONIZABLE_NEGATIVE) -> float:
        """16. pH'a bağlı net yük Q(pH), Henderson–Hasselbalch ile.

        Args:
            ph: pH değeri.
            pka: ``PKA_SETS`` içindeki kümenin adı ya da aynı anahtarlarla bir sözlük.
            positive, negative: İyonlaşan pozitif ve negatif yan zincirler.

        Raises:
            ValueError: pka geçersiz bir küme adıysa.
        """
        pk = self._pka_set(pka)
        counts = self.cub.amino_acid_counts(gene, positive + negative)

        def plus(value: float) -> float:
            return 1 / (1 + 10 ** (ph - value))

        def minus(value: float) -> float:
            return 1 / (1 + 10 ** (value - ph))

        q = sum(counts[aa] * plus(pk[aa]) for aa in positive)
        q -= sum(counts[aa] * minus(pk[aa]) for aa in negative)
        n_term, c_term = pk.get("N_term_residue", {}), pk.get("C_term_residue", {})
        for first, last in self.termini(gene):
            q += plus(n_term.get(first, pk["N_term"]))
            q -= minus(c_term.get(last, pk["C_term"]))
        return q

    def isoelectric_point(self, gene: Optional[str] = None,
                          pka: str | Mapping = DEFAULT_PKA_SET,
                          tolerance: float = PI_TOLERANCE,
                          ph_range: Tuple[float, float] = PH_RANGE) -> float:
        """17. Teorik pI: Q(pH) = 0 olan pH, ``ph_range`` içinde ikiye bölme ile.

        Q pH ile azaldığından kök tektir. Q aralığın iki ucunda aynı
        işaretliyse yakın uç döner.

        Args:
            tolerance: Aramanın durduğu pH aralığı genişliği.
        """
        low, high = ph_range
        if self.net_charge(gene, low, pka) <= 0:
            return low
        if self.net_charge(gene, high, pka) >= 0:
            return high
        while high - low > tolerance:
            middle = (low + high) / 2
            if self.net_charge(gene, middle, pka) > 0:
                low = middle
            else:
                high = middle
        return (low + high) / 2

    def titration_curve(self, gene: Optional[str] = None,
                        pka: str | Mapping = DEFAULT_PKA_SET,
                        step: float = 0.1,
                        ph_range: Tuple[float, float] = PH_RANGE) -> List[Tuple[float, float]]:
        """18. Titrasyon eğrisi: ``ph_range`` boyunca ``step`` aralıklarla (pH, Q(pH)) çiftleri.

        Raises:
            ValueError: step 0 ya da negatifse.
        """
        if step <= 0:
            raise ValueError(f"step 0'dan büyük olmalı, {step!r} verildi")
        low, high = ph_range
        count = int(round((high - low) / step))
        return [(round(low + i * step, 10), self.net_charge(gene, low + i * step, pka))
                for i in range(count + 1)]

    def __repr__(self) -> str:
        return f"ProteinCharge({self.cub!r})"
