"""Protein özellikleri (E-A … E-G): kütle, atom bileşimi, yük, hidropati, optik,
stabilite ve düzensizlik.

``ProteinProperties`` iki tür kaynakla çalışır:

- bir ``CUB``: proteinler genlerin CDS'lerinden okunan amino asitlerdir
  (``CUB``'da her gen proteinle birebir örtüşen ORF'ye indirgenir);
- bir ya da birden çok ``PRT``: dizideki 20 standart amino asit sayılır
  (``*``, ``-`` ve belirsiz kodlar atlanır).

Metotlar ``gene`` (CUB'da gen adı, PRT'de protein adı) verilmezse tüm
proteinleri birlikte, verilirse yalnızca onu hesaplar. Birlikte hesapta her
protein ayrı bir zincirdir: kütleye her biri için bir H2O eklenir, uç gruplar
ve disülfit eşleri zincir başına sayılır.

==========  ====  ==============================================================  ===========================
 Bölüm      No    Ölçüm                                                           Metot
==========  ====  ==============================================================  ===========================
 E-A        1     Uzunluk (N): Σ n_i                                              ``length``
 sayım,     2     Amino asit bileşimi (%): n_i / N × 100                          ``composition``
 kütle      3     Amino asit grubu bileşimi: Σ_grup n_i / N                       ``group_composition``
            4     Ortalama moleküler ağırlık: Σ n_i·M_i + H2O                     ``molecular_weight``
            5     Monoizotopik kütle                                              ``monoisotopic_mass``
            6     Ortalama kalıntı kütlesi: MW / N                                ``mean_residue_mass``
            7     Cys, Trp, Tyr sayısı; disülfit eşi ⌊nCys / 2⌋                   ``residue_counts``,
                                                                                  ``disulfide_pairs``
 E-B        8     Atom formülü: Σ n_i × formül − (N − 1) H2O                      ``formula``, ``formula_string``
 atomlar    9     Toplam atom sayısı                                              ``atom_count``
            10    Kalıntı başına atom sayısı                                      ``atoms_per_residue``
            11    Ortalama karbon oksidasyon durumu (Dick)                        ``carbon_oxidation_state``
 E-C        12    Pozitif (R + K) ve negatif (D + E) kalıntı sayısı               ``charged_counts``
 yük        13    Kaba net yük (pH 7): (R + K) − (D + E)                          ``simple_net_charge``
            14    FCR: (R + K + D + E) / N                                        ``fcr``
            15    NCPR: [(R + K) − (D + E)] / N                                   ``ncpr``
            16    Q(pH): Henderson–Hasselbalch                                    ``net_charge``
            17    Teorik pI: Q(pH) = 0 (ikiye bölme)                              ``isoelectric_point``
            18    Titrasyon eğrisi                                                ``titration_curve``
 E-D        19    GRAVY (Kyte–Doolittle ortalaması)                               ``gravy``
 hidropati  20    Diğer ölçeklerle ortalama                                       ``hydropathy``
            21    Aromatiklik: (F + Y + W) / N                                    ``aromaticity``
            22    Alifatik indeks (Ikai)                                          ``aliphatic_index``
            23    Boman indeksi                                                   ``boman_index``
 E-E        24    ε280: 5500·nTrp + 1490·nTyr + 125·nCystine                      ``extinction_coefficient``
 optik      25    %0,1 absorbans: ε / MW                                          ``absorbance_01``
            26    Kısmi spesifik hacim (Cohn–Edsall)                              ``partial_specific_volume``
            27    Ortalama kalıntı hacmi / yüzey alanı                            ``mean_residue_volume``,
                                                                                  ``mean_residue_surface``
 E-F        28    IVYWREL oranı                                                   ``ivywrel``
 stabilite  29    CvP yanlılığı (%)                                               ``cvp_bias``
            30    Wilkinson–Harrison çözünürlüğü                                  ``solubility_cv``,
                                                                                  ``solubility_probability``
            31    Bileşimin Shannon entropisi                                     ``composition_entropy``
 E-G        32    Uversky grafiği                                                 ``uversky``
 düzensiz-  33    TOP-IDP ortalaması                                              ``top_idp``
 lik        34    Bileşim tabanlı doğrusal skor                                   ``composition_features``,
                                                                                  ``linear_score``
==========  ====  ==============================================================  ===========================

Tablolar ``crn.tables`` içindedir ve parametreyle değiştirilebilir. pKa
kümeleri: ``"Bjellqvist"`` (varsayılan; ExPASy), ``"EMBOSS"``,
``"Lehninger"``; Bjellqvist'te uç pKa'ları uçtaki kalıntıya göre değişir. Q(pH)
hesabında uç gruplar serbest, tüm Cys'ler indirgenmiş sayılır.

**Wilkinson–Harrison** sabitleri Davis ve ark. (1999) modelinindir; CV − CV' > 0
çözünmez, < 0 çözünür tahmin edilir. **Madde 34**: yayımlanmış modellerin
(ör. Protein-Sol) ağırlıkları burada yoktur; ``linear_score``'a verilmelidir.

Kullanım::

    from crn import ProteinProperties

    props = ProteinProperties(cub)                 # ya da ProteinProperties(prt)
    props.molecular_weight(gene="MEFV_NM_000243.3")
    props.isoelectric_point(gene="MEFV_NM_000243.3", pka="EMBOSS")
    props.summary(gene="MEFV_NM_000243.3")
"""
from __future__ import annotations

import math
from collections import Counter
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from crn.constants import AMINO_ACIDS
from crn.seq import PRT
from crn.tables import (AMINO_ACID_GROUPS, RESIDUE_AVERAGE_MASS, RESIDUE_MONOISOTOPIC_MASS,
                        WATER_AVERAGE_MASS, WATER_MONOISOTOPIC_MASS, ELEMENTS,
                        AMINO_ACID_FORMULAS, WATER_FORMULA, POSITIVE_RESIDUES, NEGATIVE_RESIDUES,
                        IONIZABLE_POSITIVE, IONIZABLE_NEGATIVE, PKA_SETS, HYDROPATHY_SCALES,
                        BOMAN_SCALE, EXTINCTION_COEFFICIENTS_280, RESIDUE_SPECIFIC_VOLUMES,
                        RESIDUE_VOLUMES, RESIDUE_MAX_ASA, IVYWREL_RESIDUES, CVP_CHARGED,
                        CVP_POLAR, TURN_FORMING_RESIDUES, TOP_IDP_SCALE)

__all__ = ["ProteinProperties"]

# Varsayılan pKa kümesi ve pI aramasının sınırları, duyarlılığı.
DEFAULT_PKA_SET = "Bjellqvist"
PH_RANGE = (0.0, 14.0)
PI_TOLERANCE = 1e-4
# Ikai (1980) alifatik indeksindeki Val ve Ile/Leu katsayıları.
ALIPHATIC_V = 2.9
ALIPHATIC_IL = 3.9
# Wilkinson–Harrison modelinin sabitleri (Davis ve ark. 1999).
WH_LAMBDA_TURN = 15.43
WH_LAMBDA_CHARGE = -29.56
WH_CHARGE_OFFSET = 0.03
WH_DISCRIMINANT = 1.71
WH_PROBABILITY = (0.4934, 0.276, -0.0392)
# Uversky (2000) sınır çizgisi: ⟨H⟩ = (⟨R⟩ + UVERSKY_OFFSET) / UVERSKY_SLOPE.
UVERSKY_OFFSET = 1.151
UVERSKY_SLOPE = 2.785
# Kyte–Doolittle'ın 0-1 aralığına çekilmesi: (KD + 4,5) / 9.
_KD = HYDROPATHY_SCALES["Kyte-Doolittle"]
NORMALIZED_KYTE_DOOLITTLE = {aa: (v - min(_KD.values())) / (max(_KD.values()) - min(_KD.values()))
                             for aa, v in _KD.items()}


class _Proteins:
    """PRT listesini CUB ile aynı sayım arayüzüne (amino_acid_counts, termini, chain_count) uyarlar."""

    def __init__(self, proteins: Iterable[PRT]):
        self._counts: Dict[str, Counter] = {}
        self._termini: Dict[str, Tuple[str, str]] = {}
        for prt in proteins:
            if not isinstance(prt, PRT):
                raise TypeError(f"Protein bir PRT olmalı, {type(prt).__name__} verildi")
            if prt.name in self._counts:
                raise ValueError(f"Aynı adlı proteinler var: {prt.name}")
            residues = [aa for aa in prt.seq if aa in AMINO_ACIDS]
            self._counts[prt.name] = Counter(residues)
            if residues:
                self._termini[prt.name] = (residues[0], residues[-1])
        if not self._counts:
            raise ValueError("En az bir protein verilmeli")
        self._pooled = sum(self._counts.values(), Counter())

    @property
    def names(self) -> List[str]:
        return list(self._counts)

    def _selected(self, gene: Optional[str]) -> List[str]:
        if gene is not None and gene not in self._counts:
            raise KeyError(f"'{gene}' proteini yok")
        return self.names if gene is None else [gene]

    def amino_acid_counts(self, gene: Optional[str] = None,
                          amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, int]:
        counter = self._pooled if gene is None else self._counts[self._selected(gene)[0]]
        return {aa: counter.get(aa, 0) for aa in amino_acids}

    def termini(self, gene: Optional[str] = None) -> List[Tuple[str, str]]:
        return [self._termini[n] for n in self._selected(gene) if n in self._termini]

    def chain_count(self, gene: Optional[str] = None) -> int:
        return len(self._selected(gene))


class ProteinProperties:
    """Proteinlerin bileşimden hesaplanan fizikokimyasal özellikleri.

    Args:
        source: Bir ``CUB``, bir ``PRT`` ya da ``PRT`` listesi.

    Raises:
        TypeError: source bunlardan biri değilse.
        ValueError: PRT listesi boşsa ya da adlar tekrarlıyorsa.
    """

    def __init__(self, source):
        from crn.cub.core import CUB

        if isinstance(source, CUB):
            self.source = source
        elif isinstance(source, PRT):
            self.source = _Proteins([source])
        elif isinstance(source, (list, tuple)):
            self.source = _Proteins(source)
        else:
            raise TypeError(f"source bir CUB, PRT ya da PRT listesi olmalı, "
                            f"{type(source).__name__} verildi")

    @property
    def names(self) -> List[str]:
        return self.source.names

    def _counts(self, gene: Optional[str], amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, int]:
        return self.source.amino_acid_counts(gene, amino_acids)

    # --- E-A. Sayım ve kütle ------------------------------------------------------

    def chain_count(self, gene: Optional[str] = None) -> int:
        """Hesaba giren zincir sayısı: gene None ise protein sayısı, değilse 1."""
        return self.source.chain_count(gene)

    def length(self, gene: Optional[str] = None, amino_acids: Sequence[str] = AMINO_ACIDS) -> int:
        """1. Uzunluk (N): amino asit sayılarının toplamı (stop hariç)."""
        return sum(self._counts(gene, amino_acids).values())

    def composition(self, gene: Optional[str] = None,
                    amino_acids: Sequence[str] = AMINO_ACIDS) -> Dict[str, float]:
        """2. Amino asit bileşimi: her amino asidin yüzdesi, n_i / N × 100."""
        counts = self._counts(gene, amino_acids)
        total = sum(counts.values())
        return {aa: (n / total * 100 if total else 0.0) for aa, n in counts.items()}

    def group_composition(self, gene: Optional[str] = None,
                          groups: Mapping[str, str] = AMINO_ACID_GROUPS) -> Dict[str, float]:
        """3. Amino asit grubu bileşimi: gruptaki amino asitlerin toplam oranı (0-1)."""
        counts = self._counts(gene)
        total = sum(counts.values())
        return {name: (sum(counts.get(aa, 0) for aa in members) / total if total else 0.0)
                for name, members in groups.items()}

    def scale_average(self, scale: Mapping[str, float], gene: Optional[str] = None) -> float:
        """Kalıntı ölçeğinin bileşimle ağırlıklı ortalaması: Σ(n_i × ölçek_i) / N.

        GRAVY, Boman, kalıntı hacmi, TOP-IDP gibi ölçümlerin ortak formülü;
        yalnızca ``scale`` içindeki amino asitler sayılır.
        """
        counts = self._counts(gene, tuple(scale))
        n = sum(counts.values())
        return sum(count * scale[aa] for aa, count in counts.items()) / n if n else 0.0

    def molecular_weight(self, gene: Optional[str] = None,
                         masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS,
                         water: float = WATER_AVERAGE_MASS) -> float:
        """4. Ortalama moleküler ağırlık (Da): Σ(n_i × M_i) + zincir başına H2O."""
        counts = self._counts(gene, tuple(masses))
        if not sum(counts.values()):
            return 0.0
        return sum(n * masses[aa] for aa, n in counts.items()) + water * self.chain_count(gene)

    def monoisotopic_mass(self, gene: Optional[str] = None,
                          masses: Mapping[str, float] = RESIDUE_MONOISOTOPIC_MASS,
                          water: float = WATER_MONOISOTOPIC_MASS) -> float:
        """5. Monoizotopik kütle (Da)."""
        return self.molecular_weight(gene, masses, water)

    def mean_residue_mass(self, gene: Optional[str] = None,
                          masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS,
                          water: float = WATER_AVERAGE_MASS) -> float:
        """6. Ortalama kalıntı kütlesi (Da): MW / N."""
        n = self.length(gene, tuple(masses))
        return self.molecular_weight(gene, masses, water) / n if n else 0.0

    def residue_counts(self, gene: Optional[str] = None, residues: Sequence[str] = "CWY") -> Dict[str, int]:
        """7. Seçilen kalıntıların sayısı; varsayılan Cys, Trp, Tyr."""
        return self._counts(gene, residues)

    def disulfide_pairs(self, gene: Optional[str] = None) -> int:
        """7. Olası disülfit bağı sayısı: zincir başına ⌊nCys / 2⌋."""
        if gene is None:
            return sum(self.disulfide_pairs(name) for name in self.names)
        return self.residue_counts(gene, "C")["C"] // 2

    # --- E-B. Atom bileşimi -------------------------------------------------------

    def formula(self, gene: Optional[str] = None,
                formulas: Mapping[str, Tuple[int, ...]] = AMINO_ACID_FORMULAS,
                water: Tuple[int, ...] = WATER_FORMULA,
                elements: Sequence[str] = ELEMENTS) -> Dict[str, int]:
        """8. Atom formülü: Σ(n_i × amino asit formülü) − (N − zincir sayısı) × H2O."""
        counts = self._counts(gene, tuple(formulas))
        n = sum(counts.values())
        bonds = n - self.chain_count(gene) if n else 0
        return {element: sum(count * formulas[aa][i] for aa, count in counts.items()) - bonds * water[i]
                for i, element in enumerate(elements)}

    def formula_string(self, gene: Optional[str] = None, **kwargs) -> str:
        """8. Atom formülü metin olarak, ör. ``"C3779H6038N1104O1167S26"`` (0 olanlar yazılmaz)."""
        return "".join(f"{element}{n}" for element, n in self.formula(gene, **kwargs).items() if n)

    def atom_count(self, gene: Optional[str] = None, **kwargs) -> int:
        """9. Toplam atom sayısı."""
        return sum(self.formula(gene, **kwargs).values())

    def atoms_per_residue(self, gene: Optional[str] = None, **kwargs) -> float:
        """10. Kalıntı başına atom sayısı: toplam atom / N."""
        n = self.length(gene, tuple(kwargs.get("formulas", AMINO_ACID_FORMULAS)))
        return self.atom_count(gene, **kwargs) / n if n else 0.0

    def carbon_oxidation_state(self, gene: Optional[str] = None, charge: float = 0.0, **kwargs) -> float:
        """11. Ortalama karbon oksidasyon durumu (Dick): (Z − nH + 3nN + 2nO + 2nS) / nC.

        Args:
            charge: Proteinin net yükü Z; nötr protein için 0.
        """
        f = self.formula(gene, **kwargs)
        if not f.get("C"):
            return 0.0
        return (charge - f.get("H", 0) + 3 * f.get("N", 0) + 2 * f.get("O", 0)
                + 2 * f.get("S", 0)) / f["C"]

    # --- E-C. Yük ve iyonlaşma -----------------------------------------------------

    def termini(self, gene: Optional[str] = None) -> List[Tuple[str, str]]:
        """Her zincirin (ilk, son) kalıntısı."""
        return self.source.termini(gene)

    def charged_counts(self, gene: Optional[str] = None, positive: str = POSITIVE_RESIDUES,
                       negative: str = NEGATIVE_RESIDUES) -> Dict[str, int]:
        """12. Yüklü kalıntı sayıları: ``{"positive": R + K, "negative": D + E}``."""
        counts = self._counts(gene, positive + negative)
        return {"positive": sum(counts[aa] for aa in positive),
                "negative": sum(counts[aa] for aa in negative)}

    def simple_net_charge(self, gene: Optional[str] = None, positive: str = POSITIVE_RESIDUES,
                          negative: str = NEGATIVE_RESIDUES) -> int:
        """13. Kaba net yük (pH 7): (R + K) − (D + E); His ve uçlar ihmal edilir."""
        counts = self.charged_counts(gene, positive, negative)
        return counts["positive"] - counts["negative"]

    def fcr(self, gene: Optional[str] = None, positive: str = POSITIVE_RESIDUES,
            negative: str = NEGATIVE_RESIDUES) -> float:
        """14. FCR (yüklü kalıntı oranı): (R + K + D + E) / N."""
        n = self.length(gene)
        counts = self.charged_counts(gene, positive, negative)
        return (counts["positive"] + counts["negative"]) / n if n else 0.0

    def ncpr(self, gene: Optional[str] = None, positive: str = POSITIVE_RESIDUES,
             negative: str = NEGATIVE_RESIDUES) -> float:
        """15. NCPR (kalıntı başına net yük): [(R + K) − (D + E)] / N."""
        n = self.length(gene)
        return self.simple_net_charge(gene, positive, negative) / n if n else 0.0

    @staticmethod
    def _pka_set(pka: str | Mapping) -> Mapping:
        if isinstance(pka, str):
            if pka not in PKA_SETS:
                raise ValueError(f"pka {', '.join(map(repr, PKA_SETS))} ya da bir sözlük "
                                 f"olmalı, {pka!r} verildi")
            return PKA_SETS[pka]
        return pka

    def net_charge(self, gene: Optional[str] = None, ph: float = 7.0,
                   pka: str | Mapping = DEFAULT_PKA_SET,
                   positive: str = IONIZABLE_POSITIVE,
                   negative: str = IONIZABLE_NEGATIVE) -> float:
        """16. pH'a bağlı net yük Q(pH), Henderson–Hasselbalch ile.

        Args:
            pka: ``PKA_SETS`` içindeki kümenin adı ya da aynı anahtarlarla sözlük.

        Raises:
            ValueError: pka geçersiz bir küme adıysa.
        """
        pk = self._pka_set(pka)
        counts = self._counts(gene, positive + negative)

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

        Q aralığın iki ucunda aynı işaretliyse yakın uç döner.
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
        """18. Titrasyon eğrisi: ``step`` aralıklarla (pH, Q(pH)) çiftleri.

        Raises:
            ValueError: step 0 ya da negatifse.
        """
        if step <= 0:
            raise ValueError(f"step 0'dan büyük olmalı, {step!r} verildi")
        low, high = ph_range
        count = int(round((high - low) / step))
        return [(round(low + i * step, 10), self.net_charge(gene, low + i * step, pka))
                for i in range(count + 1)]

    # --- E-D. Hidrofobiklik ve polarite ----------------------------------------------

    def gravy(self, gene: Optional[str] = None) -> float:
        """19. GRAVY: Kyte–Doolittle değerlerinin ortalaması."""
        return self.hydropathy(gene, "Kyte-Doolittle")

    def hydropathy(self, gene: Optional[str] = None,
                   scale: str | Mapping[str, float] = "Kyte-Doolittle") -> float:
        """20. Seçilen ölçeğin ortalaması (``HYDROPATHY_SCALES`` adı ya da sözlük).

        Raises:
            ValueError: scale geçersiz bir ölçek adıysa.
        """
        if isinstance(scale, str):
            if scale not in HYDROPATHY_SCALES:
                raise ValueError(f"scale {', '.join(map(repr, HYDROPATHY_SCALES))} ya da bir "
                                 f"sözlük olmalı, {scale!r} verildi")
            scale = HYDROPATHY_SCALES[scale]
        return self.scale_average(scale, gene)

    def aromaticity(self, gene: Optional[str] = None,
                    aromatic: str = AMINO_ACID_GROUPS["aromatic"]) -> float:
        """21. Aromatiklik: (F + Y + W) / N."""
        return self.group_composition(gene, {"aromatic": aromatic})["aromatic"]

    def aliphatic_index(self, gene: Optional[str] = None,
                        v: float = ALIPHATIC_V, il: float = ALIPHATIC_IL) -> float:
        """22. Alifatik indeks (Ikai): X(A) + a·X(V) + b·[X(I) + X(L)], X mol yüzdesi."""
        x = self.composition(gene)
        return x["A"] + v * x["V"] + il * (x["I"] + x["L"])

    def boman_index(self, gene: Optional[str] = None,
                    scale: Mapping[str, float] = BOMAN_SCALE) -> float:
        """23. Boman indeksi (kcal/mol): Σ(n_i × ΔG_i) / N; ~2,48 üstü yüksek bağlanma potansiyeli."""
        return self.scale_average(scale, gene)

    # --- E-E. Optik ve hidrodinamik özellikler -----------------------------------------

    def extinction_coefficient(self, gene: Optional[str] = None, reduced: bool = False,
                               coefficients: Mapping[str, float] = EXTINCTION_COEFFICIENTS_280) -> float:
        """24. 280 nm molar ekstinksiyon katsayısı (M⁻¹ cm⁻¹).

        Args:
            reduced: True ise tüm Cys'ler serbest (cystine terimi 0); False ise
                tümü eşleşmiş varsayılır. ExPASy ProtParam iki değeri de verir.
        """
        counts = self.residue_counts(gene, "WY")
        cystines = 0 if reduced else self.disulfide_pairs(gene)
        return (coefficients["W"] * counts["W"] + coefficients["Y"] * counts["Y"]
                + coefficients["cystine"] * cystines)

    def absorbance_01(self, gene: Optional[str] = None, reduced: bool = False,
                      coefficients: Mapping[str, float] = EXTINCTION_COEFFICIENTS_280) -> float:
        """25. %0,1 (1 mg/mL) çözeltinin 280 nm absorbansı: ε / MW."""
        mw = self.molecular_weight(gene)
        return self.extinction_coefficient(gene, reduced, coefficients) / mw if mw else 0.0

    def partial_specific_volume(self, gene: Optional[str] = None,
                                volumes: Mapping[str, float] = RESIDUE_SPECIFIC_VOLUMES,
                                masses: Mapping[str, float] = RESIDUE_AVERAGE_MASS) -> float:
        """26. Kısmi spesifik hacim v̄ (cm³/g), Cohn–Edsall: Σ(n_i·M_i·v_i) / Σ(n_i·M_i)."""
        weights = {aa: masses[aa] for aa in volumes}
        weighted = {aa: masses[aa] * volumes[aa] for aa in volumes}
        total = self.scale_average(weights, gene)
        return self.scale_average(weighted, gene) / total if total else 0.0

    def mean_residue_volume(self, gene: Optional[str] = None,
                            volumes: Mapping[str, float] = RESIDUE_VOLUMES) -> float:
        """27. Ortalama kalıntı hacmi (Å³, Zamyatnin)."""
        return self.scale_average(volumes, gene)

    def mean_residue_surface(self, gene: Optional[str] = None,
                             surfaces: Mapping[str, float] = RESIDUE_MAX_ASA) -> float:
        """27. Ortalama kalıntı yüzey alanı (Å², en büyük ASA, Tien ve ark. 2013)."""
        return self.scale_average(surfaces, gene)

    # --- E-F. Stabilite ve çözünürlük ------------------------------------------------

    def ivywrel(self, gene: Optional[str] = None, residues: str = IVYWREL_RESIDUES) -> float:
        """28. IVYWREL oranı: (I + V + Y + W + R + E + L) / N."""
        return self.group_composition(gene, {"ivywrel": residues})["ivywrel"]

    def cvp_bias(self, gene: Optional[str] = None,
                 charged: str = CVP_CHARGED, polar: str = CVP_POLAR) -> float:
        """29. CvP yanlılığı (%): [(D + E + K + R) − (N + Q + S + T)] / N × 100."""
        groups = self.group_composition(gene, {"charged": charged, "polar": polar})
        return (groups["charged"] - groups["polar"]) * 100

    def solubility_cv(self, gene: Optional[str] = None,
                      turn: str = TURN_FORMING_RESIDUES,
                      lambda_turn: float = WH_LAMBDA_TURN,
                      lambda_charge: float = WH_LAMBDA_CHARGE,
                      offset: float = WH_CHARGE_OFFSET) -> float:
        """30. Wilkinson–Harrison kanonik değişkeni: λ1·(N + G + P + S)/N + λ2·|NCPR − 0,03|."""
        fraction = self.group_composition(gene, {"turn": turn})["turn"]
        return lambda_turn * fraction + lambda_charge * abs(self.ncpr(gene) - offset)

    def solubility_probability(self, gene: Optional[str] = None,
                               discriminant: float = WH_DISCRIMINANT) -> Dict[str, float | bool]:
        """30. Wilkinson–Harrison tahmini.

        Returns:
            ``{"cv": CV, "soluble": CV − CV' <= 0, "probability": tahmin edilen
            sınıfın olasılığı}``.
        """
        cv = self.solubility_cv(gene)
        d = cv - discriminant
        a, b, c = WH_PROBABILITY
        return {"cv": cv, "soluble": d <= 0, "probability": a + b * abs(d) + c * d * d}

    def composition_entropy(self, gene: Optional[str] = None) -> float:
        """31. Bileşimin Shannon entropisi (bit): −Σ p_i log2 p_i; en çok log2 20 ≈ 4,32."""
        return -sum(p / 100 * math.log2(p / 100) for p in self.composition(gene).values() if p > 0)

    # --- E-G. Düzensizlik -----------------------------------------------------------

    def uversky(self, gene: Optional[str] = None,
                hydropathy: Mapping[str, float] = NORMALIZED_KYTE_DOOLITTLE) -> Dict[str, float | bool]:
        """32. Uversky grafiğindeki konum.

        Returns:
            ``{"mean_net_charge": ⟨R⟩ = |NCPR|, "mean_hydropathy": ⟨H⟩,
            "boundary": sınırdaki ⟨H⟩, "disordered": ⟨H⟩ < sınır}``.
        """
        r = abs(self.ncpr(gene))
        h = self.scale_average(hydropathy, gene)
        boundary = (r + UVERSKY_OFFSET) / UVERSKY_SLOPE
        return {"mean_net_charge": r, "mean_hydropathy": h,
                "boundary": boundary, "disordered": h < boundary}

    def top_idp(self, gene: Optional[str] = None, scale: Mapping[str, float] = TOP_IDP_SCALE) -> float:
        """33. TOP-IDP ortalaması; pozitif değer düzensizlik eğilimi."""
        return self.scale_average(scale, gene)

    def composition_features(self, gene: Optional[str] = None) -> Dict[str, float]:
        """34. Bileşim tabanlı çözünürlük/toplanma modellerinin kullandığı özellikler.

        Anahtarlar: amino asit kodları (oran, 0-1), ``"K-R"``, ``"D-E"``,
        ``"K+R"``, ``"D+E"``, ``"K+R-D-E"``, ``"K+R+D+E"``, ``"F+W+Y"``, ``"pI"``,
        ``"abs_charge"`` (|Q(pH 7)| / N), ``"gravy"``, ``"entropy"``,
        ``"top_idp"``, ``"fcr"``, ``"ncpr"``.
        """
        x = {aa: p / 100 for aa, p in self.composition(gene).items()}
        n = self.length(gene)
        result = dict(x)
        result.update({
            "K-R": x["K"] - x["R"],
            "D-E": x["D"] - x["E"],
            "K+R": x["K"] + x["R"],
            "D+E": x["D"] + x["E"],
            "K+R-D-E": x["K"] + x["R"] - x["D"] - x["E"],
            "K+R+D+E": x["K"] + x["R"] + x["D"] + x["E"],
            "F+W+Y": x["F"] + x["W"] + x["Y"],
            "pI": self.isoelectric_point(gene),
            "abs_charge": abs(self.net_charge(gene, 7.0)) / n if n else 0.0,
            "gravy": self.gravy(gene),
            "entropy": self.composition_entropy(gene),
            "top_idp": self.top_idp(gene),
            "fcr": self.fcr(gene),
            "ncpr": self.ncpr(gene),
        })
        return result

    def linear_score(self, gene: Optional[str] = None,
                     weights: Optional[Mapping[str, float]] = None,
                     intercept: float = 0.0) -> float:
        """34. Doğrusal bileşim skoru: sabit + Σ ağırlık × özellik.

        Raises:
            ValueError: weights verilmezse.
            KeyError: weights içinde ``composition_features``'ta olmayan anahtar varsa.
        """
        if not weights:
            raise ValueError("weights verilmeli: composition_features() anahtarı -> regresyon katsayısı")
        values = self.composition_features(gene)
        unknown = [key for key in weights if key not in values]
        if unknown:
            raise KeyError(f"Bilinmeyen özellikler: {', '.join(unknown)}; "
                           f"seçenekler: {', '.join(values)}")
        return intercept + sum(w * values[key] for key, w in weights.items())

    # --- Özet -------------------------------------------------------------------------

    def summary(self, gene: Optional[str] = None) -> Dict[str, float]:
        """Başlıca skaler özellikler tek sözlükte (ProtParam benzeri)."""
        return {
            "length": self.length(gene),
            "MW": self.molecular_weight(gene),
            "pI": self.isoelectric_point(gene),
            "net_charge_7": self.net_charge(gene, 7.0),
            "GRAVY": self.gravy(gene),
            "aromaticity": self.aromaticity(gene),
            "aliphatic_index": self.aliphatic_index(gene),
            "extinction_280": self.extinction_coefficient(gene),
            "IVYWREL": self.ivywrel(gene),
            "CvP": self.cvp_bias(gene),
            "TOP_IDP": self.top_idp(gene),
        }

    def __repr__(self) -> str:
        return f"ProteinProperties({len(self.names)} protein)"
