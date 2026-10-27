"""Kodon kullanım yanlılığı (codon usage bias) — E. Protein özellikleri, G. Bileşimden düzensizlik tahmini.

Bir ``CUB`` nesnesindeki proteinlerin bileşiminden yapısal düzensizlik
eğilimini hesaplar. Ölçek ortalamaları E-A'daki ``ProteinMass``'tan, NCPR ve
pI E-C'deki ``ProteinCharge``'tan, GRAVY E-D'deki ``ProteinHydropathy``'den,
entropi E-F'deki ``ProteinStability``'den gelir.

E-G. Bileşimden düzensizlik tahmini:

====  ==============================================================  ======================
 No   Ölçüm                                                           Metot
====  ==============================================================  ======================
 32   Uversky grafiği: ⟨R⟩ = |NCPR|, ⟨H⟩ = normalize KD ortalaması;   ``uversky``
      sınır ⟨H⟩ = (⟨R⟩ + 1,151) / 2,785, altı düzensiz
 33   TOP-IDP ortalaması: Σ(n_i × TOP-IDP_i) / N                      ``top_idp``
 34   Bileşim tabanlı doğrusal skor: Σ ağırlık × özellik + sabit      ``features``,
                                                                      ``linear_score``
====  ==============================================================  ======================

Uversky'de Kyte–Doolittle değerleri 0-1 aralığına çekilir:
(KD + 4,5) / 9. TOP-IDP ölçeği ``crn.io.crn_io_const.TOP_IDP_SCALE``'dir
(Campen ve ark. 2008); pozitif ortalama düzensizlik, negatif ortalama düzenli
yapı eğilimi demektir.

**Madde 34.** Protein-Sol gibi bileşim tabanlı çözünürlük ve toplanma
skorları, bileşimden türetilen özelliklerin regresyonla bulunmuş
ağırlıklarla toplamıdır. ``features`` bu tür modellerin kullandığı
özellikleri (amino asit oranları, yük kombinasyonları, pI, GRAVY, entropi,
TOP-IDP ...) hesaplar; ``linear_score`` verilen ağırlıklarla skoru bulur.
Yayımlanmış bir modelin ağırlıkları bu modülde yoktur; özgün makalesinden
alınıp verilmelidir.

Kullanım::

    from crn.cub import CUB, ProteinDisorder

    dis = ProteinDisorder(cub)
    dis.uversky(gene="MEFV_NM_000243.3")
    dis.top_idp(gene="MEFV_NM_000243.3")
    dis.linear_score(gene="MEFV_NM_000243.3", weights={"K-R": 0.5, "gravy": -1.2})
"""
from __future__ import annotations

from typing import (Dict, Mapping, Optional)

from crn.io.crn_io_const import (HYDROPATHY_SCALES, TOP_IDP_SCALE)
from crn.cub.crn_cub_counts import CUB
from crn.cub.crn_cub_protein_mass import ProteinMass
from crn.cub.crn_cub_protein_charge import ProteinCharge
from crn.cub.crn_cub_protein_hydropathy import ProteinHydropathy
from crn.cub.crn_cub_protein_stability import ProteinStability


# Uversky (2000) sınır çizgisi: ⟨H⟩ = (⟨R⟩ + UVERSKY_OFFSET) / UVERSKY_SLOPE.
UVERSKY_OFFSET = 1.151
UVERSKY_SLOPE = 2.785
# Kyte–Doolittle'ın 0-1 aralığına çekilmesi: (KD − en küçük) / (en büyük − en küçük).
_KD = HYDROPATHY_SCALES["Kyte-Doolittle"]
NORMALIZED_KYTE_DOOLITTLE = {aa: (v - min(_KD.values())) / (max(_KD.values()) - min(_KD.values()))
                             for aa, v in _KD.items()}


class ProteinDisorder:
    """Bir ``CUB`` nesnesindeki proteinlerin bileşimden düzensizlik eğilimi.

    Metotlar ``CUB`` ile aynı biçimde ``gene`` verilmezse tüm genleri
    birlikte, verilirse yalnızca o geni hesaplar.

    Args:
        cub: Kodonları hazırlanmış ``CUB`` nesnesi.

    Raises:
        TypeError: cub bir CUB değilse.
    """

    def __init__(self, cub: CUB):
        self.__mass = ProteinMass(cub)
        self.__charge = ProteinCharge(cub)
        self.__hydropathy = ProteinHydropathy(cub)
        self.__stability = ProteinStability(cub)

    @property
    def cub(self) -> CUB:
        return self.__mass.cub

    # --- E-G. Düzensizlik -----------------------------------------------------------

    def uversky(self, gene: Optional[str] = None,
                hydropathy: Mapping[str, float] = NORMALIZED_KYTE_DOOLITTLE) -> Dict[str, float | bool]:
        """32. Uversky grafiğindeki konum.

        Returns:
            ``{"mean_net_charge": ⟨R⟩, "mean_hydropathy": ⟨H⟩,
            "boundary": sınırdaki ⟨H⟩, "disordered": ⟨H⟩ < sınır}``.
        """
        r = abs(self.__charge.ncpr(gene))
        h = self.__mass.scale_average(hydropathy, gene)
        boundary = (r + UVERSKY_OFFSET) / UVERSKY_SLOPE
        return {"mean_net_charge": r, "mean_hydropathy": h,
                "boundary": boundary, "disordered": h < boundary}

    def top_idp(self, gene: Optional[str] = None,
                scale: Mapping[str, float] = TOP_IDP_SCALE) -> float:
        """33. TOP-IDP ortalaması: Σ(n_i × TOP-IDP_i) / N; pozitif değer düzensizlik eğilimi."""
        return self.__mass.scale_average(scale, gene)

    def features(self, gene: Optional[str] = None) -> Dict[str, float]:
        """34. Bileşim tabanlı modellerde kullanılan özellikler.

        Anahtarlar: tek harfli amino asit kodları (oran, 0-1), ``"K-R"``,
        ``"D-E"``, ``"K+R"``, ``"D+E"``, ``"K+R-D-E"``, ``"K+R+D+E"``,
        ``"F+W+Y"`` (oranlar), ``"pI"``, ``"abs_charge"`` (|Q(pH 7)| / N),
        ``"gravy"``, ``"entropy"``, ``"top_idp"``, ``"fcr"``, ``"ncpr"``.
        """
        x = {aa: p / 100 for aa, p in self.__mass.composition(gene).items()}
        n = self.__mass.length(gene)
        result = dict(x)
        result.update({
            "K-R": x["K"] - x["R"],
            "D-E": x["D"] - x["E"],
            "K+R": x["K"] + x["R"],
            "D+E": x["D"] + x["E"],
            "K+R-D-E": x["K"] + x["R"] - x["D"] - x["E"],
            "K+R+D+E": x["K"] + x["R"] + x["D"] + x["E"],
            "F+W+Y": x["F"] + x["W"] + x["Y"],
            "pI": self.__charge.isoelectric_point(gene),
            "abs_charge": abs(self.__charge.net_charge(gene, 7.0)) / n if n else 0.0,
            "gravy": self.__hydropathy.gravy(gene),
            "entropy": self.__stability.composition_entropy(gene),
            "top_idp": self.top_idp(gene),
            "fcr": self.__charge.fcr(gene),
            "ncpr": self.__charge.ncpr(gene),
        })
        return result

    def linear_score(self, gene: Optional[str] = None,
                     weights: Optional[Mapping[str, float]] = None,
                     intercept: float = 0.0) -> float:
        """34. Doğrusal bileşim skoru: sabit + Σ ağırlık × özellik.

        Args:
            weights: ``features`` anahtarı -> ağırlık (regresyon katsayısı).
            intercept: Modelin sabit terimi.

        Raises:
            ValueError: weights verilmezse.
            KeyError: weights içinde ``features``'ta olmayan bir anahtar varsa.
        """
        if not weights:
            raise ValueError("weights verilmeli: features() anahtarı -> regresyon katsayısı")
        values = self.features(gene)
        unknown = [key for key in weights if key not in values]
        if unknown:
            raise KeyError(f"Bilinmeyen özellikler: {', '.join(unknown)}; "
                           f"seçenekler: {', '.join(values)}")
        return intercept + sum(w * values[key] for key, w in weights.items())

    def __repr__(self) -> str:
        return f"ProteinDisorder({self.cub!r})"
