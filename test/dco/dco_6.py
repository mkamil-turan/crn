"""Örnek 6: Bütün hastalıkların karşılaştırması ve yöntemin kalibrasyonu.

Önceki örnekler hastalıkları tek tek inceledi. Burada iki büyük soru sorulur:

Soru 1: Fisher testine ne kadar güvenebiliriz?
    Kontrol havuzu, biyolojik olarak anlamlı bir ayrım içermeyen rastgele iki gruba
    bölünür. Doğru çalışan bir test burada neredeyse hiç kodon bulmamalıdır, çünkü
    gruplar arasında gerçek fark yoktur. Fisher'ın burada bulduğu her şey yanlış
    pozitiftir. Aynı bölünmeler gen düzeyinde permütasyona da verilir.

Soru 2: Hangi hastalık kontrolden gerçekten ayrılıyor?
    Her hastalık aynı akıştan geçirilir: eşleştirilmiş kontrol, Fisher, permütasyon, CUB
    parametreleri. Sonuçlar tek tablo ve ısı haritasında özetlenir. Fisher sayıları
    yalnızca karşılaştırma içindir; bir bulgu ancak permütasyondan ya da gen düzeyi
    parametre testinden geçerse sayılır.

Isı haritası nasıl okunur: Her hücre Cliff deltadır. Kırmızı, parametrenin hastalık
genlerinde yüksek olduğunu, mavi düşük olduğunu gösterir. Yıldız q<0.05 demektir.

Çıktılar: test/dco/cikti/dco_6/
Not: Hastalık başına ayrı analiz klasörleri ve ayrıntılı tablolar da üretilir.
"""
import contextlib
import io
import random
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

from crn.cub import DifferentialCodonUsage
from dco_ortak import CIKTI, ILISKI, PARAMETRELER, TOHUM, Hpo, analiz, gen_al, hazirla, paketle

KLASOR = CIKTI / "dco_6"
KLASOR.mkdir(parents=True, exist_ok=True)
BOLUNME = 12

k = hazirla()

# --- Soru 1: kalibrasyon ------------------------------------------------------------
print(f"KALİBRASYON: kontrol havuzu {BOLUNME} kez rastgele ikiye bölünüyor (grup başına 70 gen)")
print("Gerçek fark yoktur, yani bulunan her anlamlı kodon yanlış pozitiftir.\n")
fisher_sayilari, perm_sayilari = [], []
for tohum in range(BOLUNME):
    secim = random.Random(tohum).sample(k.havuz, 140)
    paket = KLASOR / "kalibrasyon.crn"
    paketle({"a": secim[:70], "b": secim[70:]}, paket)
    dcu = DifferentialCodonUsage.from_crn(paket, "a", "b")
    fisher_sayilari.append(sum(s["significant"] for s in dcu.table("family")))
    perm_sayilari.append(sum(s["perm_q"] < .05 for s in dcu.permutation("family", n=1000, seed=TOHUM)))
paket.unlink()
print(f"Fisher, aile düzeyi anlamlı kodon: ortalama {np.mean(fisher_sayilari):.1f}/59, "
      f"en çok {max(fisher_sayilari)}/59, bölünmeler: {fisher_sayilari}")
print(f"Permütasyon:                      ortalama {np.mean(perm_sayilari):.1f}/59, "
      f"en çok {max(perm_sayilari)}/59\n")

# --- Soru 2: hastalıklar -------------------------------------------------------------
GRUPLAR = {
    "SCA": dict(disease=re.compile("spinocerebellar ataxia", re.I), association=ILISKI),
    "CMT": dict(disease=re.compile("charcot-marie-tooth", re.I), association=ILISKI),
    "RP": dict(disease=re.compile("retinitis pigmentosa", re.I), association=ILISKI),
    "ND": dict(disease=re.compile("alzheimer|parkinson|amyotrophic lateral|frontotemporal|huntington", re.I),
               association=ILISKI),
    "Otizm": dict(term="Autistic behavior", association=ILISKI, evidence=["PCS", "TAS"], min_frequency=0.3),
    "Sizofreni": dict(term="Schizophrenia"),
}

sonuclar = {}
for ad, filtre in GRUPLAR.items():
    genler = gen_al(k, **filtre)
    with contextlib.redirect_stdout(io.StringIO()):     # ayrıntılar kendi klasörlerinde
        sonuclar[ad] = analiz(k, ad, genler, f"dco_6_{ad}")

print(f"{'Hastalık':<11}{'Gen':>5}{'Hastalık':>10}{'Fisher':>8}{'Sağlam':>8}{'Anlamlı param':>15}  En büyük etki")
tablo = []
for ad, s in sonuclar.items():
    anlamli = [p["parametre"] for p in s.param if p["q"] < .05]
    enb = max(s.param, key=lambda p: abs(p["delta"]))
    print(f"{ad:<11}{len(s.hasta):>5}{s.klinik['hastalik']:>10}{s.fisher_aile:>8}{len(s.saglam):>8}"
          f"{len(anlamli):>15}  {enb['parametre']} delta {enb['delta']:+.2f} (q {enb['q']:.2f})")
    tablo.append({"hastalik": ad, "gen": len(s.hasta), "hastalik_kaydi": s.klinik["hastalik"],
                  "fisher_aile": s.fisher_aile, "permutasyon_saglam": len(s.saglam),
                  "saglam_kodonlar": "|".join(x["codon"] for x in s.saglam),
                  "anlamli_parametreler": "|".join(anlamli), "en_buyuk_etki": enb["parametre"],
                  "delta": round(enb["delta"], 3), "q": round(enb["q"], 4)})
Hpo.write_tsv(tablo, KLASOR / "ozet.tsv")

# --- Isı haritası --------------------------------------------------------------------
# Renk polariteyi taşır: mavi düşük, kırmızı yüksek, gri orta nokta "fark yok" demektir.
# Sayılar her hücrede yazılıdır, yani renk görülemese de tablo okunabilir.
adlar, parametreler = list(sonuclar), list(PARAMETRELER)
degerler = np.array([[next(p["delta"] for p in sonuclar[a].param if p["parametre"] == q) for a in adlar]
                     for q in parametreler])
qler = np.array([[next(p["q"] for p in sonuclar[a].param if p["parametre"] == q) for a in adlar]
                 for q in parametreler])
renk = LinearSegmentedColormap.from_list("mavi_kirmizi", ["#2a78d6", "#f0efec", "#e34948"])
fig, ax = plt.subplots(figsize=(8.5, 6.5))
ax.imshow(degerler, cmap=renk, vmin=-.6, vmax=.6, aspect="auto")
ax.set_xticks(range(len(adlar)), [f"{a}\n{len(sonuclar[a].hasta)} gen" for a in adlar])
ax.set_yticks(range(len(parametreler)), parametreler)
ax.tick_params(length=0)
ax.xaxis.tick_top()
for kenar in ax.spines.values():
    kenar.set_visible(False)
for i in range(len(parametreler)):
    for j in range(len(adlar)):
        ax.text(j, i, f"{degerler[i, j]:+.2f}{' *' if qler[i, j] < .05 else ''}", ha="center", va="center",
                fontsize=10, color="#1a1a19")
ax.set_title("Hastalık genleri ve eşleştirilmiş kontroller: Cliff delta", pad=44, fontsize=12, loc="left")
fig.text(.01, .01, "Kırmızı: hastalıkta yüksek. Mavi: hastalıkta düşük. *: q<0.05. L ve GC eşleştirildiği için ~0 beklenir.",
         fontsize=9, color="#555555")
fig.tight_layout(rect=(0, .03, 1, 1))
fig.savefig(KLASOR / "etki_isi_haritasi.png", dpi=150)
print(f"\nÇıktılar: {KLASOR}")
