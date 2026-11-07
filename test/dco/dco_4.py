"""Örnek 4: Yetişkin başlangıçlı nörodejenerasyon.

Klinik arka plan: Alzheimer, Parkinson, ALS, frontotemporal demans ve Huntington, farklı
klinik tablolara rağmen protein katlanma bozukluğu, kümeleşme ve proteostaz yükünü
paylaşır. Sinonim kodonlar çeviri hızını değiştirerek eş-çevirisel katlanmayı
etkileyebildiği için bu hastalıklarda kodon kullanımı son yıllarda ilgi çekmektedir.

Soru: Bu hastalıkların bilinen genleri ortak bir kodon profili taşıyor mu? Örneğin
katlanması zor proteinlerde yavaş kodon ya da düşük tRNA uyumu var mı?

Bu örnekte:
  1. Hastalık adıyla beş hastalığın genleri tek grupta toplanır.
  2. Yedi adımlı analiz yapılır.
  3. Genler hastalığa göre ayrılıp tAI ve GC3s medyanları hastalık hastalık yazılır.

Çıktılar: test/dco/cikti/dco_4/
Not: Grup Mendel genleri ve HPO'nun çok genli etiketlediği risk genlerinden oluşur (LRRK2, GBA1).
Sporadik vakaların büyük çoğunluğu bu genlerle açıklanmaz.
"""
import re
from statistics import median

from dco_ortak import ILISKI, analiz, gen_adi, gen_al, hazirla

HASTALIKLAR = {"Alzheimer": "alzheimer", "Parkinson": "parkinson",
               "ALS": "amyotrophic lateral", "FTD": "frontotemporal", "Huntington": "huntington"}
k = hazirla()
filtre = {"association": ILISKI}
hasta = gen_al(k, disease=re.compile("|".join(HASTALIKLAR.values()), re.I), **filtre)
sonuc = analiz(k, "ND", hasta, "dco_4")

ad_tablo = {gen_adi(g): s for g, s in zip(sonuc.dcu.a.genes, sonuc.tablo_a)}
print(f"{'Hastalık':<13}{'Gen':>5}{'tAI':>8}{'GC3s':>8}{'ENC':>8}  Genlerden örnekler")
for ad, desen in HASTALIKLAR.items():
    semboller = [s for s in k.hpo.symbols(disease=re.compile(desen, re.I), **filtre) if s in ad_tablo]
    satirlar = [ad_tablo[s] for s in semboller]
    if satirlar:
        print(f"{ad:<13}{len(satirlar):>5}{median(s['tAI'] for s in satirlar):>8.3f}"
              f"{median(s['GC3s'] for s in satirlar):>8.3f}{median(s['ENC'] for s in satirlar):>8.1f}"
              f"  {', '.join(semboller[:5])}")
