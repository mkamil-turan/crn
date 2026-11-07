"""Örnek 5: Otizm spektrum bozukluğu ve şizofreni.

Klinik arka plan: Otizm erken çocuklukta, şizofreni genellikle ergenlik sonrasında
başlayan, genetik yatkınlığı yüksek nörogelişimsel bozukluklardır. Ortak sinaps ve
kromatin genleri paylaştıkları bilinir. Otizm genlerinin belirgin biçimde uzun olduğu da
bildirilmiştir.

Soru: İki bozukluğun genleri kodon kullanımında ayrılıyor mu? Ayrıca bir tuzağı
gösteriyoruz: eşleştirme yapılmazsa uzunluk ve GC, biyolojik bulgu gibi görünür.
Doğrudan karşılaştırmada bütün parametreler anlamlı çıkar, ama her grup kendi eşleşmiş
kontrolüyle kıyaslandığında bu farklar kaybolur.

Bu örnekte:
  1. Otizm grubu: sıklığı en az %30 ve yayın/uzman kanıtlı genler.
  2. Şizofreni grubu: HPO terimi, Mendel filtresi uygulanmadan (yoksa 15 gen kalır).
  3. Eşleştirme öncesi ve sonrası uzunluk ile GC3s karşılaştırılır.
  4. Her grup kontrolüyle, sonra iki grup doğrudan karşılaştırılır.
  5. Otizmde bulunan kodonlar, farklı kontrol seçimleriyle tekrar sınanır. Bir bulgu tek
     bir rastgele kontrol seçimine bağlıysa güvenilmez.

Çıktılar: test/dco/cikti/dco_5_otizm, dco_5_sizofreni, dco_5_dogrudan
Not: HPO nadir hastalık odaklıdır. Şizofreninin yaygın risk genleri (GWAS) burada yoktur
ve liste 22q11 gibi sendromların sinir sistemiyle ilgisiz genlerini de içerir. Şizofreni
sonucu bu nedenle en zayıf kanıt düzeyindedir. Daha iyi gen seti için SCHEMA ve SFARI
gibi özel kaynaklar gerekir.
"""
from statistics import median

from scipy.stats import mannwhitneyu

from crn.cub import CUB, DifferentialCodonUsage
from dco_ortak import CIKTI, ILISKI, analiz, eslestir, gen_adi, gen_al, hazirla, paketle

k = hazirla()
otizm = gen_al(k, "Autistic behavior", association=ILISKI, evidence=["PCS", "TAS"], min_frequency=0.3)
sizofreni = gen_al(k, "Schizophrenia")
ortak = {gen_adi(g) for g in otizm} & {gen_adi(g) for g in sizofreni}
print(f"Otizm {len(otizm)} gen, şizofreni {len(sizofreni)} gen, ortak {len(ortak)}: {sorted(ortak)}\n")


def uzunluk_gc3(genler):
    tablo = CUB(genler).table(protein=False)
    return [s["L"] for s in tablo], [s["GC3s"] for s in tablo]


havuz_l, havuz_g = uzunluk_gc3(k.havuz)
print("Eşleştirmenin neden gerekli olduğu (medyan kodon sayısı, medyan GC3s):")
print(f"  {'Tüm kontrol havuzu':<26}{median(havuz_l):>7.0f}{median(havuz_g):>8.3f}")
for ad, genler in (("Otizm", otizm), ("Şizofreni", sizofreni)):
    l, g = uzunluk_gc3(genler)
    print(f"  {ad + ', eşleştirme öncesi':<26}{median(l):>7.0f}{median(g):>8.3f}"
          f"   havuza karşı p: L {mannwhitneyu(l, havuz_l).pvalue:.1e}, GC3s {mannwhitneyu(g, havuz_g).pvalue:.1e}")
    kl, kg = uzunluk_gc3(eslestir(genler, k.havuz)[0])
    print(f"  {ad + ' kontrolü, eşleşmiş':<26}{median(kl):>7.0f}{median(kg):>8.3f}"
          f"   gruba karşı p: L {mannwhitneyu(l, kl).pvalue:.2f}, GC3s {mannwhitneyu(g, kg).pvalue:.2f}")
print()

sonuc_otizm = analiz(k, "Otizm", otizm, "dco_5_otizm")
analiz(k, "Sizofreni", sizofreni, "dco_5_sizofreni")
analiz(k, "Otizm", [g for g in otizm if gen_adi(g) not in ortak], "dco_5_dogrudan",
       karsi=[g for g in sizofreni if gen_adi(g) not in ortak], kontrol_adi="Sizofreni")

# Duyarlılık: Otizm bulguları kontrol seçimi değişince de kalıyor mu?
kodonlar = [s["codon"] for s in sonuc_otizm.saglam]
if kodonlar:
    print(f"Duyarlılık: {', '.join(kodonlar)} için 10 farklı kontrol seçimi")
    gecici = CIKTI / "dco_5_otizm" / "duyarlilik.crn"
    gecen = dict.fromkeys(kodonlar, 0)
    for tohum in range(10):
        kontrol, _ = eslestir(otizm, k.havuz, tohum)
        paketle({"Otizm": otizm, "kontrol": kontrol}, gecici)
        dcu = DifferentialCodonUsage.from_crn(gecici, "Otizm", "kontrol")
        for satir in dcu.permutation("family", n=2000, seed=2026):
            if satir["codon"] in gecen:
                gecen[satir["codon"]] += satir["perm_q"] < .05
    gecici.unlink()
    for kodon, n in gecen.items():
        print(f"  {kodon}: 10 seçimin {n} tanesinde permütasyon q<0.05")
