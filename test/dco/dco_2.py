"""Örnek 2: Poliglutamin hastalıkları ve CAG tekrarı.

Klinik arka plan: Huntington hastalığı, SCA1/2/3/6/7/17, DRPLA ve spinal bulbar
müsküler atrofi aynı mekanizmayı paylaşır: proteinde glutamin (Q) dizisini kodlayan
CAG tekrarı genişler, belli bir eşiği aşınca hastalık başlar. Tekrar uzadıkça
hastalık genellikle daha erken başlar (antisipasyon).

Glutaminin iki kodonu vardır: CAG ve CAA. İkisi aynı amino asidi verir, yani protein
düzeyinde fark yoktur. Ama DNA düzeyinde tekrar içine giren CAA "kesinti" gibi davranır
ve tekrarın kararlılığını artırır. ATXN2 ve HTT'de kesinti kaybının hastalık
riskini ve başlangıç yaşını etkilediği bildirilmiştir. Bu yüzden burada sinonim kodon
seçimi tıbbi açıdan doğrudan anlamlıdır.

Soru: Referans genomdaki tekrar ne kadar uzun, eşiğe ne kadar yakın, içinde CAA var
mı? Ve hastalık genlerinin CAG tercihi yalnızca tekrardan mı geliyor, yoksa genin
geri kalanında da mı var?

Bu örnekte:
  1. Her genin en uzun glutamin bölgesi kodon kodon bulunur (tek kodonluk kesintilerle).
  2. Tekrar yapısı, örneğin (CAG)13-CAA-(CAG)8 biçiminde yazılır.
  3. Hastalık eşiğine kalan mesafe hesaplanır.
  4. Genin CAG payı tekrar dahil ve hariç hesaplanır, eşleştirilmiş kontrolle kıyaslanır.

Çıktılar: test/dco/cikti/dco_2/
Not: Tekrar uzunluğu burada kesintiler dahil kodon sayısıdır; laboratuvar raporları farklı
sayabilir (ör. yalnızca CAG). Eşikler kaynaklarda sık anılan yaklaşık değerlerdir ve kurala göre değişir, yalnızca
örnek içindir. Referans genom sağlıklı bir aleli gösterir, hastalık aleli değildir. Genin
kullanılan transkripti tekrarın kodlandığı transkript olmayabilir. Grup yalnızca 9 gen
olduğu için istatistiksel güç düşüktür.
"""
import re
from collections import Counter
from itertools import groupby

from scipy.stats import mannwhitneyu

from crn.cub import CUB
from crn.cub.differential import fisher_log10_p
from dco_ortak import CIKTI, Hpo, eslestir, gen_adi, gen_al, hazirla, paketle

CIKTI_YOLU = CIKTI / "dco_2"
CIKTI_YOLU.mkdir(parents=True, exist_ok=True)

# gen -> (hastalık, patojenik sayılan yaklaşık tekrar eşiği)
HASTALIKLAR = {
    "HTT": ("Huntington", 36), "ATXN1": ("SCA1", 39), "ATXN2": ("SCA2", 33), "ATXN3": ("SCA3", 55),
    "CACNA1A": ("SCA6", 20), "ATXN7": ("SCA7", 37), "TBP": ("SCA17", 43),
    "ATN1": ("DRPLA", 48), "AR": ("SBMA", 38),
}


def tekrar_bolgesi(kodonlar):
    """En uzun glutamin bölgesi: (başlangıç, bitiş). Tek kodonluk kesintilere izin verir.

    Kesinti tıbbi açıdan önemlidir. ATXN1'de CAT (histidin), ATXN3'te AAG (lizin)
    tekrarı böler. Yalnızca kesintisiz dizi aranırsa ATXN1 tekrarı yarıya iner ve
    kesintinin kendisi gözden kaçar.
    """
    harita = "".join("Q" if kodon in ("CAG", "CAA") else "." for kodon in kodonlar)
    bolgeler = [m.span() for m in re.finditer(r"Q+(?:\.Q+)*", harita)]
    return max(bolgeler, key=lambda b: b[1] - b[0], default=(0, 0))


def yapi(kodonlar):
    """Tekrarı okunur yazar: CAG CAG CAA CAG -> (CAG)2-CAA-CAG."""
    parcalar = []
    for kodon, grup in groupby(kodonlar):
        n = len(list(grup))
        parcalar.append(kodon if n == 1 else f"({kodon}){n}")
    return "-".join(parcalar)


def cag_payi(kodonlar):
    cag, caa = kodonlar.count("CAG"), kodonlar.count("CAA")
    return cag, caa, cag / (cag + caa) if cag + caa else 0.0


k = hazirla()
# association filtresi kullanılmaz: HPO, SCA2'yi "POLYGENIC" etiketler ve ATXN2 düşer.
hasta = gen_al(k, gene=list(HASTALIKLAR))
kontrol, denge = eslestir(hasta, k.havuz)
paketle({"poliQ": hasta, "kontrol": kontrol}, CIKTI_YOLU / "poliq_kontrol.crn")
print(f"Eşleştirme dengesi: uzunluk {denge['uzunluk']:.3f}, GC {denge['GC']:.3f}\n")

cub = CUB(hasta)
print(f"{'Gen':<9}{'Hastalık':<12}{'Tekrar':>7}{'Eşik':>6}{'Fark':>6}  {'Kesintiler':<16}Yapı")
satirlar, disarida = [], []
for gen in sorted(hasta, key=lambda g: HASTALIKLAR[gen_adi(g)][1]):
    ad = gen_adi(gen)
    kodonlar = cub.codons(gen.name)
    bas, bit = tekrar_bolgesi(kodonlar)
    tekrar = kodonlar[bas:bit]
    kesinti = Counter(c for c in tekrar if c != "CAG")
    kesinti_yazi = ", ".join(f"{c}x{n}" for c, n in kesinti.most_common()) or "yok"
    hastalik, esik = HASTALIKLAR[ad]
    geri_kalan = kodonlar[:bas] + kodonlar[bit:]
    tum_cag, tum_caa, pay_tum = cag_payi(kodonlar)
    cag_d, caa_d, pay_dis = cag_payi(geri_kalan)
    disarida.append((cag_d, caa_d))
    satirlar.append({"gen": ad, "hastalik": hastalik, "q_tekrari": len(tekrar), "esik": esik,
                     "esige_kalan": esik - len(tekrar), "kesintiler": kesinti_yazi,
                     "yapi": yapi(tekrar), "cag_payi_tum": round(pay_tum, 3),
                     "cag_payi_tekrar_haric": round(pay_dis, 3)})
    print(f"{ad:<9}{hastalik:<12}{len(tekrar):>7}{esik:>6}{esik - len(tekrar):>6}  "
          f"{kesinti_yazi:<16}{yapi(tekrar)[:64]}")
Hpo.write_tsv(satirlar, CIKTI_YOLU / "tekrar_yapilari.tsv")

# Tekrar çıkarıldığında CAG tercihi kalıyor mu? Kontrol genleri en uzun Q bölgelerinden
# aynı biçimde arındırılır, böylece iki grup eşit koşulda karşılaştırılır.
kontrol_cub = CUB(kontrol)
k_tum, k_dis = [], []
for gen in kontrol:
    kodonlar = kontrol_cub.codons(gen.name)
    bas, bit = tekrar_bolgesi(kodonlar)
    k_tum.append(cag_payi(kodonlar))
    k_dis.append(cag_payi(kodonlar[:bas] + kodonlar[bit:]))


def karsilastir(ad, hasta_paylari, kontrol_paylari, h_toplam, k_toplam):
    test = mannwhitneyu(hasta_paylari, kontrol_paylari, alternative="two-sided")
    a, b = h_toplam
    c, d = k_toplam
    print(f"{ad:<16}poliQ medyan {sorted(hasta_paylari)[len(hasta_paylari) // 2]:.3f}, "
          f"kontrol medyan {sorted(kontrol_paylari)[len(kontrol_paylari) // 2]:.3f}, "
          f"gen düzeyi Mann-Whitney p={test.pvalue:.3f}, havuzlu Fisher p={10 ** fisher_log10_p(a, b, c, d):.1e}")


print("\nCAG payı (CAG / CAG + CAA):")
tum = [cag_payi(cub.codons(g.name)) for g in hasta]
karsilastir("Tekrar dahil", [x[2] for x in tum], [x[2] for x in k_tum],
            (sum(x[0] for x in tum), sum(x[1] for x in tum)), (sum(x[0] for x in k_tum), sum(x[1] for x in k_tum)))
karsilastir("Tekrar hariç", [c / (c + a) for c, a in disarida if c + a], [x[2] for x in k_dis],
            (sum(c for c, _ in disarida), sum(a for _, a in disarida)),
            (sum(x[0] for x in k_dis), sum(x[1] for x in k_dis)))
print(f"\nÇıktılar: {CIKTI_YOLU}")
