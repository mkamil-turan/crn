"""Örnek 1: Spinoserebellar ataksi (SCA) genlerinde kodon kullanımı.

Klinik arka plan: SCA, beyincik ve bağlantılarını etkileyen, ilerleyici yürüme ve
konuşma bozukluğuyla seyreden kalıtsal ataksi grubudur. Altta yatan genler çok
çeşitlidir: iyon kanalları, poliglutamin proteinleri, mitokondri proteinleri.

Soru: Bu kadar çeşitli genler ortak bir kodon kullanım özelliği taşıyor mu? Örneğin
beyincik nöronlarının sürekli çalışması gerektiği düşünülürse, bu genler tRNA havuzuna
daha iyi uyuyor olabilir mi? Ya da tam tersine yavaş kodonlar mı taşıyor?

Akış (ayrıntı için dco_ortak.py başındaki açıklama):
  1. Adında "spinocerebellar ataxia" geçen hastalıkların genleri alınır.
  2. Her gene uzunluk ve genomik GC'si en yakın bir kontrol geni seçilir.
  3. İki grup tek CRN dosyasına yazılır.
  4. Kodonlar küresel, aile ve alt aile düzeyinde test edilir.
  5. Aile düzeyi bulgular gen etiketleri karıştırılarak sınanır.
  6. ENC, GC3s, CAI, tAI gibi parametreler gen başına karşılaştırılır.
  7. Klinik profil ve kurallı yorum yazılır, uç genler listelenir.

Çıktılar: test/dco/cikti/dco_1/
"""
import re

from dco_ortak import ILISKI, analiz, gen_adi, gen_al, hazirla

k = hazirla()
hasta = gen_al(k, disease=re.compile("spinocerebellar ataxia", re.I), association=ILISKI)
sonuc = analiz(k, "SCA", hasta, "dco_1")

# Grup ortalamasının altında kalan tek tek genler de klinik açıdan ilginçtir: tRNA
# havuzuna en az uyumlu SCA genleri hangileri? Bu yalnızca hipotez üretir.
adlar = {g.name: gen_adi(g) for g in sonuc.hasta}
print("tRNA havuzuna en az uyan 6 SCA geni (tAI), yanına genin uzunluğu ve ENC:")
for s in sorted(sonuc.tablo_a, key=lambda x: x["tAI"])[:6]:
    print(f"  {adlar[s['gene']]:<10} tAI {s['tAI']:.3f}  L {s['L']:>5}  ENC {s['ENC']:.1f}  GC3s {s['GC3s']:.2f}")
