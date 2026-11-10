"""Örnek 1: ClinVar'ı kurmak ve elinizdeki sürümü tanımak.

Senaryo: Bir varyant yorumlama ya da gen grubu çalışmasına başlıyorsunuz. Yöntem
bölümüne hangi ClinVar sürümünü, hangi derlemeyi kullandığınızı ve verinin genel
dağılımını yazmanız gerekecek.

Bu örnekte:
  1. Veri ilk çağrıda indirilir ve SQLite veritabanına hazırlanır (tqdm ilerleme çubuğu).
  2. Sürüm, derleme, dosya özeti (SHA-256) ve kayıt sayıları okunur.
  3. Klinik anlam sınıflarının dağılımı yüzde olarak verilir.
  4. İnceleme durumlarının yıldız karşılıkları listelenir.

Çıktı: yalnızca ekran.
Not: İlk çalıştırma yaklaşık 450 MB indirir ve hazırlık birkaç dakika sürer
(GRCh38 için 1 GB veritabanı). Sonraki çalıştırmalar internet istemez.
"""
import json

from crn.fetch import Clinvar
from crn.fetch.clinvar import REVIEW_STARS, SIGNIFICANCES

cv = Clinvar()                       # sürüm verilmezse sys/ içindeki en yeni hazır sürüm
ozet = cv.summary()                  # veri yoksa burada indirilir ve hazırlanır
manifest = json.loads((cv.release_dir / "manifest.json").read_text())

print(f"ClinVar {ozet['release']} ({ozet['assembly']})")
print(f"  Kaynak : {manifest['source']}")
print(f"  SHA-256: {manifest['sha256']}")
print(f"  Veritabanı: {cv.database}  ({cv.database.stat().st_size / 1e9:.2f} GB)")
print(f"  {ozet['variants']:,} varyant, {ozet['genes']:,} gen, {ozet['conditions']:,} hastalık\n")

print("Klinik anlam sınıfları (birleşik değerler birden çok sınıfa girer):")
for sinif in SIGNIFICANCES:
    sayi = ozet[sinif]
    print(f"  {sinif:<18} {sayi:>10,}  %{100 * sayi / ozet['variants']:5.1f}")

print("\nTanınmayan terim:", ozet["unknown_significance"] or "yok")

print("\nİnceleme durumu -> yıldız:")
for durum, yildiz in sorted(REVIEW_STARS.items(), key=lambda x: -x[1]):
    print(f"  {'*' * yildiz:<4} {yildiz}  {durum}")

print("\nYöntem cümlesi için:")
print(f"  ClinVar {ozet['release']} aylık arşivi (variant_summary, SHA-256 "
      f"{manifest['sha256'][:12]}...) {ozet['assembly']} koordinatlarıyla kullanılmıştır.")
