"""Örnek 4: Çok kayıtlı FASTA için kalite kontrolü.

Senaryo: Bir dizileme merkezinden çok sayıda contig içeren FASTA dosyası aldınız.
Analizden önce kısa ve belirsiz bazı (N) yüksek kayıtları ayıklamak, sonucu
raporlamak istiyorsunuz.

Bu örnekte:
  1. Sentetik bir FASTA üretilir (.fasta).
  2. Dosya okunur; her kayıt için uzunluk, GC yüzdesi ve N yüzdesi hesaplanır.
  3. Eşiklere göre her kaydın geçip geçmediği ve nedeni belirlenir.
  4. Geçenler yeni FASTA'ya, tüm sonuçlar TSV raporuna yazılır.
  5. Bozuk bir dosyanın okunurken hangi satırda yakalandığı gösterilir.

Çıktılar: test/sequence/cikti/seq_orn_4/
Not: Eşikler örnek amaçlıdır, kendi deneyinize göre değiştirin.
"""
import csv
import random
from pathlib import Path

from crn import DNA, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_4"
CIKTI.mkdir(parents=True, exist_ok=True)

EN_KISA = 500        # nt
EN_FAZLA_N = 5.0     # yüzde
GC_ARALIGI = (30.0, 70.0)  # yüzde

# 1. Sentetik veri: (ad, uzunluk, GC oranı, N oranı).
rastgele = random.Random(7)
profil = [("contig_1", 2400, 0.42, 0.00), ("contig_2", 1800, 0.55, 0.01),
          ("contig_3", 320, 0.48, 0.00), ("contig_4", 1500, 0.46, 0.18),
          ("contig_5", 3100, 0.81, 0.00), ("contig_6", 950, 0.38, 0.02)]


def uret(uzunluk: int, gc: float, n: float) -> str:
    agirlik = [(1 - n) * (1 - gc) / 2, (1 - n) * gc / 2, (1 - n) * gc / 2, (1 - n) * (1 - gc) / 2, n]
    return "".join(rastgele.choices("ACGTN", weights=agirlik, k=uzunluk))


girdi = CIKTI / "ham_contigler.fasta"
write_fasta([DNA(uret(u, g, n), name=ad) for ad, u, g, n in profil], girdi)


# 2-3. Okuma ve değerlendirme.
def degerlendir(kayit: DNA) -> dict:
    uzunluk = len(kayit)
    n_sayisi = kayit.seq.count("N")
    belirli = uzunluk - n_sayisi
    gc = 100 * (kayit.seq.count("G") + kayit.seq.count("C")) / belirli if belirli else 0.0
    n_yuzde = 100 * n_sayisi / uzunluk

    nedenler = []
    if uzunluk < EN_KISA:
        nedenler.append("kısa")
    if n_yuzde > EN_FAZLA_N:
        nedenler.append("çok N")
    if not GC_ARALIGI[0] <= gc <= GC_ARALIGI[1]:
        nedenler.append("GC aykırı")
    return {"ad": kayit.name, "uzunluk": uzunluk, "gc_yuzde": round(gc, 1),
            "n_yuzde": round(n_yuzde, 1), "durum": "GEÇTİ" if not nedenler else "ELENDİ",
            "neden": ", ".join(nedenler)}


kayitlar = DNA.from_fasta_file(girdi)
rapor = [degerlendir(k) for k in kayitlar]

print(f"{'Kayıt':<12}{'nt':>6}{'GC%':>7}{'N%':>7}  Durum")
for r in rapor:
    print(f"{r['ad']:<12}{r['uzunluk']:>6}{r['gc_yuzde']:>7}{r['n_yuzde']:>7}  "
          f"{r['durum']} {r['neden']}".rstrip())

# 4. Sonuçları yazma.
gecenler = [k for k, r in zip(kayitlar, rapor) if r["durum"] == "GEÇTİ"]
write_fasta(gecenler, CIKTI / "temiz_contigler.fasta")
with open(CIKTI / "kalite_raporu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(rapor[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(rapor)

toplam = sum(len(k) for k in kayitlar)
kalan = sum(len(k) for k in gecenler)
print(f"\n{len(gecenler)}/{len(kayitlar)} kayıt geçti; "
      f"{kalan:,}/{toplam:,} nt korundu (%{100 * kalan / toplam:.1f})")

# 5. Bozuk dosya: DNA olarak okunurken geçersiz harf, satır numarasıyla bildirilir.
bozuk = CIKTI / "bozuk.fasta"
bozuk.write_text(">iyi_kayit\nATGCATGC\n>kotu_kayit\nATGC!!GC\n", encoding="utf-8")
try:
    DNA.from_fasta_file(bozuk)
except ValueError as hata:
    print("\nBozuk dosya yakalandı:", hata)

print("\nÇıktılar:", CIKTI)
