"""Örnek 5: Kayan pencereyle GC profili ve CpG adası taraması.

Senaryo: Bir promotör bölgesinde metilasyon çalışması planlıyorsunuz. CpG adası
adaylarını bulmak ve GC içeriğinin dizi boyunca değişimini görmek istiyorsunuz.

Bu örnekte:
  1. Ortasında CG'ce zengin bir bölge bulunan sentetik DNA üretilir.
  2. sliding ile 200 nt pencereler 10 nt adımla çıkarılır.
  3. Her pencerede GC yüzdesi ve gözlenen/beklenen CpG oranı hesaplanır.
  4. Gardiner-Garden & Frommer (1987) ölçütünü geçen pencereler birleştirilir:
     uzunluk >= 200 nt, GC > %50, gözlenen/beklenen CpG > 0.6.
  5. Profil terminalde çizilir, sonuçlar TSV ve FASTA olarak kaydedilir.

Çıktılar: test/sequence/cikti/seq_orn_5/
"""
import csv
import random
from pathlib import Path

from crn import DNA, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_5"
CIKTI.mkdir(parents=True, exist_ok=True)

PENCERE = 200
ADIM = 10
GC_ESIK = 50.0
OE_ESIK = 0.6

# 1. Sentetik dizi: 2 kb AT zengini + 600 nt CpG zengini + 1.4 kb AT zengini.
rastgele = random.Random(11)
sol = "".join(rastgele.choices("ACGT", weights=[3, 2, 2, 3], k=2000))
ada = "".join(rastgele.choices(["CG", "GC", "A", "T", "G", "C"],
                               weights=[4, 2, 1, 1, 1, 1], k=420))[:600]
sag = "".join(rastgele.choices("ACGT", weights=[3, 2, 2, 3], k=1400))
promotor = DNA(sol + ada + sag, name="promotor_bolgesi")


# 3. Pencere başına ölçümler.
def gc_yuzde(pencere: DNA) -> float:
    dizi = pencere.seq
    return 100 * (dizi.count("G") + dizi.count("C")) / len(dizi)


def cpg_oe(pencere: DNA) -> float:
    dizi = pencere.seq
    c, g = dizi.count("C"), dizi.count("G")
    return dizi.count("CG") * len(dizi) / (c * g) if c and g else 0.0


# 2. Her eleman yine DNA'dır; bu yüzden pencereler üzerinde aynı metotlar çalışır.
pencereler = promotor.sliding(PENCERE, step_size=ADIM)
olcumler = [{"baslangic": i * ADIM + 1, "bitis": i * ADIM + PENCERE,
             "gc_yuzde": round(gc_yuzde(p), 1), "cpg_oe": round(cpg_oe(p), 3)}
            for i, p in enumerate(pencereler)]
print(f"{promotor.name}: {len(promotor)} nt -> {len(pencereler)} pencere "
      f"({type(pencereler[0]).__name__} sınıfından)")

# 4. Ölçütü geçen ardışık pencereleri birleştir. Pencerenin tamamı yerine merkezi
# kullanılır; yoksa adanın kenarına değen her pencere bölgeyi dışa doğru uzatır.
adalar = []
for olcum in olcumler:
    if olcum["gc_yuzde"] > GC_ESIK and olcum["cpg_oe"] > OE_ESIK:
        merkez = olcum["baslangic"] + PENCERE // 2
        if adalar and merkez <= adalar[-1][1] + ADIM:
            adalar[-1][1] = merkez
        else:
            adalar.append([merkez, merkez])

# 5a. Terminal profili.
print("\nGC profili (her satır 200 nt pencere, 200 nt aralıkla)")
for olcum in olcumler[::20]:
    isaret = " <- CpG adası" if any(a <= olcum["baslangic"] + PENCERE // 2 <= b for a, b in adalar) else ""
    cubuk = "#" * int(olcum["gc_yuzde"] / 2)
    print(f"{olcum['baslangic']:>5} {olcum['gc_yuzde']:>5.1f}% {cubuk}{isaret}")

# 5b. Sonuçları kaydet.
with open(CIKTI / "gc_profili.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(olcumler[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(olcumler)

print(f"\nBulunan CpG adası adayı: {len(adalar)}")
adalar_dna = []
for no, (bas, bit) in enumerate(adalar, start=1):
    bolge = promotor[bas - 1:bit]
    bolge.name = f"CpG_adasi_{no}"
    bolge.features.set("konum", (bas, bit))
    adalar_dna.append(bolge)
    print(f"  {bolge.name}: {bas}-{bit} ({len(bolge)} nt), "
          f"GC %{gc_yuzde(bolge):.1f}, O/E {cpg_oe(bolge):.2f}")
if adalar_dna:
    write_fasta(adalar_dna, CIKTI / "cpg_adalari.fasta")
print(f"Yerleştirilen bölge: {len(sol) + 1}-{len(sol) + len(ada)} "
      f"(sınırlar yaklaşık yarım pencere şaşabilir)")
print("\nÇıktılar:", CIKTI)
