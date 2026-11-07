"""Örnek 13: Tüm MANE transkriptlerinde start, stop ve Kozak bağlamı.

Senaryo: İnsan mRNA'larında translasyonun başlangıç ve bitiş bağlamını, bütün
genomdan tek bir özet olarak çıkarmak istiyorsunuz.

Bu örnekte:
  1. Tüm MANE Select protein kodlayan transkriptler (RNA) alınır.
  2. Stop kodonu kullanımı sayılır.
  3. Kozak bağlamındaki -3 ve +4 konumları ölçülür.
  4. UTR uzunluklarının ortancaları hesaplanır.
  5. Ortak başlangıç bağlamı konsensüsü çıkarılır.

Çıktılar: test/gencode/cikti/gencode_orn_13/
"""
import csv
from collections import Counter
from pathlib import Path
from statistics import median

from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_13"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)
kayitlar = gc.transcripts(tag="MANE_Select", transcript_type="protein_coding")

baslangic, stoplar, utr5, utr3 = Counter(), Counter(), [], []
konumlar = {i: Counter() for i in range(-6, 7)}
for t in kayitlar:
    rna = gc.rna(t["transcript_id"])
    cds = rna.features.first("cds")
    if cds is None:
        continue
    bas, bit = cds
    baslangic[rna.seq[bas - 1:bas + 2]] += 1
    stoplar[rna.seq[bit - 3:bit]] += 1
    utr5.append(bas - 1)
    utr3.append(len(rna) - bit)
    if bas > 6:
        # -1 ATG'nin hemen öncesi, +1 A'dır (0 konumu kullanılmaz).
        for i, baz in zip(list(range(-6, 0)) + list(range(1, 7)), rna.seq[bas - 7:bas + 5]):
            konumlar[i][baz] += 1

n = sum(baslangic.values())
print(f"MANE Select transkripti: {n:,}")
print("Başlangıç kodonu:", ", ".join(f"{k} (%{100 * v / n:.2f})" for k, v in baslangic.most_common(3)))
print("Stop kodonu     :", ", ".join(f"{k} (%{100 * v / n:.1f})" for k, v in stoplar.most_common(3)))
print(f"Ortanca 5'UTR: {median(utr5):.0f} nt, 3'UTR: {median(utr3):.0f} nt")

print("\nBaşlangıç bağlamı (baz yüzdesi):")
print(f"{'Konum':>6}" + "".join(f"{b:>7}" for b in "ACGU") + "  Baskın")
satirlar = []
for i, sayac in konumlar.items():
    if i == 0:
        continue
    toplam = sum(sayac.values())
    yuzde = {b: 100 * sayac[b] / toplam for b in "ACGU"}
    baskin = max(yuzde, key=yuzde.get)
    satirlar.append({"konum": i, **{b: round(yuzde[b], 1) for b in "ACGU"}})
    print(f"{i:>6}" + "".join(f"{yuzde[b]:>7.1f}" for b in "ACGU") + f"  {baskin}")

print("\n-3 konumunda purin (A/G):", f"%{sum(100 * konumlar[-3][b] / sum(konumlar[-3].values()) for b in 'AG'):.1f}")
print("+4 konumunda G          :", f"%{100 * konumlar[4]['G'] / sum(konumlar[4].values()):.1f}")

with open(CIKTI / "kozak_profili.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)
print("Çıktılar:", CIKTI)
