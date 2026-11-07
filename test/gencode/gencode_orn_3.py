"""Örnek 3: Protein kodlayan genlerde genom ölçeğinde derin sorgular.

Senaryo: Protein kodlayan genlerin istatistiğini çıkarmak istiyorsunuz: en uzun
genler, en büyük proteinler, en çok ekzonlu genler ve kromozomların gen yoğunluğu.

Bu örnekte:
  1. protein_coding genler sorgulanır, genomik uzunluğa göre sıralanır.
  2. Kromozom başına gen yoğunluğu (gen/Mb) hesaplanır.
  3. MANE Select transkriptleri üzerinden protein ve ekzon istatistikleri çıkarılır.
  4. En uzun 5 protein (PRT) tek CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_3/
"""
from collections import Counter
from pathlib import Path
from statistics import median

from crn import write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_3"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)
genler = gc.genes(gene_type="protein_coding")
print(f"Protein kodlayan gen: {len(genler):,}")

uzunluklar = sorted(genler, key=lambda g: g["end"] - g["start"], reverse=True)
print(f"Ortanca genomik uzunluk: {median(g['end'] - g['start'] + 1 for g in genler):,.0f} nt")
print("\nEn uzun 5 gen:")
for g in uzunluklar[:5]:
    print(f"  {g['gene_name']:<10}{g['chromosome']:<7}{(g['end'] - g['start'] + 1) / 1e6:>6.2f} Mb")

kromozom_boyu = {}
with open(gc.release_dir / "genome.fai") as dosya:
    for satir in dosya:
        ad, boy = satir.split("\t")[:2]
        kromozom_boyu[ad] = int(boy)
sayac = Counter(g["chromosome"] for g in genler)
yogunluk = {k: 1e6 * sayac[k] / kromozom_boyu[k] for k in sayac if k.startswith("chr") and k != "chrM"}
sirali = sorted(yogunluk.items(), key=lambda x: x[1], reverse=True)
print("\nEn yoğun kromozomlar :", ", ".join(f"{k} ({v:.1f})" for k, v in sirali[:3]), "gen/Mb")
print("En seyrek kromozomlar:", ", ".join(f"{k} ({v:.1f})" for k, v in sirali[-3:]), "gen/Mb")

mane = gc.transcripts(tag="MANE_Select", transcript_type="protein_coding")
gen_adi = {g["gene_id"]: g["gene_name"] for g in genler}
kayitlar = []
for t in mane:
    protein = gc.protein(t["transcript_id"])
    kayitlar.append((gen_adi.get(t["gene_id"], "?"), protein, len(gc.details(t["transcript_id"])["exons"])))

print(f"\nMANE Select proteini: {len(kayitlar):,}")
print(f"Ortanca protein: {median(len(p) for _, p, _ in kayitlar):.0f} aa, "
      f"ortanca ekzon: {median(e for _, _, e in kayitlar):.0f}")
print(f"Tek ekzonlu gen: {sum(e == 1 for _, _, e in kayitlar):,}")

en_uzun = sorted(kayitlar, key=lambda k: len(k[1]), reverse=True)[:5]
print("\nEn uzun proteinler:")
for ad, p, e in en_uzun:
    print(f"  {ad:<10}{len(p):>8,} aa{e:>6} ekzon")
print("En çok ekzonlu:", ", ".join(f"{a} ({e})" for a, _, e in sorted(kayitlar, key=lambda k: k[2], reverse=True)[:3]))

paket = []
for ad, p, _ in en_uzun:
    p.name = f"{ad}_MANE"
    paket.append(p)
write_crn_list(paket, CIKTI / "en_uzun_proteinler.crn")
print("\nÇıktılar:", CIKTI)
