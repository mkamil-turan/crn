"""Örnek 7: Genomik bölge sorgusu ve yapısal varyantların gen sayımına etkisi.

Senaryo: 22q11.2 delesyon sendromu bölgesinde çalışıyorsunuz. Bölgedeki varyantları
listeleyecek, kopya sayısı değişimlerini (CNV) tek nükleotit değişimlerinden ayıracak
ve gen düzeyindeki sayımlarda büyük delesyonların nasıl bir şişme yarattığını
göreceksiniz.

Bu örnekte:
  1. Bölge kromozom ve GRCh38 koordinatlarıyla sorgulanır.
  2. Varyant tipleri ve gen sayısı dağılımı çıkarılır.
  3. Patojen CNV'ler uzunluk ve kapsadıkları gen sayısıyla listelenir.
  4. ClinVar'ın büyük CNV'lerde yalnızca doz duyarlı genleri listelediği gösterilir.
  5. max_genes=1 (varsayılan) ile max_genes=None arasındaki gen sayımı farkı gösterilir.

Çıktılar: test/clinvar/cikti/clinvar_orn_7/22q11_cnv.tsv
Not: Bölge filtresi seçili derlemenin (burada GRCh38) koordinatlarına uygulanır.
GRCh37 koordinatlarıyla çalışmak için Clinvar(assembly="GRCh37") kullanılır.
"""
from collections import Counter
from pathlib import Path

from crn.fetch import Clinvar

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_7"
CIKTI.mkdir(parents=True, exist_ok=True)

BOLGE = {"chromosome": "22", "start": 18_900_000, "end": 21_500_000}   # LCR22 A-D, GRCh38
cv = Clinvar(verbose=False)

hepsi = cv.variants(max_genes=None, **BOLGE)
print(f"chr22:{BOLGE['start']:,}-{BOLGE['end']:,} ({cv.assembly}): {len(hepsi):,} varyant")
print("Tip:", dict(Counter(r["type"] for r in hepsi).most_common(6)))
gen_sayisi = Counter(r["n_genes"] for r in hepsi)
print(f"Tek genli: {gen_sayisi[1]:,}, 10'dan çok gen kapsayan: "
      f"{sum(n for k, n in gen_sayisi.items() if k > 10):,}\n")

cnv = [r for r in cv.variants(max_genes=None, type="copy number loss",
                              significance=["pathogenic", "likely_pathogenic"], **BOLGE)]
cnv.sort(key=lambda r: -(r["stop"] - r["start"]))
print(f"Patojen kopya kaybı: {len(cnv)}; en uzun beşi:")
for r in cnv[:5]:
    print(f"  {(r['stop'] - r['start']) / 1e6:5.2f} Mb  {r['n_genes']:>3} gen  "
          f"{r['chromosome']}:{r['start']:,}-{r['stop']:,}  doz duyarlı: {r['gene_symbol'] or '-'}")

print("\nGen düzeyinde patojen sayısı, max_genes etkisi:")
for gen in ("TBX1", "COMT", "HIRA", "CRKL"):
    tek = cv.genes(gene=gen, significance="pathogenic")
    tum = cv.genes(gene=gen, significance="pathogenic", max_genes=None)
    print(f"  {gen:<5} max_genes=1: {tek[0]['variants'] if tek else 0:>4}   "
          f"max_genes=None: {tum[0]['variants'] if tum else 0:>4}")
print("Sınırsız sayımda her büyük delesyon, listelediği doz duyarlı genlerin (TBX1, CRKL)\n"
      "her birine bir patojen varyant ekler;")
print("gen-hastalık ilişkisini değerlendirirken varsayılan (max_genes=1) kullanılmalıdır.")

cv.write_tsv(cnv, CIKTI / "22q11_cnv.tsv")
print(f"\nYazıldı: {CIKTI / '22q11_cnv.tsv'}")
