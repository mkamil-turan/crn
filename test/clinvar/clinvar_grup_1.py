"""Grup 1: Hastalık kimliklerinden güvenilir gen listeleri.

Senaryo: Birkaç kalıtsal hastalık için "kanıtı sağlam genler" listesini hazırlayıp
başka bir programa (ör. bir panel tasarım aracına) vereceksiniz. Gruplar OMIM
kimliğiyle tanımlanır, yalnızca en az iki yıldızlı patojen/olası patojen varyantı
olan genler alınır.

Bu örnekte:
  1. Her hastalık için gen satırları alınır (gen başına P/LP sayısı ve en yüksek yıldız).
  2. Sembol listeleri metin dosyasına, ayrıntı tablosu TSV'ye yazılır.
  3. Aynı hastalıkta yıldız eşiğinin liste büyüklüğüne etkisi gösterilir.

Çıktılar: test/clinvar/cikti/clinvar_grup_1/
Not: Tek bir varyantla listeye giren gen, iyi kurulmuş bir gen-hastalık ilişkisi
anlamına gelmez; ClinGen gen-hastalık geçerliliği gibi kaynaklarla karşılaştırın.
"""
from pathlib import Path

from crn.fetch import Clinvar

CIKTI = Path(__file__).parent / "cikti" / "clinvar_grup_1"
CIKTI.mkdir(parents=True, exist_ok=True)

HASTALIKLAR = {
    "kistik_fibrozis": "OMIM:219700",
    "marfan": "OMIM:154700",
    "fenilketonuri": "OMIM:261600",
    "hcm": "MONDO:0005045",            # hipertrofik kardiyomiyopati
    "lynch": "MONDO:0005835",          # Lynch sendromu
}
FILTRE = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 2}
cv = Clinvar(verbose=False)

tablo = []
for ad, kimlik in HASTALIKLAR.items():
    genler = cv.genes(kimlik, **FILTRE)
    semboller = sorted(g["gene_symbol"] for g in genler)
    (CIKTI / f"{ad}.txt").write_text("\n".join(semboller) + "\n")
    tablo += [{"group": ad, "condition": kimlik, **g} for g in genler]
    ilk = ", ".join(f"{g['gene_symbol']}({g['pathogenic']})" for g in genler[:6])
    print(f"{ad:<16} {kimlik:<14} {len(genler):>3} gen  {ilk}")

cv.write_tsv(tablo, CIKTI / "genler.tsv")

print("\nYıldız eşiğinin etkisi (hcm, P/LP):")
for yildiz in range(4):
    genler = cv.genes(HASTALIKLAR["hcm"], significance=FILTRE["significance"], min_stars=yildiz)
    print(f"  en az {yildiz} yıldız: {len(genler):>3} gen")

print(f"\nYazıldı: {CIKTI}")
