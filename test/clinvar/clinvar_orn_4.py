"""Örnek 4: Hastalık adından kesin hastalık kimliğine.

Senaryo: "Kistik fibrozis" ve "Marfan sendromu" genlerini çıkaracaksınız. ClinVar'daki
hastalık adları gönderildiği biçimde tutulduğundan aynı hastalık birden çok adla
geçer, kısa bir ad parçası da ilgisiz hastalıklara uyar. Yayımlanacak bir analizde
grubu OMIM ya da MONDO kimliğiyle tanımlamak gerekir.

Bu örnekte:
  1. Ad parçasıyla hastalık aranır, kimlikler ve varyant sayıları görülür.
  2. Çok geniş bir arama ("fever") ile dar bir arama karşılaştırılır.
  3. Düzenli ifadeyle tam ad eşleşmesi yapılır.
  4. Ad parçası, tam ad ve kimlikle yapılan sorguların gen listeleri karşılaştırılır.

Çıktı: yalnızca ekran.
Not: Bir kimlik yazımı "OMIM:219700" biçimindedir; MONDO, MedGen ve Orphanet
kimlikleri de aynı biçimde kullanılabilir.
"""
import re

from crn.fetch import Clinvar

cv = Clinvar(verbose=False)
PLP = ["pathogenic", "likely_pathogenic"]

for ad in ("cystic fibrosis", "Marfan"):
    print(f"'{ad}' araması:")
    for satir in cv.conditions(ad, limit=5):
        print(f"  {satir['variants']:>6} varyant {satir['genes']:>3} gen  {satir['condition_name']}")
        print(f"         {satir['condition_ids']}")
    print()

genis = cv.conditions("fever")
print(f"'fever' {len(genis)} hastalığa uyar; ilk beşi:",
      "; ".join(s["condition_name"] for s in genis[:5]), "\n")

tam = cv.conditions(re.compile(r"^Marfan syndrome$", re.IGNORECASE))
print("Düzenli ifade ile tam ad:", [(s["condition_name"], s["variants"]) for s in tam], "\n")

print("Aynı soru, üç sorgu biçimi (patojen/olası patojen, en az 1 yıldız):")
for etiket, kosul in [("ad parçası 'Marfan'", "Marfan"),
                      ("tam ad", re.compile(r"^Marfan syndrome$", re.IGNORECASE)),
                      ("OMIM:154700", "OMIM:154700")]:
    genler = cv.genes(kosul, significance=PLP, min_stars=1)
    ilk = ", ".join(f"{g['gene_symbol']}({g['pathogenic']})" for g in genler[:6])
    print(f"  {etiket:<22} {len(genler):>3} gen  {ilk}")

print("\nAd parçası Marfan'a benzeyen sendromları (Marfanoid, Loeys-Dietz eşlikli adlar)")
print("de toplar; kimlik yalnızca OMIM 154700 kaydına bağlı varyantları verir.")
