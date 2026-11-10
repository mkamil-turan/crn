"""Örnek 2: Tek bir genin klinik varyant profili.

Senaryo: Ailesel Akdeniz ateşi (FMF) çalışan bir laboratuvardasınız. MEFV geninde
ClinVar'da ne kadar varyant olduğunu, bunların ne kadarının güvenilir biçimde
patojen sayıldığını ve hangi değişimlerin öne çıktığını görmek istiyorsunuz.

Bu örnekte:
  1. Genin bütün varyantları alınır; sınıf ve yıldız çapraz tablosu kurulur.
  2. Patojen/olası patojen varyantların protein etkisi dağılımı çıkarılır.
  3. En çok laboratuvarın bildirdiği patojen değişimler listelenir.
  4. Genin ilişkili olduğu hastalıklar sıralanır.
  5. Patojen varyantlar TSV olarak yazılır.

Çıktılar: test/clinvar/cikti/clinvar_orn_2/
Not: Bir varyant birden çok sınıfa girebilir (ör. "Pathogenic/Likely pathogenic"),
bu yüzden tablodaki satırların toplamı varyant sayısından büyük olabilir.
"""
from collections import Counter
from pathlib import Path

from crn.fetch import Clinvar
from crn.fetch.clinvar import SIGNIFICANCES

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_2"
CIKTI.mkdir(parents=True, exist_ok=True)

GEN = "MEFV"
PLP = ["pathogenic", "likely_pathogenic"]
cv = Clinvar(verbose=False)

hepsi = cv.variants(GEN)
print(f"{GEN}: ClinVar {cv.release} içinde {len(hepsi)} varyant (tek genli)\n")

print(f"{'sınıf':<18}" + "".join(f"{'*' * y or '0':>6}" for y in range(5)) + f"{'toplam':>8}")
for sinif in SIGNIFICANCES:
    satirlar = [r for r in hepsi if sinif in r["classes"].split("|")]
    if not satirlar:
        continue
    yildiz = Counter(r["stars"] for r in satirlar)
    print(f"{sinif:<18}" + "".join(f"{yildiz.get(y, 0):>6}" for y in range(5)) + f"{len(satirlar):>8}")

plp = cv.variants(GEN, significance=PLP)
guvenilir = [r for r in plp if r["stars"] >= 1]
print(f"\nPatojen/olası patojen: {len(plp)}, en az 1 yıldızlı: {len(guvenilir)}")
print("Protein etkisi:", dict(Counter(r["effect"] or "(yazım yok)" for r in plp).most_common()))

print("\nEn çok gönderen laboratuvarı olan patojen değişimler:")
for r in sorted(guvenilir, key=lambda r: -r["submitters"])[:8]:
    print(f"  p.{r['hgvs_p'] or '?':<12} {r['hgvs_c']:<16} {'*' * r['stars']:<4} "
          f"{r['submitters']:>3} gönderen  rs{r['rs'] or '-'}")

print(f"\n{GEN} ile ilişkili hastalıklar (bilgilendirici adlar):")
for satir in cv.conditions(gene=GEN, limit=6):
    print(f"  {satir['variants']:>5} varyant  {satir['condition_name']}")

cv.write_tsv(plp, CIKTI / f"{GEN}_patojen.tsv")
print(f"\nYazıldı: {CIKTI / f'{GEN}_patojen.tsv'}")
