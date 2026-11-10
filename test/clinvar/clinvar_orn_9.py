"""Örnek 9: Kendi gen listenizi ClinVar özetiyle işaretlemek.

Senaryo: Bir RNA-seq ya da kodon kullanımı analizinden çıkan genleriniz var ve her
birinin klinik önemini sonuç tablosunda göstermek istiyorsunuz: hangi hastalıklarla
ilişkili, kaç patojen varyantı var. Bilgiyi genlerin CRN dosyasına yazarsanız sonraki
analizlerde de taşınır.

Bu örnekte:
  1. Genler GENCODE'dan ana transkriptleriyle alınır.
  2. Clinvar.annotate ile features'a ClinVar özeti yazılır (iki yıldız ve üstü).
  3. Liste CRN dosyasına yazılır, geri okunup değerler denetlenir.
  4. Özet tablo TSV olarak çıkarılır.

Çıktılar: test/clinvar/cikti/clinvar_orn_9/
Not: ClinVar'da hiç varyantı olmayan gen KeyError verir; bu örnekte yakalanıp
"kayıt yok" olarak işaretlenir. Bu, genin hastalıkla ilgisiz olduğu anlamına gelmez.
"""
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Clinvar, Gencode

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_9"
CIKTI.mkdir(parents=True, exist_ok=True)

GENLER = ["MEFV", "GPX4", "SOD1", "HTT", "ACTB", "RPL3", "OR1A1", "GAPDH"]
cv = Clinvar(verbose=False)
gc = Gencode(verbose=False)

isaretli, tablo = [], []
for sembol in GENLER:
    gen = gc.canonical(sembol)
    try:
        cv.annotate(gen, min_stars=2)
    except KeyError:
        gen.features.set("clinvar_release", cv.release)
        gen.features.set("clinvar_variants", 0)
    isaretli.append(gen)
    hastaliklar = gen.features.get("clinvar_condition") or []
    tablo.append({"gene": sembol,
                  "variants_2star": gen.features.first("clinvar_variants"),
                  "pathogenic_2star": gen.features.first("clinvar_pathogenic") or 0,
                  "conditions": "; ".join(hastaliklar[:3])})

dosya = CIKTI / "isaretli.crn"
write_crn_list(isaretli, dosya)
geri = read_crn_list(dosya)
print(f"{len(geri)} gen yazıldı ve okundu: {dosya}")
print("ClinVar sürümü korunmuş:", {g.features.first("clinvar_release") for g in geri}, "\n")

print(f"{'gen':<7}{'≥2* varyant':>12}{'P/LP':>6}  hastalıklar")
for r in tablo:
    print(f"{r['gene']:<7}{r['variants_2star']:>12}{r['pathogenic_2star']:>6}  "
          f"{r['conditions'][:70] or '-'}")

cv.write_tsv(tablo, CIKTI / "isaretli.tsv")
print(f"\nYazıldı: {CIKTI / 'isaretli.tsv'}")
