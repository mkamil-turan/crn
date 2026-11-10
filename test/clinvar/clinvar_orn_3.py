"""Örnek 3: Hasta raporundaki varyantları ClinVar'da bulmak.

Senaryo: Bir genetik test raporunda birkaç varyant var. Kimi HGVS kodlama yazımıyla
(c.), kimi protein yazımıyla (p.), kimi yalnızca dbSNP numarasıyla verilmiş. Her
birinin ClinVar sınıfını, inceleme düzeyini ve son değerlendirme tarihini, ayrıca
gnomAD ya da IGV'de bakmak için VCF koordinatını çıkaracaksınız.

Bu örnekte:
  1. Varyant gene göre alınır, c., p. ya da rs ile eşleştirilir.
  2. Sınıf, yıldız, gönderen sayısı, tarih ve hastalıklar yazılır.
  3. VCF biçiminde konum (kromozom, konum, ref, alt) üretilir.
  4. Bulunamayan girdi açıkça bildirilir.

Çıktılar: test/clinvar/cikti/clinvar_orn_3/rapor.tsv
Not: Protein yazımı tek başına bir varyantı belirlemeyebilir; aynı amino asit
değişimine birden çok nükleotit değişimi yol açabilir. Bu durumda bütün eşleşmeler
listelenir. CFTR p.Phe508del ClinVar'da iki kayıtla geçer (c.1521_1523del ve
c.1522_1524del): aynı delesyon, farklı hizalanmış yazım. HBB c.20A>T (orak hücre)
eski adlandırmada Glu6Val, HGVS'te Glu7Val'dir.
"""
from pathlib import Path

from crn.fetch import Clinvar
from crn.fetch.clinvar import UNINFORMATIVE_CONDITIONS

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_3"
CIKTI.mkdir(parents=True, exist_ok=True)

RAPOR = [
    ("CFTR", "p.Phe508del"),
    ("MEFV", "p.Met694Val"),
    ("HBB", "c.20A>T"),
    ("BRCA1", "c.68_69del"),
    ("PAH", "rs5030858"),
    ("APOE", "p.Cys130Arg"),
    ("LDLR", "c.9999A>G"),           # rapor hatası: böyle bir varyant yok
]
cv = Clinvar(verbose=False)


def eslesir(satir, sorgu):
    if sorgu.startswith("rs"):
        return satir["rs"] == sorgu[2:]
    if sorgu.startswith("p."):
        return satir["hgvs_p"] == sorgu[2:]
    return satir["hgvs_c"] == sorgu


tablo = []
for gen, sorgu in RAPOR:
    bulunan = [s for s in cv.variants(gen, max_genes=None) if eslesir(s, sorgu)]
    print(f"{gen} {sorgu}")
    if not bulunan:
        print("  ClinVar'da bulunamadı; yazımı ve transkripti denetleyin\n")
        continue
    for s in bulunan:
        vcf = f"{s['chromosome']}-{s['position_vcf']}-{s['ref']}-{s['alt']}"
        hastalik = [h for h in s["conditions"].split("|")
                    if h and h.lower() not in UNINFORMATIVE_CONDITIONS]
        print(f"  {s['name']}")
        print(f"  {s['significance']}  {'*' * s['stars'] or '0 yıldız'}  ({s['review_status']})")
        print(f"  {s['submitters']} gönderen, son değerlendirme {s['last_evaluated'] or '-'}, "
              f"köken {s['origin']}")
        print(f"  {s['assembly']} {vcf}   rs{s['rs'] or '-'}   VariationID {s['variation_id']}")
        print(f"  Hastalık: {'; '.join(hastalik[:3]) or '-'}\n")
        tablo.append({"rapor_gen": gen, "rapor_sorgu": sorgu, "vcf": vcf, **s})

cv.write_tsv(tablo, CIKTI / "rapor.tsv")
print(f"{len(tablo)} eşleşme yazıldı: {CIKTI / 'rapor.tsv'}")
