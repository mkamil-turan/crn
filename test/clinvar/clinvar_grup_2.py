"""Grup 2: Aynı hastalık için HPO ve ClinVar gen listelerini karşılaştırmak.

Senaryo: Hipertrofik kardiyomiyopati gen grubunu kuracaksınız. HPO'da genler
fenotip notlarından ("bu gendeki hastalıkta bu bulgu görülür"), ClinVar'da patojen
varyant bildirimlerinden gelir. İki kaynağın ortak genleri en güçlü adaylardır;
yalnızca birinde olan genler ayrı ayrı incelenmelidir.

Bu örnekte:
  1. HPO'dan "Hypertrophic cardiomyopathy" terimiyle Mendel genleri alınır.
  2. ClinVar'dan aynı hastalığın en az 1 yıldızlı P/LP varyantı olan genleri alınır.
  3. Kesişim ve farklar hesaplanır; yalnızca HPO'da olanların ClinVar kanıtı gösterilir.
  4. Her gen için iki kaynağın bilgisi tek TSV'de birleştirilir.

Çıktılar: test/clinvar/cikti/clinvar_grup_2/hcm_kaynaklar.tsv
Not: HPO'da HCM, bir sendromun (ör. Noonan, Pompe) bulgularından biri olarak da
geçer; bu genler ClinVar'da çoğunlukla başka hastalık adıyla kayıtlıdır. Fark her
zaman eksik veri anlamına gelmez, sorunun nasıl sorulduğunu gösterir.
"""
from pathlib import Path

from crn.fetch import Clinvar, Hpo

CIKTI = Path(__file__).parent / "cikti" / "clinvar_grup_2"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
cv = Clinvar(verbose=False)
PLP = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 1}

hpo_genler = {g["gene_symbol"]: g for g in hpo.genes("Hypertrophic cardiomyopathy", association="MENDELIAN")}
cv_genler = {g["gene_symbol"]: g for g in cv.genes("hypertrophic cardiomyopathy", **PLP)}

ortak = sorted(set(hpo_genler) & set(cv_genler))
yalniz_hpo = sorted(set(hpo_genler) - set(cv_genler))
yalniz_cv = sorted(set(cv_genler) - set(hpo_genler))
print(f"HPO {hpo.release}: {len(hpo_genler)} gen,  ClinVar {cv.release}: {len(cv_genler)} gen")
print(f"Ortak: {len(ortak)}  yalnız HPO: {len(yalniz_hpo)}  yalnız ClinVar: {len(yalniz_cv)}\n")

print("Ortak genler, ClinVar P/LP sayısına göre:")
print("  " + ", ".join(f"{g}({cv_genler[g]['pathogenic']})"
                       for g in sorted(ortak, key=lambda g: -cv_genler[g]["pathogenic"])[:15]))

print("\nYalnız ClinVar'da (HCM adıyla P/LP bildirimi var, HPO notu yok):")
print("  " + ", ".join(f"{g}({cv_genler[g]['pathogenic']})"
                       for g in sorted(yalniz_cv, key=lambda g: -cv_genler[g]["pathogenic"])[:15]))

print("\nYalnız HPO'da olanların ClinVar'daki en sık hastalık adı (ilk 8):")
for gen in yalniz_hpo[:8]:
    try:
        ad = cv.conditions(gene=gen, limit=1)
    except KeyError:
        ad = []
    print(f"  {gen:<9} {ad[0]['condition_name'] if ad else '(ClinVar kaydı yok)'}")

tablo = [{"gene": g, "in_hpo": g in hpo_genler, "in_clinvar": g in cv_genler,
          "clinvar_plp": cv_genler[g]["pathogenic"] if g in cv_genler else 0,
          "clinvar_max_stars": cv_genler[g]["max_stars"] if g in cv_genler else "",
          "hpo_disease_count": hpo_genler[g]["diseases"] if g in hpo_genler else 0}
         for g in sorted(set(hpo_genler) | set(cv_genler))]
cv.write_tsv(tablo, CIKTI / "hcm_kaynaklar.tsv")
print(f"\nYazıldı: {CIKTI / 'hcm_kaynaklar.tsv'}")
