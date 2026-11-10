"""Örnek 5: Bir tanı panelinde anlamı belirsiz (VUS) ve çelişkili varyant yükü.

Senaryo: Kalıtsal kardiyomiyopati ve aritmi paneli raporlayan bir laboratuvardasınız.
Hangi genlerde raporlanan varyantların çoğunun belirsiz kaldığını, hangilerinde
laboratuvarlar arası çelişkinin yüksek olduğunu bilmek, yorum ve yeniden sınıflandırma
iş yükünü planlamanızı sağlar.

Bu örnekte:
  1. Panel genleri için gen başına sınıf sayıları alınır.
  2. VUS oranı, çelişki oranı ve uzman paneli (3+ yıldız) desteği hesaplanır.
  3. Genler VUS oranına göre sıralanır ve TSV olarak yazılır.

Çıktılar: test/clinvar/cikti/clinvar_orn_5/panel.tsv
Not: Oranlar ClinVar'a gönderilen varyantlar üzerindendir, popülasyondaki taşıyıcı
sıklığını göstermez. Büyük genler (TTN) yalnızca uzunlukları nedeniyle çok VUS toplar.
"""
from pathlib import Path

from crn.fetch import Clinvar

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_5"
CIKTI.mkdir(parents=True, exist_ok=True)

PANEL = ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL2", "MYL3",
         "LMNA", "TTN", "PKP2", "DSP", "DSG2", "KCNQ1", "KCNH2", "SCN5A", "RYR2", "CASQ2"]
cv = Clinvar(verbose=False)

genler = {s["gene_symbol"]: s for s in cv.genes(gene=PANEL)}
uzman = {s["gene_symbol"]: s["variants"] for s in cv.genes(gene=PANEL, min_stars=3)}
plp2 = {s["gene_symbol"]: s["variants"]
        for s in cv.genes(gene=PANEL, significance=["pathogenic", "likely_pathogenic"], min_stars=2)}

tablo = []
for gen in PANEL:
    s = genler[gen]
    tablo.append({
        "gene": gen, "variants": s["variants"], "pathogenic": s["pathogenic"],
        "pathogenic_2star": plp2.get(gen, 0), "uncertain": s["uncertain"],
        "conflicting": s["conflicting"], "benign": s["benign"],
        "vus_pct": round(100 * s["uncertain"] / s["variants"], 1),
        "conflict_pct": round(100 * s["conflicting"] / s["variants"], 1),
        "expert_panel": uzman.get(gen, 0),
    })
tablo.sort(key=lambda r: -r["vus_pct"])

print(f"{'gen':<8}{'varyant':>8}{'P/LP':>7}{'P/LP≥2*':>9}{'VUS%':>7}{'çelişki%':>10}{'3+*':>6}")
for r in tablo:
    print(f"{r['gene']:<8}{r['variants']:>8}{r['pathogenic']:>7}{r['pathogenic_2star']:>9}"
          f"{r['vus_pct']:>7}{r['conflict_pct']:>10}{r['expert_panel']:>6}")

toplam = sum(r["variants"] for r in tablo)
vus = sum(r["uncertain"] for r in tablo)
print(f"\nPanel toplamı: {toplam:,} varyant, %{100 * vus / toplam:.1f} VUS")
uzmanli = [r["gene"] for r in tablo if r["expert_panel"]]
print("Uzman paneli (ClinGen VCEP) değerlendirmesi olan genler:", ", ".join(uzmanli) or "yok")

cv.write_tsv(tablo, CIKTI / "panel.tsv")
print(f"Yazıldı: {CIKTI / 'panel.tsv'}")
