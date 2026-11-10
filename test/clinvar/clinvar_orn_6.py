"""Örnek 6: Patojen varyantların protein etkisinden hastalık mekanizmasına.

Senaryo: Bir genin hastalığa işlev kaybıyla mı (haployetmezlik; kesen varyantlar)
yoksa değişmiş protein işleviyle mi (yanlış anlamlı varyantlar, çoğunlukla baskın
negatif ya da işlev kazanımı) yol açtığını ilk bakışta görmek istiyorsunuz. Bu ayrım
tedavi yaklaşımını ve yeni bir varyantın nasıl yorumlanacağını belirler.

Bu örnekte:
  1. Seçili genlerde patojen varyantların etki dağılımı çıkarılır.
  2. Kesen (nonsense + frameshift) ve yanlış anlamlı oranları hesaplanır.
  3. Bütün genomda kesen varyantlarla ve yanlış anlamlı varyantlarla öne çıkan genler
     ayrı ayrı sıralanır.

Çıktılar: test/clinvar/cikti/clinvar_orn_6/mekanizma.tsv
Not: Etki HGVS protein yazımından çıkarılır. Ekleme bölgesi varyantlarının protein
yazımı çoğunlukla yoktur ve bu sayımlara girmez; işlev kaybı oranı bu yüzden alt sınırdır.
"""
from collections import Counter
from pathlib import Path

from crn.fetch import Clinvar

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_6"
CIKTI.mkdir(parents=True, exist_ok=True)

PLP = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 1}
KESEN = ["nonsense", "frameshift"]
cv = Clinvar(verbose=False)

print("Seçili genler (P/LP, en az 1 yıldız):")
print(f"{'gen':<8}{'P/LP':>7}{'kesen%':>8}{'yanlış%':>9}  baskın etki")
for gen in ("BRCA2", "TTN", "NF1", "SCN1A", "KCNQ1", "MYH7", "FGFR3", "RYR1"):
    etki = Counter(r["effect"] for r in cv.variants(gen, **PLP) if r["effect"])
    n = sum(etki.values())
    kesen = etki["nonsense"] + etki["frameshift"]
    print(f"{gen:<8}{n:>7}{100 * kesen / n:>8.0f}{100 * etki['missense'] / n:>9.0f}  "
          f"{etki.most_common(1)[0][0]}")

# Genom ölçeğinde iki sorgu yeter: genes() her gen için süzgeçten geçen sayıyı verir.
kesen = {s["gene_symbol"]: s["variants"] for s in cv.genes(effect=KESEN, **PLP)}
yanlis = {s["gene_symbol"]: s["variants"] for s in cv.genes(effect="missense", **PLP)}
tablo = []
for gen in set(kesen) | set(yanlis):
    k, y = kesen.get(gen, 0), yanlis.get(gen, 0)
    if k + y >= 20:
        tablo.append({"gene": gen, "truncating": k, "missense": y,
                      "truncating_share": round(k / (k + y), 3)})
tablo.sort(key=lambda r: (-r["truncating_share"], -(r["truncating"] + r["missense"])))

print(f"\nEn az 20 P/LP kesen ya da yanlış anlamlı varyantı olan {len(tablo)} gen.")
lof = [r for r in tablo if r["truncating_share"] >= 0.9]
mis = [r for r in tablo if r["truncating_share"] <= 0.1]
print(f"Kesen payı ≥ %90: {len(lof)} gen, ör. "
      + ", ".join(r["gene"] for r in sorted(lof, key=lambda r: -r["truncating"])[:10]))
print(f"Yanlış anlamlı payı ≥ %90: {len(mis)} gen, ör. "
      + ", ".join(r["gene"] for r in sorted(mis, key=lambda r: -r["missense"])[:10]))

cv.write_tsv(tablo, CIKTI / "mekanizma.tsv")
print(f"\nYazıldı: {CIKTI / 'mekanizma.tsv'}  (clinvar_grup_3 bu tabloyu gruplara çevirir)")
