"""Grup 6: İki semptom aynı hastalıkta mı birlikte görülüyor?

Senaryo: Nöbet ve mikrosefalinin birlikte görüldüğü genleri arıyorsunuz. İki gen
listesini kesiştirmek yanıltabilir: bir gen, nöbeti bir hastalığında, mikrosefaliyi
başka bir hastalığında taşıyor olabilir. Daha sıkı olan, ikisinin aynı gen-hastalık
çiftinde bulunmasını istemektir.

Bu örnekte:
  1. Liste kesişimi hesaplanır.
  2. annotations ile aynı gen-hastalık çiftindeki birliktelik hesaplanır.
  3. Yalnızca ilk yöntemde kalan genler, yani yalancı birliktelikler listelenir.

Çıktılar: test/hpo/cikti/hpo_grup_6/
Not: Bu iki semptomda fark küçüktür, çünkü genlerin çoğu tek hastalıkla ilişkilidir.
Birden çok hastalığa yol açan genlerde ya da nadir semptomlarda fark büyüyebilir.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_6"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
A, B = "Seizure", "Microcephaly"
KALITE = {"association": "MENDELIAN"}

satirlar_a = hpo.annotations(A, **KALITE)
satirlar_b = hpo.annotations(B, **KALITE)

liste_kesisimi = {s["gene_symbol"] for s in satirlar_a} & {s["gene_symbol"] for s in satirlar_b}

ciftler_a = {(s["gene_symbol"], s["disease_id"]) for s in satirlar_a}
ciftler_b = {(s["gene_symbol"], s["disease_id"]) for s in satirlar_b}
ortak_ciftler = ciftler_a & ciftler_b
ayni_hastalik = {gen for gen, _ in ortak_ciftler}

yalanci = sorted(liste_kesisimi - ayni_hastalik)
print(f"{A}: {len({s['gene_symbol'] for s in satirlar_a})} gen, "
      f"{B}: {len({s['gene_symbol'] for s in satirlar_b})} gen")
print(f"Liste kesişimi: {len(liste_kesisimi)} gen")
print(f"Aynı hastalıkta birlikte: {len(ayni_hastalik)} gen ({len(ortak_ciftler)} gen-hastalık çifti)")
print(f"Yalnızca kesişimde kalan: {len(yalanci)} gen "
      f"(%{100 * len(yalanci) / len(liste_kesisimi):.1f})")

adlar = {s["disease_id"]: s["disease_name"] for s in satirlar_a + satirlar_b}
for gen in yalanci[:3]:
    h_a = sorted({adlar[d] for g, d in ciftler_a if g == gen})
    h_b = sorted({adlar[d] for g, d in ciftler_b if g == gen})
    print(f"\n  {gen}")
    print(f"    {A}: {h_a[0][:60]}")
    print(f"    {B}: {h_b[0][:60]}")

Hpo.write_tsv([{"gene_symbol": g, "disease_id": d, "disease_name": adlar[d]}
               for g, d in sorted(ortak_ciftler)], CIKTI / "ayni_hastalikta.tsv")
(CIKTI / "yalanci_birliktelik.txt").write_text("\n".join(yalanci) + "\n")
