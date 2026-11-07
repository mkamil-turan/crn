"""Örnek 5: Semptomların gen örtüşmesi.

Senaryo: Epilepsi ile otizm klinikte sık birlikte görülür. İki semptomun gen
setleri ne kadar örtüşüyor? Birine özgü genler hangileri? Bu soru, gen gruplarını
karşılaştırmadan önce sorulmalıdır, çünkü ortak genler grup farkını gizler.

Bu örnekte:
  1. Üç nörogelişimsel semptomun gen setleri aynı filtreyle çıkarılır.
  2. İkili örtüşmeler ve Jaccard benzerliği yazılır.
  3. Nöbete özgü genler ``exclude`` ile bulunur.

Çıktılar: test/hpo/cikti/hpo_orn_5/
Not: Aynı hastalık birçok semptomu taşıdığı için örtüşme yüksek çıkar. Örtüşme,
ortak mekanizmayı değil, HPO'daki ortak hastalık notlarını gösterir.
"""
from itertools import combinations
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_orn_5"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
FILTRE = {"association": "MENDELIAN", "evidence": ["PCS", "TAS"]}
SEMPTOMLAR = {"Nöbet": "Seizure", "Otizm": "Autistic behavior", "Ataksi": "Ataxia"}

setler = {ad: set(hpo.symbols(terim, **FILTRE)) for ad, terim in SEMPTOMLAR.items()}
for ad, genler in setler.items():
    print(f"{ad:<8}{len(genler):>6} gen")

print("\nİkili örtüşme:")
for a, b in combinations(setler, 2):
    ortak = setler[a] & setler[b]
    birlesim = setler[a] | setler[b]
    print(f"  {a} ve {b}: {len(ortak)} ortak gen, Jaccard {len(ortak) / len(birlesim):.2f}")

ozgu = hpo.genes("Seizure", exclude=["Autistic behavior", "Ataxia"], **FILTRE)
Hpo.write_tsv(ozgu, CIKTI / "nobete_ozgu.tsv")
print(f"\nOtizm ve ataksi notu olmayan nöbet genleri: {len(ozgu)}")
print("  En çok hastalıkla ilişkili:", ", ".join(s["gene_symbol"] for s in ozgu[:8]))

uclu = set.intersection(*setler.values())
print(f"\nÜç semptomda ortak: {len(uclu)} gen")
