"""Grup 1: Tek semptomdan gen listesi (en basit senaryo).

Senaryo: Parkinsonizm ile ilişkili genlerin listesini bir arkadaşınıza ya da
başka bir programa vereceksiniz. Tek gereken semptomun adı.

Bu örnekte:
  1. Terim adıyla sorgu yapılır; alt terimler kendiliğinden dahildir.
  2. Semboller alfabetik liste olarak metin dosyasına yazılır.

Çıktılar: test/hpo/cikti/hpo_grup_1/
Not: Bu ham listedir, hiçbir güven filtresi uygulanmamıştır. Grup 2 bu listeyi daraltır.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_1"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
genler = hpo.symbols("Parkinsonism")

(CIKTI / "parkinsonizm.txt").write_text("\n".join(genler) + "\n")

print(f"Parkinsonizm: {len(genler)} gen")
print("İlk 10:", ", ".join(genler[:10]))
print("Bilinen genler listede mi:")
for gen in ("SNCA", "LRRK2", "PRKN", "PINK1", "PARK7"):
    print(f"  {gen:<6} {'var' if gen in genler else 'YOK'}")
