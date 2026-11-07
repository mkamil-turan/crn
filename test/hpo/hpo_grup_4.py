"""Grup 4: Dışlama ile daha saf bir grup.

Senaryo: Spastisite genlerini arıyorsunuz, ama zihinsel yetersizlik ve nöbet eşlik
eden yaygın nörogelişimsel sendromları istemiyorsunuz. Amacınız motor sisteme daha
sınırlı görünen bir grup.

Bu örnekte:
  1. Dışlanacak semptomlar tek tek eklenir, grubun nasıl küçüldüğü görülür.
  2. Kalan grup ve dışlananlar ayrı dosyalara yazılır.

Çıktılar: test/hpo/cikti/hpo_grup_4/
Not: exclude geni, dışlanan semptom o genin başka bir hastalığında geçse bile
çıkarır ve filtrelerden etkilenmez. Sonuç "bu notu hiç taşımayan genler"dir; o
semptomun gerçekten yokluğunu kanıtlamaz, nadir hastalıklarda notlar eksik olabilir.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_4"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
TERIM = "Spasticity"
KALITE = {"association": "MENDELIAN"}

tum = set(hpo.symbols(TERIM, **KALITE))
print(f"Spastisite: {len(tum)} gen")

dislanan = []
for semptom in ("Intellectual disability", "Seizure", "Global developmental delay"):
    dislanan.append(semptom)
    kalan = hpo.symbols(TERIM, exclude=dislanan, **KALITE)
    print(f"  - {semptom:<28}: {len(kalan):>4} gen kaldı")

saf = hpo.genes(TERIM, exclude=dislanan, **KALITE)
semboller = {s["gene_symbol"] for s in saf}
Hpo.write_tsv(saf, CIKTI / "saf_spastisite.tsv")
(CIKTI / "dislananlar.txt").write_text("\n".join(sorted(tum - semboller)) + "\n")

print(f"\nSaf grup: {len(saf)} gen")
print("Örnekler:", ", ".join(s["gene_symbol"] for s in saf[:10]))
