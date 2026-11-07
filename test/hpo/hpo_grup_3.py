"""Grup 3: Birkaç semptomun birleşiminden tek grup.

Senaryo: "Hareket bozuklukları" tek bir terim değil; distoni, kore ve parkinsonizm
ayrı dallardadır. Üçünden en az birine sahip genleri tek grupta toplamak, ama her
genin kaç semptomla ilişkili olduğunu da görmek istiyorsunuz.

Bu örnekte:
  1. Terim listesi verilir; sonuç birleşimdir (VEYA).
  2. Her genin hangi semptomlarda geçtiği tablolanır.
  3. Üç semptomda birden geçen genler ayrı yazılır.

Çıktılar: test/hpo/cikti/hpo_grup_3/
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_3"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
SEMPTOMLAR = ["Dystonia", "Chorea", "Parkinsonism"]
KALITE = {"association": "MENDELIAN"}

birlesim = hpo.symbols(SEMPTOMLAR, **KALITE)
setler = {terim: set(hpo.symbols(terim, **KALITE)) for terim in SEMPTOMLAR}

print(f"Hareket bozukluğu grubu: {len(birlesim)} gen")
for terim, genler in setler.items():
    print(f"  {terim:<14}{len(genler):>5} gen")

satirlar = [{"gene_symbol": gen, **{t: int(gen in setler[t]) for t in SEMPTOMLAR},
             "toplam": sum(gen in setler[t] for t in SEMPTOMLAR)} for gen in birlesim]
satirlar.sort(key=lambda s: (-s["toplam"], s["gene_symbol"]))
Hpo.write_tsv(satirlar, CIKTI / "hareket_bozukluklari.tsv")

ucu_de = [s["gene_symbol"] for s in satirlar if s["toplam"] == 3]
print(f"\nÜç semptomda birden: {len(ucu_de)} gen")
print(" ", ", ".join(ucu_de[:15]))
