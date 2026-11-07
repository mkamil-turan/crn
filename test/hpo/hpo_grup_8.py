"""Grup 8: Kalıtım biçimine göre çok gruplu ayrım.

Senaryo: Epilepsi genlerinde kalıtım biçimi hastalık mekanizmasını düşündürür:
dominant formlarda tek kopya yetersizliği ya da işlev kazancı, resesif formlarda
işlev kaybı daha çok beklenir. Aynı semptomu kalıtıma göre bölmek istiyorsunuz.

Bu örnekte:
  1. Aynı semptom üç kalıtım grubuna bölünür.
  2. Her kalıtım grubunun her semptomdaki gen sayısı bir tabloya yazılır.
  3. Hem dominant hem resesif geçen genler ayrılır.

Çıktılar: test/hpo/cikti/hpo_grup_8/
Not: Filtre hastalık düzeyindedir. Bir gen dominant ve resesif iki ayrı hastalıkla
ilişkiliyse ya da tek hastalıkta iki kalıtım da yazılıysa iki gruba da girer.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_8"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
KALITIM = {"AD": "Autosomal dominant inheritance", "AR": "Autosomal recessive inheritance",
           "XL": "X-linked inheritance"}
SEMPTOMLAR = ["Seizure", "Hearing impairment", "Cardiomyopathy"]
KALITE = {"association": "MENDELIAN"}

print(f"{'Semptom':<22}" + "".join(f"{k:>7}" for k in KALITIM) + f"{'AD+AR':>8}")
tablo = []
for semptom in SEMPTOMLAR:
    setler = {k: set(hpo.symbols(semptom, inheritance=v, **KALITE)) for k, v in KALITIM.items()}
    ortak = setler["AD"] & setler["AR"]
    print(f"{semptom:<22}" + "".join(f"{len(s):>7}" for s in setler.values()) + f"{len(ortak):>8}")
    tablo.append({"semptom": semptom, **{k: len(s) for k, s in setler.items()}, "AD+AR": len(ortak)})
Hpo.write_tsv(tablo, CIKTI / "kalitim_tablosu.tsv")

nobet = {k: set(hpo.symbols("Seizure", inheritance=v, **KALITE)) for k, v in KALITIM.items()}
for ad, genler in nobet.items():
    (CIKTI / f"nobet_{ad}.txt").write_text("\n".join(sorted(genler)) + "\n")

ortak = sorted(nobet["AD"] & nobet["AR"])
print(f"\nNöbette hem dominant hem resesif: {len(ortak)} gen")
print(" ", ", ".join(ortak[:12]))
