"""Örnek 6: Semptoma göre gen gruplarını tek CRN dosyasına paketlemek.

Senaryo: Daha önce DisGeNET ile yaptığınız "semptomdan gen grubu" adımını HPO ile
kuruyorsunuz. İki semptom grubunu tek dosyada toplayıp diferansiyel kodon testine
doğrudan açmak istiyorsunuz.

Bu örnekte:
  1. Her grup bir terim ve isteğe bağlı kendi filtreleriyle tanımlanır.
  2. İki grupta birden olan genler çıkarılır (overlap="exclude").
  3. Grup adı her genin features'ına yazılır, liste tek CRN dosyası olur.
  4. Rapor okunur: grup büyüklüğü, ortak ve atlanan genler.
  5. Dosya DifferentialCodonUsage.from_crn ile açılır.

Çıktılar: test/hpo/cikti/hpo_orn_6/
Not: Gruplar büyüklük olarak dengeli değilse sonuç büyük gruba göre okunmalıdır.
Paket genomik DNA içerdiği için yüzlerce MB olabilir.
"""
import warnings
from collections import Counter
from pathlib import Path

from crn import read_crn_list
from crn.cub import DifferentialCodonUsage
from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_orn_6"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo()
GRUPLAR = {
    "nobet": {"term": "Seizure", "min_frequency": 0.8},
    "ataksi": {"term": "Ataxia", "min_frequency": 0.8},
}
paket = CIKTI / "gruplar.crn"

with warnings.catch_warnings():
    warnings.simplefilter("ignore")      # ortak genler raporda da yazılıdır
    rapor = hpo.pack(GRUPLAR, paket, association="MENDELIAN")

print(f"\nHPO {rapor['release']}, dosya {paket.stat().st_size / 1e6:.1f} MB")
print("Grup büyüklükleri:", rapor["groups"])
print(f"İki grupta birden olduğu için çıkarılan: {len(rapor['overlap'])} gen")
print(f"Atlanan: {len(rapor['skipped'])} gen")
for gen, neden in sorted(rapor["skipped"].items())[:5]:
    print(f"  {gen}: {neden}")

genler = read_crn_list(paket)
print("\nDosyadaki gruplar:", dict(Counter(g.features.first("group") for g in genler)))

dcu = DifferentialCodonUsage.from_crn(paket, "nobet", "ataksi")
print(dcu)
print(f"Aile düzeyi Fisher: {len(dcu.significant('family'))} kodon, "
      f"permütasyondan da geçen: {len(dcu.robust('family'))} kodon")
