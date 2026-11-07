"""Grup 2: Geniş listeden güvenilir çekirdek grup.

Senaryo: Kardiyomiyopati genleriyle çalışacaksınız. Ham liste, kardiyomiyopatinin
nadiren görüldüğü sendromların genlerini de içerir. Semptomun hastaların çoğunda
görüldüğü, Mendel tipi ve yayınla desteklenen genleri ayırmak istiyorsunuz.

Bu örnekte:
  1. Aynı terimden geniş ve çekirdek iki grup çıkarılır.
  2. Eşik değiştikçe grup büyüklüğü izlenir.
  3. Çekirdek grup gen başına özetle TSV'ye yazılır.

Çıktılar: test/hpo/cikti/hpo_grup_2/
Not: min_frequency sıklığı bilinmeyen satırları eler, yani eşik yükseldikçe yalnızca
zayıf genler değil, sıklığı raporlanmamış genler de düşer.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_2"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
TERIM = "Cardiomyopathy"
KALITE = {"association": "MENDELIAN", "evidence": ["PCS", "TAS"]}

genis = set(hpo.symbols(TERIM))
print(f"Geniş grup: {len(genis)} gen\n")

print("Sıklık eşiğine göre (Mendel, yayın/uzman kanıtı):")
for esik in (0.05, 0.3, 0.5, 0.8, 1.0):
    print(f"  en az %{esik * 100:>3.0f}: {len(hpo.symbols(TERIM, min_frequency=esik, **KALITE)):>4} gen")

cekirdek = hpo.genes(TERIM, min_frequency=0.8, **KALITE)
Hpo.write_tsv(cekirdek, CIKTI / "kardiyomiyopati_cekirdek.tsv")

print(f"\nÇekirdek grup (en az %80): {len(cekirdek)} gen")
print("En çok hastalıkla ilişkili:", ", ".join(s["gene_symbol"] for s in cekirdek[:10]))
