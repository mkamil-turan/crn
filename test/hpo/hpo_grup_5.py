"""Grup 5: Semptom yerine hastalık ailesinden grup.

Senaryo: Bir semptomla değil, adı belli bir hastalık ailesiyle çalışıyorsunuz:
spinoserebellar ataksiler, Charcot-Marie-Tooth, retinitis pigmentosa. Bu ailelerde
onlarca alt tip vardır ve her biri ayrı OMIM kaydıdır.

Bu örnekte:
  1. Hastalık adı düzenli ifadeyle eşleştirilir, kaç alt tip yakalandığı yazılır.
  2. Her aile için gen listesi çıkarılır.
  3. Aileler arası ortak genler bulunur.

Çıktılar: test/hpo/cikti/hpo_grup_5/
Not: Ad araması yalnızca adı içeren kayıtları bulur. Düzenli ifadeyi seçtikten sonra
hastalık listesini gözden geçirin; gereksiz kayıtlar araya girebilir.
"""
import re
from itertools import combinations
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_5"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
AILELER = {
    "SCA": re.compile(r"spinocerebellar ataxia", re.I),
    "CMT": re.compile(r"charcot-marie-tooth", re.I),
    "RP": re.compile(r"retinitis pigmentosa", re.I),
}

setler = {}
for ad, desen in AILELER.items():
    hastaliklar = hpo.diseases(desen)
    genler = hpo.symbols(disease=desen)
    setler[ad] = set(genler)
    (CIKTI / f"{ad}.txt").write_text("\n".join(genler) + "\n")
    Hpo.write_tsv(hastaliklar, CIKTI / f"{ad}_hastaliklar.tsv")
    print(f"{ad}: {len(hastaliklar)} hastalık kaydı, {len(genler)} gen")
    print(f"  Örnek kayıt: {hastaliklar[0]['disease_name']}")

print("\nAileler arası ortak genler:")
for a, b in combinations(setler, 2):
    ortak = sorted(setler[a] & setler[b])
    print(f"  {a} ve {b}: {len(ortak)}  {', '.join(ortak[:8])}")
