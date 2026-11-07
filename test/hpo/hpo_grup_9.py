"""Grup 9: Hastanın fenotip listesinden aday gen sıralaması.

Senaryo: Klinisyen bir hastada şu bulguları yazmış: ateşli nöbet, status epileptikus,
miyoklonik nöbet, gelişimsel gerileme ve ataksi. Tek bir semptom yerine bulguların
tümüne en iyi uyan genlerden bir aday grubu kurmak istiyorsunuz.

Bu örnekte:
  1. Her bulgunun bilgi değeri hesaplanır: nadir bulgu daha çok puan getirir.
     puan(terim) = log2(toplam gen / terimle ilişkili gen)
  2. Her gen, sahip olduğu bulguların puanları toplamıyla sıralanır.
  3. En iyi adaylar, eşleşen bulgu ve hastalıklarıyla raporlanır.

Çıktılar: test/hpo/cikti/hpo_grup_9/
Not: Bu, yöntemi gösteren basit bir sıralamadır, tanı aracı değildir. Yalnızca varlığı
sayar; hastada olmayan bulguları cezalandırmaz. Eşleşen bulgular genin farklı
hastalıklarından gelebilir (bkz. Grup 6). Eşit puanlı genler kendi aralarında sıralı
değildir; listede alfabetik görünürler.
"""
import math
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_9"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
KALITE = {"association": "MENDELIAN"}
BULGULAR = ["Febrile seizure (within the age range of 3 months to 6 years)", "Status epilepticus",
            "Myoclonic seizure", "Developmental regression", "Ataxia"]

toplam = len(hpo.genes("Phenotypic abnormality", **KALITE))
setler = {b: set(hpo.symbols(b, **KALITE)) for b in BULGULAR}
puanlar = {b: math.log2(toplam / len(s)) for b, s in setler.items()}

print(f"Toplam {toplam} gen. Bulguların bilgi değeri:")
for bulgu, puan in sorted(puanlar.items(), key=lambda x: -x[1]):
    print(f"  {puan:>5.2f}  {len(setler[bulgu]):>4} gen  {bulgu[:50]}")

adaylar = []
for gen in set().union(*setler.values()):
    eslesen = [b for b in BULGULAR if gen in setler[b]]
    adaylar.append({"gene_symbol": gen, "puan": round(sum(puanlar[b] for b in eslesen), 3),
                    "eslesen_bulgu": len(eslesen), "bulgular": "|".join(eslesen)})
adaylar.sort(key=lambda s: (-s["puan"], s["gene_symbol"]))
Hpo.write_tsv(adaylar, CIKTI / "aday_genler.tsv")

print(f"\nEn iyi adaylar ({len(adaylar)} genden):")
for aday in adaylar[:8]:
    hastaliklar = hpo.diseases(gene=aday["gene_symbol"], term=BULGULAR, **KALITE)
    print(f"  {aday['gene_symbol']:<9}{aday['puan']:>6.2f}  {aday['eslesen_bulgu']}/5  "
          f"{hastaliklar[0]['disease_name'][:50]}")

scn1a = next(a for a in adaylar if a["gene_symbol"] == "SCN1A")
onde = sum(a["puan"] > scn1a["puan"] for a in adaylar)
esit = sum(a["puan"] == scn1a["puan"] for a in adaylar)
print(f"\nDravet sendromunun geni SCN1A: önünde {onde} gen var, {esit} gen aynı puanda")
