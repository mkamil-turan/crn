"""Grup 10: Uçtan uca çalışma tasarımı; iki hastalık grubu ve kontrol grubu.

Senaryo: Nörodejeneratif ve nörogelişimsel genleri, sinir sistemi bulgusu taşımayan
bir kontrol grubuyla karşılaştıracaksınız. Grupları kurarken üç karar verilir: her
grup nasıl tanımlanacak, ortak genlere ne olacak, kontrol grubu ne kadar büyük olacak.

Bu örnekte:
  1. Gruplar yaş, semptom ve kanıt filtreleriyle tanımlanır.
  2. Kontrol, sinir sistemi notu taşımayan genlerden sabit tohumla rastgele seçilir;
     aynı tohum aynı grubu verir.
  3. Üç grup tek CRN dosyasına yazılır, ortak ve atlanan genler raporlanır.
  4. Dosya yeniden okunup grup, transkript ve HPO bilgileri denetlenir.
  5. overlap kararının grup büyüklüğüne etkisi hesaplanır.
  6. Dosyanın diferansiyel teste açıldığı doğrulanır.

Çıktılar: test/hpo/cikti/hpo_grup_10/
Not: Paket genomik DNA içerir. Kontrol grubu nadir hastalık genlerinden gelir, yani
genomun rastgele bir örneği değildir. Gen uzunluğu, GC oranı ya da ifade düzeyine
göre eşleştirilmemiştir. Sinir sistemi notu bulunmaması, genin sinir sisteminde görevi
olmadığı anlamına gelmez.
"""
import random
import warnings
from collections import Counter
from pathlib import Path

from crn import read_crn_list
from crn.cub import CubWarning, DifferentialCodonUsage
from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_10"
CIKTI.mkdir(parents=True, exist_ok=True)

TOHUM = 42
KONTROL_BOYUTU = 200
hpo = Hpo(verbose=False)
ORTAK = {"association": "MENDELIAN", "evidence": ["PCS", "TAS"]}

adaylar = hpo.genes("Phenotypic abnormality", exclude="Abnormality of the nervous system", **ORTAK)
kimlikler = sorted(satir["ncbi_gene_id"] for satir in adaylar)
kontrol = random.Random(TOHUM).sample(kimlikler, KONTROL_BOYUTU)

GRUPLAR = {
    "dejeneratif": {"term": ["Dementia", "Parkinsonism"], "onset": "Adult onset"},
    "gelisimsel": {"term": "Autistic behavior",
                   "onset": ["Congenital onset", "Infantile onset", "Childhood onset"]},
    "kontrol": {"term": "Phenotypic abnormality", "gene": kontrol},
}

ilk = {ad: len(hpo.genes(**{**ORTAK, **tanim})) for ad, tanim in GRUPLAR.items()}
print(f"Kontrol havuzu: {len(kimlikler)} gen, seçilen {KONTROL_BOYUTU} (tohum {TOHUM})")
print("Filtrelerden geçen, ortak ve atlanan genler çıkmadan önce:", ilk)

paket = CIKTI / "calisma.crn"
with warnings.catch_warnings():
    warnings.simplefilter("ignore")      # ortak genler raporda da yazılıdır
    rapor = hpo.pack(GRUPLAR, paket, **ORTAK)

print(f"\nPaket: {paket.stat().st_size / 1e6:.1f} MB, HPO {rapor['release']}")
print("Son grup büyüklükleri:", rapor["groups"])
print(f"Birden çok grupta olduğu için çıkarılan: {len(rapor['overlap'])} gen  "
      f"{', '.join(rapor['overlap'][:8])}")
print(f"Atlanan: {len(rapor['skipped'])} gen")
for gen, neden in sorted(rapor["skipped"].items())[:4]:
    print(f"  {gen}: {neden}")

genler = read_crn_list(paket)
print("\nDosyadan okunan gruplar:", dict(Counter(g.features.first("group") for g in genler)))
print("Transkript seçimi:", dict(Counter(g.features.first("canonical_tag") for g in genler)))
print("Her genin yapısı 1/1:", all(len(g.transcripts) == len(g.proteins) == 1 for g in genler))
print("Hepsinde HPO notu:", all(g.features.get("hpo_id") for g in genler))

Hpo.write_tsv([{"gen": g.name, "grup": g.features.first("group"),
                "gen_adi": g.features.first("gene_name"),
                "transkript_secimi": g.features.first("canonical_tag"),
                "hastalik_sayisi": len(g.features.get("disease_id", []))} for g in genler],
              CIKTI / "calisma_genleri.tsv")

print("\noverlap kararının etkisi (dejeneratif, gelisimsel; atlanan genler çıkmadan önce):")
dej = {s["ncbi_gene_id"] for s in hpo.genes(**{**ORTAK, **GRUPLAR["dejeneratif"]})}
gel = {s["ncbi_gene_id"] for s in hpo.genes(**{**ORTAK, **GRUPLAR["gelisimsel"]})}
print(f"  exclude: {len(dej - gel)}, {len(gel - dej)}")
print(f"  first:   {len(dej)}, {len(gel - dej)}")

with warnings.catch_warnings():
    warnings.simplefilter("ignore", CubWarning)     # üçüncü grup bu çiftte kullanılmaz
    dcu = DifferentialCodonUsage.from_crn(paket, "dejeneratif", "kontrol")
print(f"\nTeste hazır: {dcu}")
