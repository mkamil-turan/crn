"""Örnek 4: Gen panelini indir ve tek CRN dosyasında paketle.

Senaryo: Kalıtsal kanser paneli için genleri tek tek dosyalarda tutmak yerine
tek bir paket halinde saklamak ve laboratuvardaki başka bir bilgisayara taşımak
istiyorsunuz.

Bu örnekte:
  1. Gen listesi tanımlanır ve Gene nesneleri indirilir.
  2. API anahtarı varsa ortam değişkeninden alınır (koda yazılmaz).
  3. Liste write_crn_list ile tek dosyaya yazılır.
  4. Dosya read_crn_list ile okunup özgün listeyle karşılaştırılır.
  5. Panel özeti TSV olarak kaydedilir.

Çıktılar: test/ncbi/cikti/ncbi_orn_4/
Sonraki örnekler (5, 6, 8) bu paketi kullanır.
"""
import csv
import os
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Ncbi, NcbiSession

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_4"
CIKTI.mkdir(parents=True, exist_ok=True)

PANEL = ["BRCA1", "BRCA2", "TP53", "PTEN", "ATM", "CHEK2", "PALB2"]

with NcbiSession(api_key=os.environ.get("NCBI_API_KEY")) as oturum:
    ncbi = Ncbi(session=oturum, verbose=False)
    genler = [ncbi.gene(sembol) for sembol in PANEL]

yol = CIKTI / "kanser_paneli.crn"
write_crn_list(genler, yol)
paket = read_crn_list(yol)

print(f"{'Gen':<8}{'Konum':<12}{'Zincir':<8}{'DNA nt':>9}{'NM':>5}{'XM':>6}{'Protein':>9}")
satirlar = []
for gen in paket:
    bolge = gen.features.first("genomic_region")
    nm = sum(k.startswith("NM_") for k in gen.transcripts)
    xm = sum(k.startswith("XM_") for k in gen.transcripts)
    satir = {"gen": gen.name, "konum": gen.features.first("map_location"), "zincir": bolge[3],
             "dna_nt": len(gen.dna), "nm": nm, "xm": xm, "protein": len(gen.proteins)}
    satirlar.append(satir)
    print(f"{gen.name:<8}{satir['konum']:<12}{bolge[3]:<8}{len(gen.dna):>9,}{nm:>5}{xm:>6}{len(gen.proteins):>9}")

ayni = all(a.name == b.name and a.dna == b.dna and list(a.transcripts) == list(b.transcripts)
           for a, b in zip(genler, paket)) and len(genler) == len(paket)
print(f"\nPaket: {len(paket)} gen, {yol.stat().st_size / 1e6:.1f} MB")
print("Okunan paket indirilenle aynı mı:", ayni)

with open(CIKTI / "panel_ozeti.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)
print("Çıktılar:", CIKTI)
