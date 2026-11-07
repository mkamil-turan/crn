"""Örnek 13: Bir gen ailesini karşılaştır: insan akuaporinleri.

Senaryo: Akuaporin ailesinin üyelerini tek pakette toplamak, mRNA bölgelerini
ve aileye özgü NPA motiflerini yan yana görmek istiyorsunuz.

Bu örnekte:
  1. Aile üyeleri indirilir.
  2. Her üyenin ana transkriptinde 5'UTR, CDS, 3'UTR uzunluğu ve GC'si ölçülür.
  3. sliding(3) ile proteinde NPA motifleri aranır.
  4. Genler ve ana proteinler aynı CRN listesine karışık olarak yazılır.
  5. Paket okunur, isinstance ile tiplerine ayrılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_13/
"""
import csv
from pathlib import Path

from crn import PRT, Gene, read_crn_list, write_crn_list
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_13"
CIKTI.mkdir(parents=True, exist_ok=True)

AILE = ["MIP", "AQP1", "AQP2", "AQP3", "AQP4", "AQP5", "AQP6", "AQP7", "AQP8", "AQP9", "AQP10", "AQP11"]


def gc(parca) -> float:
    return 100 * (parca.seq.count("G") + parca.seq.count("C")) / len(parca) if len(parca) else 0.0


ncbi = Ncbi(verbose=False)
genler, proteinler, satirlar = [], [], []
for sembol in AILE:
    gen = ncbi.gene(sembol)
    genler.append(gen)
    adaylar = [r for k, r in gen.transcripts.items() if k.startswith("NM_") and "cds" in r.features]
    rna = max(adaylar, key=lambda r: r.features.first("cds")[1] - r.features.first("cds")[0])
    bas, bit = rna.features.first("cds")

    protein = gen.protein(rna.features.first("protein_id"))[:]
    protein.name = f"{gen.name}_protein"
    proteinler.append(protein)

    npa = [i + 1 for i, p in enumerate(protein.sliding(3)) if p == "NPA"]
    satirlar.append({"gen": gen.name, "transkript": rna.name, "utr5_nt": bas - 1, "cds_nt": bit - bas + 1,
                     "utr3_nt": len(rna) - bit, "cds_gc": round(gc(rna[bas - 1:bit]), 1),
                     "aa": len(protein), "npa": ",".join(map(str, npa)) or "-"})

print(f"{'Gen':<7}{'Transkript':<16}{'5UTR':>6}{'CDS':>6}{'3UTR':>6}{'GC%':>7}{'aa':>5}  NPA")
for s in satirlar:
    print(f"{s['gen']:<7}{s['transkript']:<16}{s['utr5_nt']:>6}{s['cds_nt']:>6}{s['utr3_nt']:>6}"
          f"{s['cds_gc']:>7}{s['aa']:>5}  {s['npa']}")

with open(CIKTI / "aqp_tablosu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)

yol = CIKTI / "aqp_ailesi.crn"
write_crn_list(genler + proteinler, yol)
paket = read_crn_list(yol)
gen_say = sum(isinstance(x, Gene) for x in paket)
prt_say = sum(isinstance(x, PRT) for x in paket)
print(f"\nKarışık paket: {gen_say} Gene + {prt_say} PRT, {yol.stat().st_size / 1e6:.2f} MB")
print("Çıktılar:", CIKTI)
