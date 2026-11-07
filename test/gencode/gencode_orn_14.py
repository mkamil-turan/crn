"""Örnek 14: Mitokondri genomunun tam envanteri ve doğru genetik kod.

Senaryo: chrM'nin tüm genlerini görmek, 13 proteini mitokondriyal kodla çevirip
GENCODE proteinleriyle karşılaştırmak istiyorsunuz.

Bu örnekte:
  1. chrM genleri tipleriyle listelenir.
  2. Protein genleri iki kodla çevrilip iç stoplar sayılır.
  3. Çeviri GENCODE proteiniyle karşılaştırılır.
  4. Örtüşen genler (ATP8/ATP6, ND4L/ND4) bulunur.
  5. Proteinler (PRT) tek CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_14/
Not: Bazı mitokondri mRNA'larında stop kodonu poliadenilasyonla tamamlanır; ayrıca
AUU, AUA gibi başlangıçlar yine metiyoninle başlar.
"""
from collections import Counter
from pathlib import Path

from crn import GeneticCode, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_14"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)
genler = sorted(gc.genes(chromosome="chrM"), key=lambda g: g["start"])
print(f"chrM: {len(genler)} gen:", dict(Counter(g["gene_type"] for g in genler)))

print(f"\n{'Gen':<10}{'Konum':>13}{'Z':>3}{'aa':>6}{'Baş':>6}{'Son 3 nt':>9}{'Std stop':>10}{'Mt stop':>9}  GENCODE ile")
proteinler = []
for g in genler:
    if g["gene_type"] != "protein_coding":
        continue
    rna = gc.rna(gc.transcripts(gene=g["gene_id"])[0]["transcript_id"])
    bas, bit = rna.features.first("cds", (1, len(rna)))
    cds = rna[bas - 1:bit]
    standart = cds.translation(GeneticCode.STANDARD).seq
    mito = cds.translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL)
    referans = gc.protein(rna.name)

    # İlk kalıntı başlangıç kodonu olduğu için her durumda M okunur.
    govde = mito.seq.rstrip("*")[1:]
    durum = "aynı" if govde == referans.seq[1:] else "FARKLI"
    mito.name = g["gene_name"].replace("-", "_")
    proteinler.append(mito)
    konum = f"{g['start']}-{g['end']}"
    print(f"{g['gene_name']:<10}{konum:>13}{g['strand']:>3}{len(referans):>6}{cds.seq[:3]:>6}{cds.seq[-3:]:>9}"
          f"{standart[:-1].count('*'):>10}{mito.seq[:-1].count('*'):>9}  {durum}")

print("\nÖrtüşen protein genleri:")
kodlayan = [g for g in genler if g["gene_type"] == "protein_coding" and g["strand"] == "+"]
for a, b in zip(kodlayan, kodlayan[1:]):
    if a["end"] >= b["start"]:
        print(f"  {a['gene_name']} / {b['gene_name']}: {a['end'] - b['start'] + 1} nt")

write_crn_list(proteinler, CIKTI / "mitokondri_proteinleri.crn")
print("\nÇıktılar:", CIKTI)
