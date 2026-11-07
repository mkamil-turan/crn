"""Örnek 2: Bir transkriptin anatomisi: 5'UTR, CDS, 3'UTR ve protein.

Senaryo: MEFV'nin NM_000243 transkriptinde hangi bölgenin proteini kodladığını
görmek ve NCBI'ın proteiniyle kendi çevirinizi karşılaştırmak istiyorsunuz.

Bu örnekte:
  1. rna() ile transkript alınır; sürüm yazmak gerekmez.
  2. features içindeki koordinatlarla bölgeler dilimlenir.
  3. CDS proteine çevrilir, protein() sonucuyla karşılaştırılır.
  4. Bölgeler FASTA'ya ve bir CRN listesine yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_2/
"""
from pathlib import Path

from crn import read_crn_list, write_crn_list, write_fasta
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_2"
CIKTI.mkdir(parents=True, exist_ok=True)

ncbi = Ncbi(verbose=False)
rna = ncbi.rna("NM_000243")
print(rna)
print("Başlık:", rna.features.first("title"))

# Koordinatlar 1 tabanlı ve iki ucu dahildir; Python dilimi için baş - 1 alınır.
bolgeler = []
for anahtar, etiket in (("utr5", "5UTR"), ("cds", "CDS"), ("utr3", "3UTR")):
    bas, bit = rna.features.first(anahtar)
    parca = rna[bas - 1:bit]
    parca.name = f"{rna.name}_{etiket}"
    bolgeler.append(parca)
    gc = 100 * (parca.seq.count("G") + parca.seq.count("C")) / len(parca)
    print(f"{etiket:<5} {bas:>5}-{bit:<5} {len(parca):>5} nt  GC %{gc:.1f}")

cds = bolgeler[1]
print("\nİlk kodon:", cds.seq[:3], "| son kodon:", cds.seq[-3:], "| 3'e bölünür mü:", len(cds) % 3 == 0)

ceviri = cds.translation()
ceviri.name = "MEFV_ceviri"
ncbi_protein = ncbi.protein("NM_000243")
print("NCBI proteini:", ncbi_protein)
print("Çeviri NCBI ile aynı mı:", ceviri.seq.rstrip("*") == ncbi_protein.seq)

write_fasta(bolgeler + [ncbi_protein], CIKTI / "mefv_bolgeleri.fasta")
yol = CIKTI / "mefv_bolgeleri.crn"
write_crn_list(bolgeler + [ncbi_protein], yol)
print("\nCRN listesinden:", [f"{x.name} ({type(x).__name__})" for x in read_crn_list(yol)])
print("Çıktılar:", CIKTI)
