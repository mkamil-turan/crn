"""Örnek 15: Mitokondri genini genomdan çek ve doğru kodla çevir.

Senaryo: MT-CO1 barkod bölgesi çalışıyorsunuz. Geni rCRS mitokondri genomundan
koordinatla almak ve doğru genetik kodla proteine çevirmek gerekiyor.

Bu örnekte:
  1. region() ile MT-CO1 bölgesi NC_012920.1'den alınır.
  2. Dizi standart ve mitokondriyal kodla çevrilip stop sayıları karşılaştırılır.
  3. Sonuç NCBI'ın MT-CO1 proteiniyle (YP_003024028) doğrulanır.
  4. DNA, RNA ve protein tek CRN listesine yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_15/
"""
from pathlib import Path

from crn import GeneticCode, read_crn_list, write_crn_list
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_15"
CIKTI.mkdir(parents=True, exist_ok=True)

ncbi = Ncbi(verbose=False)
dna = ncbi.region("NC_012920.1", 5904, 7445, "+", name="MT_CO1")
rna = dna.transcription(template=False)
rna.name = "MT_CO1_mRNA"
print(dna)
print("Uzunluk 3'e bölünüyor mu:", len(dna) % 3 == 0)

standart = rna.translation(GeneticCode.STANDARD)
mito = rna.translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL)
mito.name = "MT_CO1_protein"
print(f"\nStandart kod      : {standart.seq[:40]}... iç stop {standart.seq[:-1].count('*')}")
print(f"Mitokondriyal kod : {mito.seq[:40]}... iç stop {mito.seq[:-1].count('*')}")

ncbi_protein = ncbi.protein("YP_003024028")
print("\nNCBI proteini:", ncbi_protein)
print("Mitokondriyal çeviri NCBI ile aynı mı:", mito.seq.rstrip("*") == ncbi_protein.seq)
print("Standart çeviri NCBI ile aynı mı     :", standart.seq.rstrip("*") == ncbi_protein.seq)

yol = CIKTI / "mt_co1.crn"
write_crn_list([dna, rna, mito], yol)
print("\nPaket:", [f"{x.name} ({type(x).__name__})" for x in read_crn_list(yol)])
print("Çıktılar:", CIKTI)
