"""Örnek 10: Türler arası protein karşılaştırması ve tür paketi.

Senaryo: p53 proteininin insan, sıçan, fare ve zebra balığı arasında ne kadar
korunduğuna hızlıca bakmak istiyorsunuz.

Bu örnekte:
  1. Her organizma için ayrı Ncbi nesnesi kurulur.
  2. CDS'si en uzun NM_ transkriptinin proteini alınır.
  3. sliding ile 5 amino asitlik parçalar çıkarılır, paylaşılan oranı hesaplanır.
  4. Proteinler (PRT) tek CRN listesine ve FASTA'ya yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_10/
Not: Paylaşılan k-mer oranı hizalama değildir, kabaca benzerlik gösterir. Yüzde
kimlik için MAFFT ya da BLAST kullanın.
"""
from pathlib import Path

from crn import PRT, read_crn_list, write_crn_list, write_fasta
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_10"
CIKTI.mkdir(parents=True, exist_ok=True)

# Gen sembolleri türe göre değişir.
TURLER = {"human": "TP53", "mouse": "Trp53", "rat": "Tp53", "zebrafish": "tp53"}
K = 5


def ana_protein(organizma: str, sembol: str) -> PRT:
    gen = Ncbi(organism=organizma, verbose=False).gene(sembol)
    adaylar = [r for k, r in gen.transcripts.items() if k.startswith("NM_") and "cds" in r.features]
    rna = max(adaylar, key=lambda r: r.features.first("cds")[1] - r.features.first("cds")[0])
    protein = gen.protein(rna.features.first("protein_id"))
    protein.features.set("organism", organizma)
    protein.name = f"{sembol}_{organizma}"
    return protein


proteinler = [ana_protein(tur, sembol) for tur, sembol in TURLER.items()]
for protein in proteinler:
    print(protein)


def kmerler(protein: PRT) -> set:
    return {pencere.seq for pencere in protein.sliding(K)}


kumeler = [kmerler(p) for p in proteinler]
adlar = [p.features.first("organism") for p in proteinler]

print(f"\nPaylaşılan {K}-mer oranı (Jaccard)")
print(" " * 11 + "".join(f"{ad:>11}" for ad in adlar))
for ad, a in zip(adlar, kumeler):
    print(f"{ad:<11}" + "".join(f"{len(a & b) / len(a | b):>11.2f}" for b in kumeler))

ortak = set.intersection(*kumeler)
print(f"\nDört türde ortak {K}-mer sayısı: {len(ortak)}")
print("Örnekler:", ", ".join(sorted(ortak)[:8]))

yol = CIKTI / "p53_ortologlari.crn"
write_crn_list(proteinler, yol)
write_fasta(proteinler, CIKTI / "p53_ortologlari.fasta")
print("\nPaket:", [x.name for x in read_crn_list(yol)])
print("Çıktılar:", CIKTI)
