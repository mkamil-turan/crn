"""Örnek 1: Bir geni NCBI'dan indir, incele ve kaydet.

Senaryo: Ailesel Akdeniz Ateşi çalışıyorsunuz ve MEFV geninin genomik dizisini,
transkriptlerini ve proteinlerini tek nesnede görmek istiyorsunuz.

Bu örnekte:
  1. Ncbi nesnesi kurulur, gen sembolü Gene ID'ye çevrilir.
  2. gene() ile DNA, tüm RefSeq RNA'ları ve proteinler alınır.
  3. Gen bilgileri ve dizileri ekrana yazılır.
  4. İkinci çağrının yerel depodan (sys/) çok daha hızlı geldiği ölçülür.
  5. Gen tek CRN dosyasına yazılır, ağ olmadan geri okunur.

Çıktılar: test/ncbi/cikti/ncbi_orn_1/
Not: İlk çalıştırma internet ister. Kayıtlar sonra sys/ncbi/ içinden okunur.
"""
import time
from pathlib import Path

from crn import Gene
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_1"
CIKTI.mkdir(parents=True, exist_ok=True)

ncbi = Ncbi(organism="human", verbose=False)

print("MEFV Gene ID:", ncbi.search_gene("MEFV"))
mefv = ncbi.gene("MEFV")
print(mefv)

print("\nTanım    :", mefv.features.first("description"))
print("Organizma:", mefv.features.first("organism"))
print("Konum    :", mefv.features.first("map_location"))
print("Takma ad :", ", ".join(mefv.features.get("aliases", [])))
print("Bölge    :", mefv.features.first("genomic_region"))

print("\nGenomik DNA:", mefv.dna)
for kimlik, rna in mefv.transcripts.items():
    print(f"Transkript : {kimlik} {len(rna)} nt, protein: {rna.features.first('protein_id')}")
for kimlik, protein in mefv.proteins.items():
    print(f"Protein    : {kimlik} {len(protein)} aa, {protein.features.first('title')}")

bas = time.perf_counter()
ncbi.gene("MEFV")
print(f"\nYerel depodan okuma: {time.perf_counter() - bas:.3f} sn")

yol = CIKTI / "mefv.crn"
mefv.to_crn(yol)
geri = Gene.from_crn(yol)
print("CRN'den okunan:", geri)
print("Diziler aynı mı:", geri.dna == mefv.dna
      and list(geri.transcripts) == list(mefv.transcripts)
      and all(geri.protein(k) == p for k, p in mefv.proteins.items()))
print("\nÇıktılar:", CIKTI)
