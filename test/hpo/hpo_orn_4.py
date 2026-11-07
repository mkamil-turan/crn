"""Örnek 4: Semptom grubunu Gene ve dizi nesnelerine dönüştürmek.

Senaryo: Nörodejenerasyonla ilişkili genlerin protein dizilerini bir analize, ana
transkriptlerini başka bir analize vereceksiniz. Sembol listesiyle uğraşmak yerine
doğrudan crn nesneleri istiyorsunuz.

Bu örnekte:
  1. Terimden Gene listesi çıkarılır: her gende bir DNA, bir transkript, bir protein.
  2. Gene nesneleri doğrudan CUB'a verilir; ek izoform ayıklaması gerekmez.
  3. Protein ve RNA dizileri FASTA'ya yazılır.
  4. Var olan bir gene HPO bilgileri eklenir, CRN dosyasına yazılıp geri okunur.

Çıktılar: test/hpo/cikti/hpo_orn_4/
Not: GENCODE verisi sys/ içinde hazır değilse ilk çalıştırma onu indirir (büyük).
Protein kodlamayan genler ile mitokondri genleri atlanır ve uyarıyla bildirilir.
"""
from pathlib import Path

from crn import Gene, write_fasta
from crn.cub import CUB
from crn.fetch import Gencode, Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_orn_4"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo()
TERIM = "Neurodegeneration"

genler = hpo.to_genes(TERIM)
print(f"\n{len(genler)} Gene nesnesi. İlk gen:")
gen = genler[0]
print(f"  Ad: {gen.name}")
print(f"  DNA {len(gen.dna)} nt, {len(gen.transcripts)} transkript, {len(gen.proteins)} protein")
print(f"  Transkript seçimi: {gen.features.first('canonical_tag')}")
print(f"  Eşleşen terimler: {', '.join(gen.features.get('hpo_name')[:3])}")
print(f"  Hastalıklar: {', '.join(gen.features.get('disease_name')[:2])}")

cub = CUB(genler)
print(f"\nCUB kuruldu: {len(cub.names)} gen")

proteinler = hpo.to_sequences("protein", TERIM)
write_fasta(proteinler, CIKTI / "proteinler.fasta")
rnalar = hpo.to_sequences("rna", TERIM)
write_fasta(rnalar, CIKTI / "transkriptler.fasta")
print(f"FASTA: {len(proteinler)} protein ({type(proteinler[0]).__name__}), "
      f"{len(rnalar)} transkript ({type(rnalar[0]).__name__})")

scn1a = Gencode(verbose=False).canonical("SCN1A")
hpo.annotate(scn1a, term="Seizure")
scn1a.to_crn(CIKTI / "scn1a.crn")
okunan = Gene.from_crn(CIKTI / "scn1a.crn")
print(f"\nSCN1A dosyadan okundu: {len(okunan.features.get('hpo_id'))} nöbet terimi, "
      f"HPO {okunan.features.first('hpo_release')}")
