"""Örnek 3: Bir geni izoformlarına ayır ve tek dosyada sakla.

Senaryo: TP53'ün birçok transkripti var; aynı proteini veren kayıtları ayıklayıp
her biri bir DNA, bir transkript ve bir protein içeren temiz izoformlar istiyorsunuz.

Bu örnekte:
  1. TP53 alınır; doğrulanmış (NM_) ve tahmini (XM_) kayıtlar sayılır.
  2. isoforms() ile izoformlar ayrılır; bildirimler yakalanıp raporlanır.
  3. Her izoformun ORF konumu ve protein uzunluğu tabloya dökülür.
  4. İzoform listesi tek CRN dosyasına yazılır ve geri okunur.

Çıktılar: test/ncbi/cikti/ncbi_orn_3/
"""
import warnings
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_3"
CIKTI.mkdir(parents=True, exist_ok=True)

tp53 = Ncbi(verbose=False).gene("TP53")
print(tp53)
nm = [k for k in tp53.transcripts if k.startswith("NM_")]
xm = [k for k in tp53.transcripts if k.startswith("XM_")]
print(f"Doğrulanmış (NM_): {len(nm)}, tahmini (XM_): {len(xm)}")

with warnings.catch_warnings(record=True) as bildirimler:
    warnings.simplefilter("always")
    izoformlar = tp53.isoforms()

print(f"\nBenzersiz CDS veren izoform: {len(izoformlar)}")
print(f"{'Transkript':<16}{'Protein':<16}{'ORF':>12}{'aa':>6}")
for izoform in izoformlar:
    tid = izoform.features.first("transcript_id")
    pid = izoform.features.first("protein_id")
    rna = izoform.transcript(tid)
    orf = f"{rna.features.first('orf_start')}-{rna.features.first('orf_end')}"
    print(f"{tid:<16}{pid:<16}{orf:>12}{len(izoform.protein(pid)):>6}")

print("\nBildirimler:")
for bildirim in bildirimler:
    print(" -", str(bildirim.message)[:110], "...")

yol = CIKTI / "tp53_izoformlari.crn"
write_crn_list(izoformlar, yol)
geri = read_crn_list(yol)
print(f"\n{len(geri)} izoform okundu, dosya {yol.stat().st_size / 1e6:.1f} MB")
print("İlk kayıt:", geri[0])
print("Çıktılar:", CIKTI)
