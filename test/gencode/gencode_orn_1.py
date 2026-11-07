"""Örnek 1: GENCODE'a giriş: sürüm, gen, transkript ve protein.

Senaryo: MEFV için GENCODE'daki tüm kayıtları görmek ve klinikte referans alınan
MANE Select transkriptini bulmak istiyorsunuz.

Bu örnekte:
  1. Yerel sürüm kontrol edilir (GENCODE bir kez kurulur, sonra ağ gerekmez).
  2. gene() ile DNA, tüm transkriptler ve proteinler alınır.
  3. transcripts() ve details() ile her transkriptin tipi, ekzon sayısı,
     destek düzeyi (TSL) ve etiketleri tabloya dökülür.
  4. MANE Select transkriptinin bölgeleri çıkarılır, gen CRN'e yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_1/
Not: Kurulu sürüm yoksa ilk kullanım ~1,2 GB indirir.
"""
from pathlib import Path

from crn import Gene
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_1"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)
print("Sürüm:", gc.release, "| hazır mı:", gc.is_ready(), "| kurulu:", gc.installed_releases())

mefv = gc.gene("MEFV")
print(mefv)
print("Tip:", mefv.features.first("gene_type"), "| konum:", mefv.features.first("chromosome"),
      mefv.features.first("start"), "-", mefv.features.first("end"), mefv.features.first("chromosome_strand"))

print(f"\n{'Ad':<11}{'Tip':<32}{'nt':>6}{'Ekzon':>7}{'TSL':>5}  Etiket")
for satir in gc.transcripts(gene="MEFV"):
    ayrinti = gc.details(satir["transcript_id"])
    uzunluk = len(mefv.transcript(satir["transcript_id"]))
    etiket = "MANE_Select" if "MANE_Select" in ayrinti["tags"] else ""
    print(f"{satir['transcript_name']:<11}{satir['transcript_type']:<32}{uzunluk:>6}"
          f"{len(ayrinti['exons']):>7}{ayrinti['support_level'] or '-':>5}  {etiket}")

mane = gc.rna(gc.transcripts(gene="MEFV", tag="MANE_Select")[0]["transcript_id"])
print("\nMANE Select:", mane.name, mane.features.first("transcript_name"))
for bolge in ("utr5", "cds", "utr3"):
    print(f"  {bolge}: {mane.features.first(bolge)}")
print("  Protein:", gc.protein(mane.name))

yol = CIKTI / "mefv_gencode.crn"
mefv.to_crn(yol)
print("\nCRN'den:", Gene.from_crn(yol))
print("Çıktılar:", CIKTI)
