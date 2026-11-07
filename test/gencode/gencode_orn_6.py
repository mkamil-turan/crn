"""Örnek 6: Psödogenler: tipleri, aileleri ve PTEN/PTENP1 örneği.

Senaryo: RNA-seq okumalarının psödogenlere kayabileceğinden şüpheleniyorsunuz.
Hangi genlerin çok psödogeni olduğunu ve PTEN ile PTENP1 benzerliğini görmek
istiyorsunuz.

Bu örnekte:
  1. Psödogen tipleri sayılır, işlenmiş ve işlenmemiş ayrımı gösterilir.
  2. Ad kalıbıyla (RPL5P3 -> RPL5) kaynak genlere göre gruplanır.
  3. Çok eşleşen PTENP1 adı gene_id ile çözülür.
  4. PTEN mRNA'sı ile PTENP1 arasındaki paylaşılan 21-mer oranı hesaplanır.
  5. En uzun ORF'ler karşılaştırılır: uzun ORF kalması kodlama kanıtı değildir.
  6. RNA'lar tek CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_6/
Not: Ad kalıbı kaba bir sezgidir; kesin ebeveyn gen için hizalama gerekir.
"""
import re
from collections import Counter
from pathlib import Path

from crn import RNA, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_6"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)

tipler = {t: n for t, n in gc.gene_types().items() if "pseudogene" in t}
print(f"Psödogen: {sum(tipler.values()):,} gen")
for tip, sayi in tipler.items():
    if sayi >= 50:
        print(f"  {tip:<40}{sayi:>7,}")

psodogenler = gc.genes(gene_type=tipler)
kaynak = Counter()
for g in psodogenler:
    eslesme = re.fullmatch(r"(.+?)P\d+", g["gene_name"])
    if eslesme:
        kaynak[eslesme.group(1)] += 1
print("\nAdına göre en çok psödogeni olan genler:")
print("  " + ", ".join(f"{a} ({n})" for a, n in kaynak.most_common(10)))

cevrilen = gc.genes(gene_type="translated_processed_pseudogene")
print("\nÇevrildiği bilinen işlenmiş psödogenler:", ", ".join(g["gene_name"] for g in cevrilen))

# PTENP1 adı iki gene karşılık gelir; ad ile sorgu hata verir.
try:
    gc.gene("PTENP1")
except ValueError as hata:
    print("\nBelirsiz ad:", hata)
kayit = gc.genes(gene_type="transcribed_processed_pseudogene", name="^PTENP1$")[0]
psodo = gc.rna(gc.transcripts(gene=kayit["gene_id"])[0]["transcript_id"])
psodo.name = "PTENP1"
pten = gc.rna(gc.transcripts(gene="PTEN", tag="MANE_Select")[0]["transcript_id"])
pten.name = "PTEN_MANE"
print(f"PTEN mRNA {len(pten):,} nt | PTENP1 {len(psodo):,} nt")

K = 21


def kmerler(rna: RNA) -> set:
    return {p.seq for p in rna.sliding(K)}


pten_k, psodo_k = kmerler(pten), kmerler(psodo)
print(f"PTENP1 {K}-merlerinin PTEN'de bulunanı: %{100 * len(pten_k & psodo_k) / len(psodo_k):.1f}")
print("Bu nedenle bu bölgelere düşen kısa okumalar iki kayıt arasında ayrılamayabilir.")


def en_uzun_orf(rna: RNA) -> int:
    sens = rna.reverse_transcription().complementary()
    return max((len(o) // 3 - 1 for o in sens.find_orfs() if o.features.first("orf_strand") == "+"), default=0)


print(f"\nEn uzun ORF: PTEN {en_uzun_orf(pten)} aa, PTENP1 {en_uzun_orf(psodo)} aa")
print("PTENP1'in ORF'si uzun kalmıştır, yine de GENCODE kaydında proteini yoktur:",
      kayit["gene_type"])

write_crn_list([pten, psodo], CIKTI / "pten_ptenp1.crn")
print("Çıktılar:", CIKTI)
