"""Örnek 8: Nonsense aracılı bozunma (NMD) kuralını veride sına.

Senaryo: Klasik kurala göre stop kodonu son ekzon birleşiminin yaklaşık 50-55 nt
yukarısındaysa transkript NMD'ye gider. Bu kuralın GENCODE'un NMD tipi
transkriptlerinde ne kadar geçerli olduğunu MANE kontrolüyle ölçmek istiyorsunuz.

Bu örnekte:
  1. Seçilen kromozomdaki NMD ve MANE Select transkriptleri alınır.
  2. Ekzon koordinatlarından son birleşimin transkript konumu hesaplanır.
  3. Stop kodonunun son birleşime uzaklığı ölçülür.
  4. İki grup karşılaştırılır, ölçümler TSV'ye yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_8/
Not: GENCODE'un NMD tipi atamasında bu kural kullanıldığından sonuç biyolojik bir
kanıt değil, anotasyonun iç tutarlılığının ölçüsüdür.
"""
import csv
from pathlib import Path
from statistics import median

from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_8"
CIKTI.mkdir(parents=True, exist_ok=True)

KROMOZOM = "chr17"
ESIK = 50
EKSIK = {"mRNA_end_NF", "cds_end_NF", "cds_start_NF", "mRNA_start_NF"}

gc = Gencode(verbose=False)


def olc(satir: dict) -> int | None:
    """Stop'un son birleşime uzaklığı; artı değer stop'un birleşimden önce olduğunu gösterir."""
    ayrinti = gc.details(satir["transcript_id"])
    if EKSIK & set(ayrinti["tags"]) or len(ayrinti["exons"]) < 2:
        return None
    cds = gc.rna(satir["transcript_id"]).features.first("cds")
    if cds is None:
        return None
    # Ekzonlar transkript sırasındadır; sonuncusu hariç toplam, son birleşimin konumudur.
    son_birlesim = sum(b - a + 1 for a, b in ayrinti["exons"][:-1])
    return son_birlesim - cds[1]


gruplar = {"NMD": gc.transcripts(transcript_type="nonsense_mediated_decay", chromosome=KROMOZOM),
           "MANE": gc.transcripts(transcript_type="protein_coding", chromosome=KROMOZOM, tag="MANE_Select")}

print(f"{KROMOZOM}\n{'Grup':<7}{'Kayıt':>7}{'Ölçülen':>9}{'Ortanca nt':>12}{f'>{ESIK} nt %':>10}")
satirlar = []
for ad, liste in gruplar.items():
    olcumler = [(t, olc(t)) for t in liste]
    gecerli = [(t, m) for t, m in olcumler if m is not None]
    degerler = [m for _, m in gecerli]
    print(f"{ad:<7}{len(liste):>7,}{len(degerler):>9,}{median(degerler):>12,.0f}"
          f"{100 * sum(m > ESIK for m in degerler) / len(degerler):>10.1f}")
    satirlar += [{"grup": ad, "transkript": t["transcript_name"], "uzaklik_nt": m} for t, m in gecerli]

with open(CIKTI / "nmd_uzakliklari.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)
print("\nÇıktılar:", CIKTI)
