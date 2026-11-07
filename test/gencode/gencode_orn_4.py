"""Örnek 4: Uzun kodlamayan RNA'lar (lncRNA) nasıl görünür?

Senaryo: Çalışmanızdaki lncRNA'ların GENCODE'daki durumunu denetlemek ve bunları
genel lncRNA profiliyle karşılaştırmak istiyorsunuz.

Bu örnekte:
  1. lncRNA olarak bilinen genler alınır; güncel gen tipi ve protein kaydı denetlenir.
  2. En uzun izoformda en uzun açık okuma çerçevesi bulunur.
  3. Genom ölçeğinde lncRNA ve mRNA'nın ekzon sayıları karşılaştırılır.
  4. Genler tek CRN listesine paketlenir.

Çıktılar: test/gencode/cikti/gencode_orn_4/
Not: ORF bulunması protein kanıtı değildir. Ters yönde, adı lncRNA diye bilinen bir gen
de yeni sürümlerde protein kodlayan olarak sınıflanabilir (TUG1 buna örnektir).
"""
from collections import Counter
from pathlib import Path
from statistics import median

from crn import read_crn_list, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_4"
CIKTI.mkdir(parents=True, exist_ok=True)

LNCRNA = ["MALAT1", "NEAT1", "XIST", "HOTAIR", "H19", "MEG3", "TUG1", "KCNQ1OT1"]

gc = Gencode(verbose=False)
print(f"{'Gen':<11}{'Tip':<16}{'Krom':<7}{'İzoform':>8}{'En uzun nt':>12}{'GC%':>7}{'Protein':>9}{'En uzun ORF aa':>16}")
genler = []
for ad in LNCRNA:
    gen = gc.gene(ad)
    genler.append(gen)
    rna = max(gen.transcripts.values(), key=len)
    gc_yuzde = 100 * (rna.seq.count("G") + rna.seq.count("C")) / len(rna)

    # find_orfs iki zinciri tarar; transkript yönündekiler (+) alınır.
    sens = rna.reverse_transcription().complementary()
    orflar = [o for o in sens.find_orfs() if o.features.first("orf_strand") == "+"]
    orf_aa = max((len(o) // 3 - 1 for o in orflar), default=0)
    print(f"{ad:<11}{gen.features.first('gene_type'):<16}{gen.features.first('chromosome'):<7}{len(gen.transcripts):>8}{len(rna):>12,}"
          f"{gc_yuzde:>7.1f}{len(gen.proteins):>9}{orf_aa:>16}")


def ekzon_profili(tip: str) -> tuple:
    sayilar = [len(gc.details(t["transcript_id"])["exons"]) for t in gc.transcripts(transcript_type=tip)]
    return len(sayilar), median(sayilar), 100 * sum(s == 1 for s in sayilar) / len(sayilar)


print(f"\n{'Tip':<18}{'Transkript':>12}{'Ortanca ekzon':>15}{'Tek ekzon %':>13}")
for tip in ("lncRNA", "protein_coding"):
    n, ort, tek = ekzon_profili(tip)
    print(f"{tip:<18}{n:>12,}{ort:>15.0f}{tek:>13.1f}")

kromozom = Counter(g["chromosome"] for g in gc.genes(gene_type="lncRNA"))
print("\nEn çok lncRNA:", ", ".join(f"{k} ({v:,})" for k, v in kromozom.most_common(3)))

yol = CIKTI / "lncrna_paneli.crn"
write_crn_list(genler, yol)
print(f"\nPaket: {len(read_crn_list(yol))} gen, {yol.stat().st_size / 1e6:.1f} MB")
print("Çıktılar:", CIKTI)
