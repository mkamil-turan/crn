"""Örnek 8: mRNA üzerinde CDS ve UTR bölgelerini belirle.

Senaryo: Elinizde tam uzunluklu bir mRNA ile bilinen protein dizisi var. mRNA'nın
hangi kısmının proteini kodladığını bulmak ve 5'UTR, CDS, 3'UTR parçalarını ayrı
ayrı saklamak istiyorsunuz (ör. 3'UTR'de miRNA hedefi taramak için).

Bu örnekte:
  1. Sentetik UTR'ler ile INS kodlayan dizisinden bir mRNA kurulur.
  2. match_orf proteini kodlayan okuma çerçevesini bulur.
  3. Bölgeler dilimlenir; 3'UTR'de AAUAAA poliadenilasyon sinyali aranır.
  4. Bilgiler features'a yazılır, CRN'ye kaydedilip geri okunur.
  5. mRNA'yla uyuşmayan protein verilince sonuç gösterilir.

Çıktılar: test/sequence/cikti/seq_orn_8/
Not: UTR'ler sentetiktir, gerçek INS mRNA'sının UTR'leri değildir.
"""
from pathlib import Path

from crn import DNA, PRT, RNA, Sequence, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_8"
CIKTI.mkdir(parents=True, exist_ok=True)

KODLAYAN = (
    "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCC"
    "TTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTC"
    "TTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCT"
    "GGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGT"
    "ACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"
)
PROTEIN = ("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGG"
           "GPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")

# 1. mRNA: 5'UTR + CDS + 3'UTR (poliA sinyali ve kuyruk dahil).
utr5 = "GCAUCGUCAGCCUCAGGAUGACCUCUAGGCUUCAGCUCC"   # içinde çerçevesi uymayan bir AUG var
utr3 = "GCCUGUCCAGUCAUGCCUAACUCAUUCGCCUAAUAAAGCUUCCACAGUCUGCAAAAAAAAAAAAAAAAAAAA"
mrna = RNA(utr5 + DNA(KODLAYAN).transcription(template=False).seq + utr3, name="INS_mRNA")
print("mRNA:", mrna)

# 2. CDS'nin bulunması. Bulunan ORF, proteinin bilgisini ve konumlarını taşır.
protein = PRT(PROTEIN, name="INS_protein")
cds = mrna.match_orf(protein)
if cds is None:
    raise SystemExit("Protein bu mRNA'da kodlanmıyor.")

bas = cds.features.first("orf_start")
son = cds.features.first("orf_end")
print(f"CDS : {bas}-{son} ({len(cds)} nt)")
print("Çeviri proteinle aynı mı:", cds.translation().seq.rstrip("*") == protein.seq)

# 3. Bölgeler. Dilimler 0 tabanlı, orf_start/orf_end 1 tabanlıdır.
bolge_5 = mrna[:bas - 1]
bolge_3 = mrna[son:]
bolge_5.name, cds.name, bolge_3.name = "INS_5UTR", "INS_CDS", "INS_3UTR"

toplam = len(mrna)
print(f"\n{'Bölge':<8}{'nt':>5}{'mRNA %':>9}")
for ad, parca in (("5'UTR", bolge_5), ("CDS", cds), ("3'UTR", bolge_3)):
    print(f"{ad:<8}{len(parca):>5}{100 * len(parca) / toplam:>9.1f}")

sinyal = bolge_3.seq.find("AAUAAA")
print(f"\nKozak bağlamı (ATG -6..+4): {mrna[bas - 7:bas + 3].seq}")
print(f"PoliA sinyali: {'3UTR ' + str(sinyal + 1) + '. konumda' if sinyal >= 0 else 'yok'}")

# 4. Koordinatları mRNA üzerinde saklama ve kayıpsız kaydetme.
mrna.features.set("gene", "INS")
mrna.features.set("utr5", (1, bas - 1))
mrna.features.set("cds", (bas, son))
mrna.features.set("utr3", (son + 1, len(mrna)))
yol = CIKTI / "ins_mrna.crn"
mrna.to_crn(yol)
geri = Sequence.from_crn(yol)
a, b = geri.features.first("cds")
print("\nCRN'den okunan koordinatla CDS yine aynı mı:", geri[a - 1:b] == cds)

write_fasta([bolge_5, cds, bolge_3], CIKTI / "ins_bolgeleri.fasta")

# 5. Bu mRNA'da kodlanmayan protein için sonuç None döner.
print("Uyumsuz protein sonucu:", mrna.match_orf("MVHLTPEEK"))
print("\nÇıktılar:", CIKTI)
