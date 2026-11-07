"""Örnek 2: Santral dogma: DNA -> mRNA -> protein ve geri cDNA.

Senaryo: İnsan insülin (INS) kodlayan dizisinden proteini üretmek ve her adımı
tek bir FASTA dosyasında toplamak istiyorsunuz.

Bu örnekte:
  1. Kodlayan zincirden mRNA yapılır (template=False).
  2. Aynı mRNA'nın kalıp zincirden de üretilebildiği gösterilir (template=True).
  3. mRNA proteine çevrilir; sinyal peptidi ve A/B zincirleri dilimlenir.
  4. Ters transkripsiyonla cDNA yapılır, başlangıç dizisine dönülür.
  5. Bir mitokondri geni iki genetik kodla okunarak fark gösterilir.
  6. Tüm ürünler FASTA dosyasına yazılır.

Çıktılar: test/sequence/cikti/seq_orn_2/
Not: Dizi örnek amaçlıdır; gerçek analizde NCBI'dan güncel kaydı alın.
"""
from pathlib import Path

from crn import DNA, GeneticCode, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_2"
CIKTI.mkdir(parents=True, exist_ok=True)

kodlayan = DNA(
    "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCC"
    "TTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTC"
    "TTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCT"
    "GGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGT"
    "ACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG",
    name="INS_kodlayan",
)
kodlayan.features.set("gene", "INS")

# 1. Kodlayan zincir: mRNA ile aynı metin, yalnızca T yerine U.
mrna = kodlayan.transcription(template=False)
mrna.name = "INS_mRNA"
print("DNA :", kodlayan)
print("mRNA:", mrna)

# 2. Hücrede RNA polimeraz kalıp zinciri okur. Kalıp 3'->5' yazıldığında sonuç aynıdır.
kalip = kodlayan.complementary()
print("\nKalıp zincir yönü:", kalip.strand.value[:2])
print("Kalıptan aynı mRNA çıktı mı:", kalip.transcription(template=True) == mrna)

# 3. Translasyon ve olgunlaşma. Bu sınırlar INS için bilinen dilimlerdir.
protein = mrna.translation()
protein.name = "INS_preproinsulin"
sinyal = protein[:24]
b_zinciri = protein[24:54]
c_peptit = protein[54:89]
a_zinciri = protein[89:110]
for ad, parca in (("sinyal_peptidi", sinyal), ("B_zinciri", b_zinciri),
                  ("C_peptidi", c_peptit), ("A_zinciri", a_zinciri)):
    parca.name = f"INS_{ad}"

print("\nProtein    :", protein)
print("Sinyal     :", sinyal.seq)
print("B zinciri  :", b_zinciri.seq)
print("A zinciri  :", a_zinciri.seq)
print("Sistein sayısı (A+B):", a_zinciri.seq.count("C") + b_zinciri.seq.count("C"),
      "-> 3 disülfit bağı için 6 sistein gerekir")

# 4. Ters transkripsiyon: ilk zincir cDNA, mRNA'ya antiparalel ve tümleyicidir.
cdna = mrna.reverse_transcription()
cdna.name = "INS_cDNA_ilk_zincir"
ikinci = cdna.complementary()
ikinci.name = "INS_cDNA_ikinci_zincir"
print("\ncDNA yönü  :", cdna.strand.value[:2])
print("İkinci zincir başlangıç diziyle aynı mı:", ikinci == kodlayan)

# 5. Genetik kod: aynı kodon mitokondride başka anlam taşır. MT-CO1'in ilk kısmı
# standart kodla okunursa 6. kodondaki UGA'da erken stop görünür; oysa mitokondride
# UGA triptofan kodlar. Mitokondriyal genleri yanlış kodla çevirmek sık yapılan bir hatadır.
mt_co1 = DNA(
    "ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACCTATTATTCGGCGCA"
    "TGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGAGCCGAGCTGGGCCAGCCAGGCAACCTTCTA"
    "GGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTTGTAATAATCTTCTTCATAGTAATACCC"
    "ATCATAATCGGAGGCTTTGGCAACTGACTAGTTCCCCTAATAATCGG",
    name="MT_CO1_ilk_kisim",
).transcription(template=False)
standart = mt_co1.translation().seq
mito = mt_co1.translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL).seq
print("\nMT-CO1 standart kod     :", standart[:30], "...")
print("MT-CO1 mitokondriyal kod:", mito[:30], "...")
print("Standart kodda stop sayısı:", standart.count("*"),
      "| mitokondriyal kodda:", mito.count("*"))
for i, (a, b) in enumerate(zip(standart, mito)):
    if a != b:
        print(f"  {i + 1:>3}. kodon {mt_co1.seq[i * 3:i * 3 + 3]}: standart {a}, mitokondri {b}")

# 6. Hepsini tek dosyada topla. Molekül tipleri farklı olsa da FASTA hepsini tutar.
yol = CIKTI / "ins_urunleri.fasta"
write_fasta([kodlayan, mrna, cdna, protein, sinyal, b_zinciri, c_peptit, a_zinciri], yol)
print("\nYazıldı:", yol)
