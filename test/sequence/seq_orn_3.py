"""Örnek 3: Bir contig'de gen adaylarını bul ve proteinlerini çıkar.

Senaryo: Dizileme sonrası elinizde anotasyonsuz bir contig var. Altı okuma
çerçevesinde açık okuma çerçevelerini (ORF) tarayıp olası proteinleri bulmak
istiyorsunuz.

Bu örnekte:
  1. Sentetik bir contig FASTA'sı üretilir: rastgele arka plana bir gen düz,
     bir gen ters zincire yerleştirilir. Böylece doğru yanıt bilinir.
  2. FASTA okunur ve find_orfs ile iki zincir taranır.
  3. İç içe ORF'lerden aynı stop kodonunu paylaşanların en uzunu bırakılır.
  4. ORF'ler proteine çevrilir, FASTA ve TSV olarak kaydedilir.
  5. Yerleştirilen genlerin bulunup bulunmadığı kontrol edilir.

Çıktılar: test/sequence/cikti/seq_orn_3/
Not: ORF bulmak gen kanıtı değildir; adayları BLAST gibi araçlarla doğrulayın.
"""
import csv
import random
from pathlib import Path

from crn import DNA, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_3"
CIKTI.mkdir(parents=True, exist_ok=True)
EN_KISA_ORF = 240  # nt, stop dahil; yaklaşık 80 amino asit

HBB = ("ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTT"
       "GGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGAT"
       "CTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTT"
       "AGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAG"
       "CTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
       "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCAC"
       "AAGTATCACTAA")
INS = ("ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCC"
       "TTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTC"
       "TTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCT"
       "GGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGT"
       "ACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG")

# 1. Sentetik contig. Tohum sabit olduğu için her çalıştırmada aynı dizi oluşur.
rastgele = random.Random(42)


def arka_plan(uzunluk: int) -> str:
    return "".join(rastgele.choices("ACGT", k=uzunluk))


contig = DNA(arka_plan(500) + HBB + arka_plan(400)
             + DNA(INS).reverse_complement().seq + arka_plan(300), name="contig_1")
fasta_yolu = CIKTI / "contig.fasta"
contig.to_fasta(fasta_yolu)

# 2. Okuma ve tarama.
contig = DNA.from_fasta_file(fasta_yolu)[0]
orflar = contig.find_orfs(minimum_length_of_orf=EN_KISA_ORF)
print(f"{contig.name}: {len(contig)} nt, {len(orflar)} ORF (iç içe olanlar dahil)")


# 3. Aynı stop'u paylaşan ORF'ler aynı genin farklı ATG'leridir; en uzunu tutulur.
def stop_konumu(orf: DNA) -> tuple:
    zincir = orf.features.first("orf_strand")
    konum = orf.features.first("orf_end" if zincir == "+" else "orf_start")
    return zincir, konum


en_uzunlar: dict = {}
for orf in orflar:
    anahtar = stop_konumu(orf)
    if anahtar not in en_uzunlar or len(orf) > len(en_uzunlar[anahtar]):
        en_uzunlar[anahtar] = orf
adaylar = sorted(en_uzunlar.values(), key=len, reverse=True)

# 4. Proteine çevirme ve kaydetme. find_orfs diziyi zaten kodlayan yönde verir.
proteinler = []
satirlar = []
for orf in adaylar:
    protein = orf.transcription(template=False).translation()
    proteinler.append(protein)
    satirlar.append({
        "orf": orf.name,
        "baslangic": orf.features.first("orf_start"),
        "bitis": orf.features.first("orf_end"),
        "zincir": orf.features.first("orf_strand"),
        "cerceve": orf.features.first("orf_frame"),
        "uzunluk_nt": len(orf),
        "uzunluk_aa": len(protein) - 1,
        "protein": protein.seq,
    })

write_fasta(proteinler, CIKTI / "aday_proteinler.fasta")
with open(CIKTI / "orf_tablosu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)

print(f"Benzersiz aday: {len(adaylar)}\n")
print(f"{'ORF':<16}{'Konum':>14}  Z  Ç  {'aa':>4}  Başlangıç")
for s in satirlar:
    konum = f"{s['baslangic']}-{s['bitis']}"
    print(f"{s['orf']:<16}{konum:>14}  {s['zincir']}  {s['cerceve']}  "
          f"{s['uzunluk_aa']:>4}  {s['protein'][:12]}")

# 5. Doğrulama: yerleştirdiğimiz genler bulundu mu?
print()
for ad, gen in (("HBB (+ zincir)", HBB), ("INS (- zincir)", INS)):
    beklenen = DNA(gen).transcription(template=False).translation().seq
    print(f"{ad:<16}: {'bulundu' if any(p.seq == beklenen for p in proteinler) else 'BULUNAMADI'}")
print("\nÇıktılar:", CIKTI)
