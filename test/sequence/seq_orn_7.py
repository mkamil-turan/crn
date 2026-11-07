"""Örnek 7: Varyantların protein üzerindeki etkisini bul.

Senaryo: Bir hastada HBB geninde saptanan varyantların sessiz, yanlış anlamlı,
anlamsız ya da çerçeve kayması olup olmadığını öğrenmek istiyorsunuz.

Bu örnekte:
  1. Varyantlar (cDNA konumu, referans, alternatif) tablosu tanımlanır.
  2. Referans bazın yerinde olduğu doğrulanır.
  3. Her varyant için yeni dizi kurulur, proteine çevrilir, referansla
     karşılaştırılır ve sonuç sınıflandırılır.
  4. Sonuç HGVS benzeri gösterimle TSV'ye, varyant proteinler FASTA'ya yazılır.

Yerleştirilen varyantlar: c.20A>T orak hücre (HbS), c.19G>A HbC, c.118C>T
beta-talasemi nonsense, c.9T>C yaygın sessiz varyant, c.20del çerçeve kayması.

Çıktılar: test/sequence/cikti/seq_orn_7/
Not: Sonuç yalnızca protein dizisindeki değişimi söyler; klinik yorum için
ClinVar gibi kaynaklara bakın.
"""
import csv
from pathlib import Path

from crn import DNA, PRT, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_7"
CIKTI.mkdir(parents=True, exist_ok=True)

UC_HARF = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
           "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
           "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
           "Y": "Tyr", "V": "Val", "*": "Ter"}

referans = DNA(
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTT"
    "GGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGAT"
    "CTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTT"
    "AGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAG"
    "CTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
    "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCAC"
    "AAGTATCACTAA",
    name="HBB_CDS",
)

# (etiket, cDNA konumu [1 tabanlı], referans baz, alternatif; "" silinme demektir)
VARYANTLAR = [("HbS", 20, "A", "T"), ("HbC", 19, "G", "A"), ("beta0_CD39", 118, "C", "T"),
              ("sessiz_c9", 9, "T", "C"), ("cerceve_c20del", 20, "A", "")]


def cevir(dna: DNA) -> PRT:
    """Dizideki ilk stop kodonunda biten proteini verir (stop işareti dahil)."""
    protein = dna.transcription(template=False).translation()
    stop = protein.seq.find("*")
    return protein if stop == -1 else protein[:stop + 1]


def siniflandir(ref: str, alt: str, degisim: str) -> tuple:
    """Protein karşılaştırması: (etki sınıfı, HGVS benzeri protein gösterimi)."""
    if ref == alt:
        return "sessiz", "p.="
    i = next((i for i, (a, b) in enumerate(zip(ref, alt)) if a != b), min(len(ref), len(alt)))
    if len(degisim) % 3:
        ter = len(alt) - i
        yeni = UC_HARF[alt[i]] if i < len(alt) else "?"
        return "çerçeve kayması", f"p.{UC_HARF[ref[i]]}{i + 1}{yeni}fsTer{ter}"
    if i < len(alt) and alt[i] == "*":
        return "anlamsız (nonsense)", f"p.{UC_HARF[ref[i]]}{i + 1}Ter"
    return "yanlış anlamlı (missense)", f"p.{UC_HARF[ref[i]]}{i + 1}{UC_HARF[alt[i]]}"


ref_protein = cevir(referans)
print(f"Referans protein: {len(ref_protein) - 1} aa\n")
print(f"{'Varyant':<16}{'cDNA':<10}{'Etki':<27}{'Protein':<22}{'aa':>4}")

satirlar = []
proteinler = [ref_protein]
for etiket, konum, ref_baz, alt_baz in VARYANTLAR:
    i = konum - 1
    if referans[i] != ref_baz:
        raise ValueError(f"{etiket}: {konum}. konumda {ref_baz} beklendi, {referans[i]} bulundu")

    varyant = referans[:i]
    varyant.seq = referans.seq[:i] + alt_baz + referans.seq[i + 1:]
    varyant.name = f"HBB_{etiket}"

    protein = cevir(varyant)
    protein.name = f"HBB_{etiket}_protein"
    proteinler.append(protein)

    etki, gosterim = siniflandir(ref_protein.seq, protein.seq, ref_baz if not alt_baz else "")
    cdna = f"c.{konum}{ref_baz}>{alt_baz}" if alt_baz else f"c.{konum}del"
    aa = len(protein) - 1
    print(f"{etiket:<16}{cdna:<10}{etki:<27}{gosterim:<22}{aa:>4}")
    satirlar.append({"varyant": etiket, "cdna": cdna, "etki": etki,
                     "protein": gosterim, "protein_uzunlugu": aa})

write_fasta(proteinler, CIKTI / "varyant_proteinler.fasta")
with open(CIKTI / "varyant_etkileri.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)
print("\nÇıktılar:", CIKTI)
