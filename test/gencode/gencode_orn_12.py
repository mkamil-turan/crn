"""Örnek 12: Ekzon-intron sınırları ve GT-AG kuralı.

Senaryo: Ekzon koordinatlarından intronları kurup çoğu intronun GT ile başlayıp
AG ile bittiğini gerçek genomda doğrulamak, HBB intronlarını da saklamak istiyorsunuz.

Bu örnekte:
  1. MANE Select transkriptlerinin ekzonlarından intronlar türetilir.
  2. Her intronun ilk ve son iki bazı zincire uygun olarak region() ile alınır.
  3. Sınır dinükleotitleri sayılır.
  4. HBB ekzon ve intronları DNA olarak tek CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_12/
"""
from collections import Counter
from pathlib import Path

from crn import DNA, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_12"
CIKTI.mkdir(parents=True, exist_ok=True)

KROMOZOMLAR = ["chr21", "chr22"]

gc = Gencode(verbose=False)


def intronlar(t: dict) -> list:
    """Transkript sırasıyla (baş, bit) intronlar; bas < bit genomik koordinattır."""
    ekzonlar = gc.details(t["transcript_id"])["exons"]
    if t["strand"] == "+":
        return [(a[1] + 1, b[0] - 1) for a, b in zip(ekzonlar, ekzonlar[1:])]
    return [(b[1] + 1, a[0] - 1) for a, b in zip(ekzonlar, ekzonlar[1:])]


def uclar(krom: str, bas: int, bit: int, zincir: str) -> tuple:
    """(donör, akseptör) dinükleotitleri; her ikisi de transkript yönünde."""
    sol = gc.region(krom, bas, bas + 1, zincir).seq
    sag = gc.region(krom, bit - 1, bit, zincir).seq
    return (sol, sag) if zincir == "+" else (sag, sol)


sayac = Counter()
for krom in KROMOZOMLAR:
    for t in gc.transcripts(transcript_type="protein_coding", chromosome=krom, tag="MANE_Select"):
        for bas, bit in intronlar(t):
            donor, akseptor = uclar(krom, bas, bit, t["strand"])
            sayac[f"{donor}-{akseptor}"] += 1

toplam = sum(sayac.values())
print(f"{', '.join(KROMOZOMLAR)}: {toplam:,} intron")
for tip, sayi in sayac.most_common(6):
    print(f"  {tip}  {sayi:>7,}  %{100 * sayi / toplam:.2f}")

t = gc.transcripts(gene="HBB", tag="MANE_Select")[0]
parcalar = []
for no, (bas, bit) in enumerate(gc.details(t["transcript_id"])["exons"], start=1):
    parcalar.append(gc.region(t["chromosome"], bas, bit, t["strand"], name=f"HBB_ekzon{no}"))
for no, (bas, bit) in enumerate(intronlar(t), start=1):
    parcalar.append(gc.region(t["chromosome"], bas, bit, t["strand"], name=f"HBB_intron{no}"))

print("\nHBB parçaları:")
for p in parcalar:
    print(f"  {p.name:<12}{len(p):>5} nt  {p.seq[:6]}...{p.seq[-6:]}")

ekzon_birlesik = "".join(p.seq for p in parcalar if "ekzon" in p.name)
print("Ekzonların birleşimi mRNA ile aynı mı:",
      gc.rna(t["transcript_id"]).reverse_transcription().complementary() == DNA(ekzon_birlesik))

write_crn_list(parcalar, CIKTI / "hbb_ekzon_intron.crn")
print("Çıktılar:", CIKTI)
