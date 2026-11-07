"""Örnek 10: Bir lokusta neler var? Beta-globin bölgesi.

Senaryo: GWAS ya da kopya sayısı analizinde chr11'de bir bölge çıktı. Bölgedeki
tüm genleri tipleriyle, sıraları ve zincirleriyle görmek istiyorsunuz.

Bu örnekte:
  1. Koordinatla örtüşen genler genes() ile sorgulanır.
  2. Kromozom sırasına dizilerek metin haritası çizilir. Gen sınırı tek bir uzak
     ekzonlu transkriptle çok genişleyebildiği için (HBG2 252 kb, MANE 1,6 kb)
     haritada MANE Select sınırları kullanılır.
  3. Bölgenin DNA'sı alınır, HBB'nin konumu dizide doğrulanır.
  4. Globin genleri ve bölge DNA'sı aynı CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_10/
"""
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_10"
CIKTI.mkdir(parents=True, exist_ok=True)

KROMOZOM, BAS, BIT = "chr11", 5_220_000, 5_310_000

gc = Gencode(verbose=False)
genler = sorted(gc.genes(chromosome=KROMOZOM, start=BAS, end=BIT), key=lambda g: g["start"])
print(f"{KROMOZOM}:{BAS:,}-{BIT:,} ({(BIT - BAS) / 1000:.0f} kb), {len(genler)} gen\n")
print(f"{'Gen':<16}{'Tip':<36}{'Başlangıç':>12}{'kb':>7}  Z")
for g in genler:
    print(f"{g['gene_name']:<16}{g['gene_type']:<36}{g['start']:>12,}{(g['end'] - g['start'] + 1) / 1000:>7.1f}  {g['strand']}")

print("\nHarita (MANE Select, her karakter 1 kb, < - zincir, > + zincir):")
for g in genler:
    mane = gc.transcripts(gene=g["gene_id"], tag="MANE_Select")
    if mane and mane[0]["start"] <= BIT:
        t = mane[0]
        sol = max(0, (t["start"] - BAS) // 1000)
        uz = max(1, (min(BIT, t["end"]) - max(BAS, t["start"])) // 1000)
        isaret = ">" if t["strand"] == "+" else "<"
        print(f"{g['gene_name']:<10}{' ' * sol}{isaret * uz}  {t['end'] - t['start'] + 1:,} nt")

bolge = gc.region(KROMOZOM, BAS, BIT, "+", name="beta_globin_lokusu")
hbb = next(g for g in genler if g["gene_name"] == "HBB")
dilim = bolge[hbb["start"] - BAS:hbb["end"] - BAS + 1]
print("\nDilimin ters tümleyeni gc.dna('HBB') ile aynı mı:", dilim.reverse_complement() == gc.dna("HBB"))

globinler = [gc.gene(g["gene_id"]) for g in genler if g["gene_type"] == "protein_coding" and g["gene_name"].startswith("HB")]
yol = CIKTI / "beta_globin_lokusu.crn"
write_crn_list([bolge] + globinler, yol)
print("Paket:", [type(x).__name__ + ":" + x.name for x in read_crn_list(yol)])
print("Çıktılar:", CIKTI)
