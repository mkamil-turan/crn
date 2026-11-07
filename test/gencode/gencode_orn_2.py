"""Örnek 2: Genom envanteri: hangi kromozomda ne kadar gen var?

Senaryo: Bir derleme yazısı için insan genomundaki gen tiplerinin dağılımını
sürüm numarasıyla birlikte tablo olarak vermek istiyorsunuz.

Bu örnekte:
  1. gene_types() ile tüm gen tipleri sayılır.
  2. Tipler biyolojik gruplara toplanır.
  3. genes() ile kromozom-grup tablosu oluşturulur.
  4. Sonuçlar TSV olarak kaydedilir.

Çıktılar: test/gencode/cikti/gencode_orn_2/
"""
import csv
from collections import Counter
from pathlib import Path

from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_2"
CIKTI.mkdir(parents=True, exist_ok=True)

KUCUK_RNA = {"miRNA", "snRNA", "snoRNA", "scaRNA", "sRNA", "misc_RNA", "rRNA", "ribozyme",
             "vault_RNA", "Mt_tRNA", "Mt_rRNA"}
GRUPLAR = ["protein kodlayan", "lncRNA", "küçük RNA", "psödogen", "immün (IG/TR)", "diğer"]


def grup(gen_tipi: str) -> str:
    if gen_tipi == "protein_coding":
        return "protein kodlayan"
    if gen_tipi == "lncRNA":
        return "lncRNA"
    if gen_tipi in KUCUK_RNA:
        return "küçük RNA"
    if "pseudogene" in gen_tipi:
        return "psödogen"
    if gen_tipi.startswith(("IG_", "TR_")):
        return "immün (IG/TR)"
    return "diğer"


gc = Gencode(verbose=False)
tipler = gc.gene_types()
toplam = sum(tipler.values())
print(f"GENCODE {gc.release}: {toplam:,} gen, {len(tipler)} tip\n")

grup_sayilari = Counter()
for tip, sayi in tipler.items():
    grup_sayilari[grup(tip)] += sayi
for ad in GRUPLAR:
    print(f"{ad:<18}{grup_sayilari[ad]:>8,}  %{100 * grup_sayilari[ad] / toplam:.1f}")

print("\nEn yaygın 8 tip:")
for tip, sayi in list(tipler.items())[:8]:
    print(f"  {tip:<38}{sayi:>8,}")

tablo: dict = {}
for gen in gc.genes():
    tablo.setdefault(gen["chromosome"], Counter())[grup(gen["gene_type"])] += 1
kromozomlar = [c for c in tablo if c.startswith("chr")]

KISA = ["Protein", "lncRNA", "Küçük RNA", "Psödogen", "İmmün", "Diğer"]
print(f"\n{'Krom':<7}" + "".join(f"{k:>11}" for k in KISA))
for krom in kromozomlar:
    print(f"{krom:<7}" + "".join(f"{tablo[krom][g]:>11,}" for g in GRUPLAR))
print(f"\nTamamlanmamış contig'lerdeki gen: {sum(sum(v.values()) for k, v in tablo.items() if not k.startswith('chr')):,}")

with open(CIKTI / "kromozom_gen_tablosu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.writer(dosya, delimiter="\t")
    yazici.writerow(["kromozom", *GRUPLAR])
    for krom, sayac in tablo.items():
        yazici.writerow([krom, *(sayac[g] for g in GRUPLAR)])
print("Çıktılar:", CIKTI)
