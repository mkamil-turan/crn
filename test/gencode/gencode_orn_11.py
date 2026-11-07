"""Örnek 11: Antisens lncRNA'lar ve çift yönlü promotörler.

Senaryo: Genom genelinde iki düzenleyici mimariyi saymak istiyorsunuz: protein
kodlayan genle karşı zincirde örtüşen lncRNA'lar ve başlangıçları birbirine çok
yakın, zıt yönde okunan gen çiftleri.

Bu örnekte:
  1. Kromozom başına protein kodlayan genler ve lncRNA'lar sıralanır.
  2. Karşı zincirde örtüşen lncRNA'lar ikili aramayla bulunur.
  3. TSS'leri 1 kb içinde, ters yönlü gen çiftleri sayılır.
  4. Sonuçlar TSV'ye yazılır, TP53 komşuluğu örnek olarak gösterilir.

Çıktılar: test/gencode/cikti/gencode_orn_11/
Not: Gen sınırları en uç transkriptlere göredir; TSS değerleri yaklaşıktır.
"""
import csv
from bisect import bisect_right
from pathlib import Path

from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_11"
CIKTI.mkdir(parents=True, exist_ok=True)

MESAFE = 1000

gc = Gencode(verbose=False)
kromozomlar = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def tss(g: dict) -> int:
    return g["start"] if g["strand"] == "+" else g["end"]


antisens, cift_yonlu, satirlar = 0, 0, []
print(f"{'Krom':<7}{'Protein':>9}{'lncRNA':>9}{'Antisens':>10}{'Çift yönlü':>12}")
for krom in kromozomlar:
    kodlayan = sorted(gc.genes(gene_type="protein_coding", chromosome=krom), key=lambda g: g["start"])
    lncler = gc.genes(gene_type="lncRNA", chromosome=krom)
    baslar = [g["start"] for g in kodlayan]
    en_uzun = max((g["end"] - g["start"] for g in kodlayan), default=0)

    krom_anti = 0
    for lnc in lncler:
        ortusen = 0
        # Başlangıcı lnc'nin bitişini geçmeyen genlerden, uzaklığı en uzun gen boyunu aşanlar örtüşemez.
        i = bisect_right(baslar, lnc["end"]) - 1
        while i >= 0 and baslar[i] >= lnc["start"] - en_uzun:
            g = kodlayan[i]
            if g["end"] >= lnc["start"] and g["strand"] != lnc["strand"]:
                ortusen += 1
                satirlar.append({"tur": "antisens", "gen1": g["gene_name"], "gen2": lnc["gene_name"],
                                 "kromozom": krom, "konum": lnc["start"]})
            i -= 1
        krom_anti += ortusen > 0     # bir lncRNA birkaç genle örtüşse de bir kez sayılır

    krom_cift = 0
    tssler = sorted(kodlayan, key=tss)
    for a, b in zip(tssler, tssler[1:]):
        # Sol gen - zincirde, sağ gen + zincirde ise iki gen zıt yöne ayrılır.
        if a["strand"] == "-" and b["strand"] == "+" and tss(b) - tss(a) <= MESAFE:
            krom_cift += 1
            satirlar.append({"tur": "cift_yonlu", "gen1": a["gene_name"], "gen2": b["gene_name"],
                             "kromozom": krom, "konum": tss(a)})
    antisens += krom_anti
    cift_yonlu += krom_cift
    print(f"{krom:<7}{len(kodlayan):>9,}{len(lncler):>9,}{krom_anti:>10,}{krom_cift:>12,}")

print(f"\nToplam antisens lncRNA: {antisens:,}, çift yönlü protein kodlayan çift: {cift_yonlu:,}")
print("Örnek çift yönlü:", ", ".join(f"{s['gen1']}/{s['gen2']}" for s in satirlar if s["tur"] == "cift_yonlu")[:90])
print("TP53 ile antisens örtüşen lncRNA:", [s["gen2"] for s in satirlar if s["tur"] == "antisens" and s["gen1"] == "TP53"])

with open(CIKTI / "mimari_ciftler.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)
print("Çıktılar:", CIKTI)
