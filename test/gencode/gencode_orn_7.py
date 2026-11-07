"""Örnek 7: Transkript varyantları: bir gen kaç farklı ürün verir?

Senaryo: TP53'ün tüm varyantlarını kaliteleriyle görmek, hangi ekzonların
kurucu (hep var), hangilerinin alternatif olduğunu anlamak istiyorsunuz.

Bu örnekte:
  1. Her transkriptin tipi, TSL, ekzon sayısı ve protein uzunluğu çıkarılır.
  2. Benzersiz protein sayısı bulunur.
  3. Tam uzunluklu transkriptlerde ekzon kullanım oranları hesaplanır.
  4. Tipi protein kodlayan olan varyantlar tek CRN listesine yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_7/
"""
import csv
import sys
from collections import Counter
from pathlib import Path

from crn import write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_7"
CIKTI.mkdir(parents=True, exist_ok=True)

GEN = sys.argv[1] if len(sys.argv) > 1 else "TP53"

gc = Gencode(verbose=False)
gen = gc.gene(GEN)
print(gen)

satirlar = []
for t in gc.transcripts(gene=GEN):
    tid = t["transcript_id"]
    ayrinti = gc.details(tid)
    protein = next((p for p in gen.proteins.values() if p.features.first("transcript_id") == tid), None)
    satirlar.append({"ad": t["transcript_name"], "kimlik": tid, "tip": t["transcript_type"],
                     "tsl": ayrinti["support_level"] or "-", "duzey": ayrinti["level"],
                     "ekzon": len(ayrinti["exons"]), "nt": len(gen.transcript(tid)) if tid in gen.transcripts else 0,
                     "aa": len(protein) if protein else 0,
                     "mane": "MANE_Select" in ayrinti["tags"],
                     "tam": not any(e in ayrinti["tags"] for e in ("mRNA_start_NF", "mRNA_end_NF", "cds_start_NF", "cds_end_NF"))})

print("\nTipler:", ", ".join(f"{t} ({n})" for t, n in Counter(s["tip"] for s in satirlar).most_common()))
print("TSL   :", ", ".join(f"{t}: {n}" for t, n in sorted(Counter(s["tsl"] for s in satirlar).items())))

print(f"\n{'Ad':<11}{'Tip':<26}{'TSL':>4}{'Ekzon':>7}{'nt':>7}{'aa':>6}")
for s in sorted(satirlar, key=lambda x: (not x["mane"], {"1": 1, "2": 2, "3": 3, "4": 4, "5": 5}.get(x["tsl"], 9), -x["aa"]))[:12]:
    print(f"{s['ad']:<11}{s['tip']:<26}{s['tsl']:>4}{s['ekzon']:>7}{s['nt']:>7}{s['aa']:>6}{'  MANE' if s['mane'] else ''}")

benzersiz = {p.seq for p in gen.proteins.values()}
print(f"\n{len(gen.proteins)} protein kaydı, {len(benzersiz)} benzersiz protein dizisi")

tamlar = [s for s in satirlar if s["tam"]]
kullanim = Counter()
for s in tamlar:
    kullanim.update(set(gc.details(s["kimlik"])["exons"]))
print(f"\nTam uzunluklu transkript: {len(tamlar)}, benzersiz ekzon: {len(kullanim)}")
print("En sık kullanılan ekzonlar:")
for (bas, bit), sayi in kullanim.most_common(5):
    print(f"  {bas:,}-{bit:,}  {bit - bas + 1:>4} nt  %{100 * sayi / len(tamlar):.0f}")

with open(CIKTI / f"{GEN}_varyantlari.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(satirlar)

kodlayan = [gen.transcript(s["kimlik"]) for s in satirlar if s["tip"] == "protein_coding" and s["kimlik"] in gen.transcripts]
write_crn_list(kodlayan, CIKTI / f"{GEN}_kodlayan_varyantlar.crn")
print(f"\n{len(kodlayan)} kodlayan varyant yazıldı. Başka gen için: python gencode_orn_7.py BRCA1")
print("Çıktılar:", CIKTI)
