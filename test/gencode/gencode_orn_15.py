"""Örnek 15: Ad kalıbıyla gen ailesi sorgula ve hatalı listeyi yönet.

Senaryo: HOX genlerini kümeler halinde toplamak, ayrıca elle yazılmış bir listede
bulunmayan ya da birden çok gene karşılık gelen adları raporlamak istiyorsunuz.

Bu örnekte:
  1. Düzenli ifadeyle HOX genleri sorgulanır, kümeye ve sıraya göre gruplanır.
  2. Karışık bir gen listesi denenir; hatalar ayrı ayrı sınıflanır.
  3. Belirsiz adlar için aday gene_id'ler gösterilir.
  4. HOX genleri tek CRN listesine yazılır, rapor TSV'ye kaydedilir.

Çıktılar: test/gencode/cikti/gencode_orn_15/
"""
import csv
import re
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_15"
CIKTI.mkdir(parents=True, exist_ok=True)

gc = Gencode(verbose=False)

hox = gc.genes(gene_type="protein_coding", name=r"^HOX[A-D]\d+$")
kumeler: dict = {}
for g in sorted(hox, key=lambda g: g["start"]):
    kumeler.setdefault(g["gene_name"][:4], []).append(g)
print(f"HOX geni: {len(hox)}\n")
for kume, liste in sorted(kumeler.items()):
    zincir = "".join(sorted({g["strand"] for g in liste}))
    adlar = " ".join(re.sub(r"^HOX[A-D]", "", g["gene_name"]) for g in liste)
    genislik = (liste[-1]["end"] - liste[0]["start"]) / 1000
    print(f"{kume} {liste[0]['chromosome']:<6}{zincir:<3}{genislik:>6.0f} kb  paraloglar: {adlar}")

paket = [gc.gene(g["gene_id"]) for g in hox]
yol = CIKTI / "hox_genleri.crn"
write_crn_list(paket, yol)
print(f"\nPaket: {len(read_crn_list(yol))} gen, {yol.stat().st_size / 1e6:.1f} MB")

print("\nKarışık liste:")
rapor = []
for ad in ["TP53", "PTENP1", "BU_GEN_YOK", "ENSG00000139618", "U6", "XIST"]:
    try:
        kayit = gc.genes(name=f"^{re.escape(ad)}$") if not ad.startswith("ENSG") else []
        gen = gc.gene(ad)
    except KeyError:
        durum, ayrinti = "BULUNAMADI", ""
    except ValueError:
        durum, ayrinti = "BELİRSİZ", f"{len(kayit)} aday: " + ", ".join(g["gene_id"] for g in kayit[:3])
    else:
        durum, ayrinti = "TAMAM", f"{gen.name} ({gen.features.first('gene_type')})"
    rapor.append({"ad": ad, "durum": durum, "ayrinti": ayrinti})
    print(f"  {ad:<18}{durum:<12}{ayrinti}")

with open(CIKTI / "sorgu_raporu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(rapor[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(rapor)
print("\nÇıktılar:", CIKTI)
