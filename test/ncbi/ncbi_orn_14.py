"""Örnek 14: Pakete girmeden önce kayıtları denetle.

Senaryo: Yalnızca tutarlı transkriptlerden oluşan bir paket istiyorsunuz. Her
NM_ kaydında CDS'nin çevirisi NCBI proteiniyle aynı olmalı. Selenoproteinlerde
UGA kodonu selenosistein (U) kodladığı için standart çeviri burada yanılır.

Bu örnekte:
  1. Selenoproteinleri de içeren bir gen listesi indirilir.
  2. Her NM_ için AUG, stop, 3'e bölünme ve protein eşleşmesi denetlenir.
  3. Selenosistein kaynaklı farklar ayrı sınıfta raporlanır.
  4. Geçen kayıtlardan yeni Gene nesneleri kurulup tek pakete yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_14/
"""
import csv
from pathlib import Path

from crn import RNA, Gene, PRT, write_crn_list
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_14"
CIKTI.mkdir(parents=True, exist_ok=True)

GENLER = ["HBB", "INS", "MEFV", "GPX1", "SELENOP", "TXNRD1"]


def denetle(rna: RNA, protein: PRT) -> str:
    bas, bit = rna.features.first("cds")
    cds = rna[bas - 1:bit]
    if len(cds) % 3:
        return "3'e bölünmüyor"
    if cds.seq[:3] != "AUG":
        return f"AUG değil ({cds.seq[:3]})"
    ceviri = cds.translation().seq
    if not ceviri.endswith("*"):
        return "stop yok"
    ceviri = ceviri[:-1]
    if ceviri == protein.seq:
        return "TAMAM"
    # Selenosistein: çeviride * olan yerler NCBI proteininde U'dur.
    if len(ceviri) == len(protein.seq) and all(
            a == b or (a == "*" and b == "U") for a, b in zip(ceviri, protein.seq)):
        return "selenoprotein"
    return "protein farklı"


ncbi = Ncbi(verbose=False)
rapor, temiz = [], []
for sembol in GENLER:
    gen = ncbi.gene(sembol)
    yeni = Gene(gen.name, gen.dna)
    yeni.features = gen.features.copy()
    for tid, rna in gen.transcripts.items():
        if not tid.startswith("NM_") or "cds" not in rna.features:
            continue
        protein = gen.protein(rna.features.first("protein_id"))
        sonuc = denetle(rna, protein)
        rapor.append({"gen": gen.name, "transkript": tid, "protein": protein.name,
                      "u_sayisi": protein.seq.count("U"), "sonuc": sonuc})
        if sonuc in ("TAMAM", "selenoprotein"):
            yeni.add_transcript(rna, tid)
            yeni.add_protein(protein, protein.name)
    temiz.append(yeni)

print(f"{'Gen':<9}{'Transkript':<16}{'Protein':<16}{'Sec':>4}  Sonuç")
for s in rapor:
    print(f"{s['gen']:<9}{s['transkript']:<16}{s['protein']:<16}{s['u_sayisi']:>4}  {s['sonuc']}")

with open(CIKTI / "denetim_raporu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(rapor[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(rapor)

yol = CIKTI / "denetlenmis_paket.crn"
write_crn_list(temiz, yol)
print("\nPaket:", ", ".join(str(g) for g in temiz))
print("Çıktılar:", CIKTI)
