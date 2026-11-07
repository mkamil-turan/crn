"""Örnek 16: GENCODE ve NCBI aynı proteini mi veriyor?

Senaryo: Laboratuvarınızda iki kaynak birlikte kullanılıyor. MANE Select tanımı
gereği GENCODE ve RefSeq transkriptleri eşleşmelidir; bunu panelde sınamak istiyorsunuz.

Bu örnekte:
  1. Her gen için GENCODE MANE Select proteini alınır.
  2. NCBI'daki NM_ proteinleri arasında birebir aynı olan aranır.
  3. Eşleşen transkriptlerin mRNA dizileri karşılaştırılır.
  4. İki kaynağın kayıtları tek CRN listesine yazılır; kaynak features'tan okunur.

Çıktılar: test/gencode/cikti/gencode_orn_16/
Not: NCBI kayıtları sys/ içinde yoksa bu örnek internet ister.
"""
from collections import Counter
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Gencode, Ncbi

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_16"
CIKTI.mkdir(parents=True, exist_ok=True)

PANEL = ["MEFV", "HBB", "TP53", "BRCA1", "BRCA2", "PTEN", "ATM"]

gc = Gencode(verbose=False)
ncbi = Ncbi(verbose=False)

paket = []
print(f"{'Gen':<8}{'GENCODE':<20}{'NCBI':<16}{'aa':>6}  mRNA")
for ad in PANEL:
    tid = gc.transcripts(gene=ad, tag="MANE_Select")[0]["transcript_id"]
    g_rna, g_protein = gc.rna(tid), gc.protein(tid)
    ncbi_gen = ncbi.gene(ad)

    n_rna = None
    for k, rna in ncbi_gen.transcripts.items():
        pid = rna.features.first("protein_id")
        if k.startswith("NM_") and pid and ncbi_gen.protein(pid) == g_protein and rna.seq == g_rna.seq:
            n_rna = rna      # protein ve mRNA birebir aynı olan en güçlü eşleşme
            break
    if n_rna is None:
        for k, rna in ncbi_gen.transcripts.items():
            pid = rna.features.first("protein_id")
            if k.startswith("NM_") and pid and ncbi_gen.protein(pid) == g_protein:
                n_rna = rna
                break
    if n_rna is None:
        print(f"{ad:<8}{tid:<20}{'eşleşme yok':<16}{len(g_protein):>6}")
        continue
    durum = "birebir aynı" if n_rna.seq == g_rna.seq else f"farklı ({len(g_rna)} / {len(n_rna)} nt)"
    print(f"{ad:<8}{tid:<20}{n_rna.name:<16}{len(g_protein):>6}  {durum}")
    paket += [g_rna, n_rna]

yol = CIKTI / "iki_kaynak.crn"
write_crn_list(paket, yol)
print("\nPaket kaynakları:", dict(Counter(x.features.first("source") for x in read_crn_list(yol))))
print("Çıktılar:", CIKTI)
