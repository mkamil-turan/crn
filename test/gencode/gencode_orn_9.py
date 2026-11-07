"""Örnek 9: MANE Select ile temiz klinik panel paketi.

Senaryo: Klinik raporda her gen için tek ve standart bir transkript gerekir. "En
uzun CDS" seçimi bazen yanlış kaydı getirir; MANE Select bu sorunu çözer.

Bu örnekte:
  1. Paneldeki her gen için MANE Select ve en uzun CDS'li transkript bulunur.
  2. İkisinin farklı olduğu genler gösterilir.
  3. Her gen için yalnızca DNA, MANE transkripti ve proteinini içeren yeni Gene kurulur.
  4. Temiz genler tek CRN listesine paketlenir, boyut tam paketle karşılaştırılır.

Çıktılar: test/gencode/cikti/gencode_orn_9/
"""
from pathlib import Path

from crn import Gene, read_crn_list, write_crn_list
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_9"
CIKTI.mkdir(parents=True, exist_ok=True)

PANEL = ["BRCA1", "BRCA2", "TP53", "PTEN", "ATM", "CHEK2", "PALB2", "MLH1", "MSH2", "AQP1"]

gc = Gencode(verbose=False)


def cds_boyu(rna) -> int:
    bas, bit = rna.features.first("cds", (1, 0))
    return bit - bas + 1


tam_paket, temiz_paket = [], []
print(f"{'Gen':<8}{'Transkript':>11}  {'MANE':<12}{'CDS':>7}{'aa':>6}  En uzun CDS farklı mı")
for ad in PANEL:
    gen = gc.gene(ad)
    tam_paket.append(gen)
    tid = gc.transcripts(gene=ad, tag="MANE_Select")[0]["transcript_id"]
    mane = gen.transcript(tid)
    protein = next(p for p in gen.proteins.values() if p.features.first("transcript_id") == tid)

    en_uzun = max(gen.transcripts.values(), key=cds_boyu)
    fark = "aynı" if cds_boyu(en_uzun) == cds_boyu(mane) else f"{en_uzun.features.first('transcript_name')} ({cds_boyu(en_uzun):,})"
    print(f"{ad:<8}{len(gen.transcripts):>11}  {mane.features.first('transcript_name'):<12}"
          f"{cds_boyu(mane):>7,}{len(protein):>6,}  {fark}")

    temiz = Gene(gen.name, gen.dna)
    temiz.features = gen.features.copy()
    temiz.features.set("transcript_id", tid)
    temiz.add_transcript(mane, tid)
    temiz.add_protein(protein, protein.name)
    temiz_paket.append(temiz)

tam, sade = CIKTI / "panel_tum_transkriptler.crn", CIKTI / "panel_mane.crn"
write_crn_list(tam_paket, tam)
write_crn_list(temiz_paket, sade)
print(f"\nTam paket {tam.stat().st_size / 1e6:.1f} MB -> MANE paketi {sade.stat().st_size / 1e6:.1f} MB")
print("Okunan:", read_crn_list(sade)[0])
print("Çıktılar:", CIKTI)
