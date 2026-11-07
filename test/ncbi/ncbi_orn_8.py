"""Örnek 8: Gen paketinden DNA, RNA ve protein paketleri üret.

Senaryo: Analiz araçlarınız yalnızca CDS ya da yalnızca protein dizisi istiyor.
Her gen için ana transkripti seçip molekül tipine göre ayrı paketler yapmak ve
başka programlar için FASTA dışa aktarmak istiyorsunuz.

Bu örnekte:
  1. Örnek 4'ün paketi okunur.
  2. Her gen için CDS'si en uzun NM_ transkripti seçilir.
  3. CDS (RNA), protein (PRT) ve genomik (DNA) paketleri ayrı yazılır.
  4. Paketler okunup sınıfları doğrulanır, FASTA çıktıları üretilir.

Önce ncbi_orn_4.py çalıştırılmalıdır.
Çıktılar: test/ncbi/cikti/ncbi_orn_8/
Not: En uzun CDS her zaman klinik referans transkript (MANE Select) değildir.
"""
from collections import Counter
from pathlib import Path

from crn import RNA, Gene, read_crn_list, write_crn_list, write_fasta

KOK = Path(__file__).parent / "cikti"
GIRDI = KOK / "ncbi_orn_4" / "kanser_paneli.crn"
CIKTI = KOK / "ncbi_orn_8"
CIKTI.mkdir(parents=True, exist_ok=True)

if not GIRDI.is_file():
    raise SystemExit("Paket yok: önce test/ncbi/ncbi_orn_4.py dosyasını çalıştırın.")


def ana_transkript(gen: Gene) -> RNA:
    """CDS'si en uzun NM_ transkripti."""
    adaylar = [r for k, r in gen.transcripts.items() if k.startswith("NM_") and "cds" in r.features]
    if not adaylar:
        raise ValueError(f"{gen.name}: protein kodlayan NM_ kaydı yok")
    return max(adaylar, key=lambda r: r.features.first("cds")[1] - r.features.first("cds")[0])


cds_listesi, protein_listesi, dna_listesi = [], [], []
print(f"{'Gen':<8}{'Transkript':<16}{'CDS nt':>8}{'aa':>7}")
for gen in read_crn_list(GIRDI):
    rna = ana_transkript(gen)
    bas, bit = rna.features.first("cds")

    cds = rna[bas - 1:bit]
    cds.name = f"{gen.name}_CDS"
    cds_listesi.append(cds)

    protein = gen.protein(rna.features.first("protein_id"))
    protein.name = f"{gen.name}_protein"
    protein_listesi.append(protein)

    dna = gen.dna[:]
    dna.name = f"{gen.name}_genomik"
    dna_listesi.append(dna)
    print(f"{gen.name:<8}{rna.name:<16}{len(cds):>8,}{len(protein):>7,}")

for ad, liste in (("panel_cds", cds_listesi), ("panel_protein", protein_listesi),
                  ("panel_genomik", dna_listesi)):
    yol = CIKTI / f"{ad}.crn"
    write_crn_list(liste, yol)
    tipler = Counter(type(x).__name__ for x in read_crn_list(yol))
    print(f"{yol.name:<20}{dict(tipler)}  {yol.stat().st_size / 1e6:.2f} MB")

write_fasta(cds_listesi, CIKTI / "panel_cds.fasta")
write_fasta(protein_listesi, CIKTI / "panel_protein.fasta")
print("\nÇıktılar:", CIKTI)
