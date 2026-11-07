"""Örnek 5: Paketi inceltmek: yalnızca doğrulanmış RefSeq kayıtlarını tut.

Senaryo: Örnek 4'teki pakette tahmini (XM_) transkriptler de var. Raporda yalnızca küratörlü NM_ ve NR_ kayıtlarını kullanmak istiyorsunuz.

Bu örnekte:
  1. Paket internet olmadan okunur.
  2. Tahmini transkriptler ve onlara bağlı proteinler çıkarılır.
  3. İnceltilmiş paket ayrı dosyaya yazılır, boyutlar karşılaştırılır.

Önce ncbi_orn_4.py çalıştırılmalıdır.
Çıktılar: test/ncbi/cikti/ncbi_orn_5/
"""
from pathlib import Path

from crn import read_crn_list, write_crn_list

KOK = Path(__file__).parent / "cikti"
GIRDI = KOK / "ncbi_orn_4" / "kanser_paneli.crn"
CIKTI = KOK / "ncbi_orn_5"
CIKTI.mkdir(parents=True, exist_ok=True)

if not GIRDI.is_file():
    raise SystemExit("Paket yok: önce test/ncbi/ncbi_orn_4.py dosyasını çalıştırın.")

genler = read_crn_list(GIRDI)
print(f"{'Gen':<8}{'Transkript':>12}{'Kalan':>8}{'Protein':>10}{'Kalan':>8}")
for gen in genler:
    once_t, once_p = len(gen.transcripts), len(gen.proteins)
    for tid in gen.transcripts:
        if not tid.startswith(("NM_", "NR_")):
            gen.remove_transcript(tid)
    for pid, protein in gen.proteins.items():
        if protein.features.first("transcript_id") not in gen.transcripts:
            gen.remove_protein(pid)
    print(f"{gen.name:<8}{once_t:>12}{len(gen.transcripts):>8}{once_p:>10}{len(gen.proteins):>8}")

yol = CIKTI / "kanser_paneli_nm.crn"
write_crn_list(genler, yol)
once, sonra = GIRDI.stat().st_size, yol.stat().st_size
print(f"\nBoyut: {once / 1e6:.1f} MB -> {sonra / 1e6:.1f} MB (%{100 * (1 - sonra / once):.0f} küçüldü)")
print("Çıktılar:", CIKTI)
