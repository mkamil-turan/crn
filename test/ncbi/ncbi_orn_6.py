"""Örnek 6: Paketleri birleştir, gen ekle ve çıkar.

Senaryo: Kanser panelinize Lynch sendromu genlerini eklemek istiyorsunuz. Bazı
genler zaten pakette varsa yeniden indirilmemeli; bir gen de panelden çıkacak.

Bu örnekte:
  1. Mevcut paket okunup ada göre sözlüğe alınır.
  2. Yalnızca pakette olmayan genler indirilir.
  3. Çıkarılacak gen silinir, liste konuma göre sıralanır.
  4. Yeni paket yazılır, paylaşılmak üzere bir gen ayrı dosyaya alınır.

Önce ncbi_orn_4.py çalıştırılmalıdır.
Çıktılar: test/ncbi/cikti/ncbi_orn_6/
"""
from pathlib import Path

from crn import read_crn_list, write_crn_list
from crn.fetch import Ncbi

KOK = Path(__file__).parent / "cikti"
GIRDI = KOK / "ncbi_orn_4" / "kanser_paneli.crn"
CIKTI = KOK / "ncbi_orn_6"
CIKTI.mkdir(parents=True, exist_ok=True)

if not GIRDI.is_file():
    raise SystemExit("Paket yok: önce test/ncbi/ncbi_orn_4.py dosyasını çalıştırın.")

EKLENECEK = ["MLH1", "MSH2", "MSH6", "PMS2", "TP53"]   # TP53 pakette zaten var
CIKARILACAK = ["PALB2"]

paket = {gen.name: gen for gen in read_crn_list(GIRDI)}
print("Mevcut paket:", ", ".join(paket))

ncbi = Ncbi(verbose=False)
eksik = [sembol for sembol in EKLENECEK if sembol not in paket]
print("Atlanan (zaten var):", ", ".join(s for s in EKLENECEK if s in paket))
print("İndirilen:", ", ".join(eksik))
for sembol in eksik:
    paket[sembol] = ncbi.gene(sembol)

for sembol in CIKARILACAK:
    print("Çıkarılan:", paket.pop(sembol).name)


def konum_sirasi(gen) -> tuple:
    kromozom = gen.features.first("chromosome")
    return (int(kromozom) if kromozom.isdigit() else 99, gen.features.first("genomic_region")[1])


genler = sorted(paket.values(), key=konum_sirasi)
yol = CIKTI / "genisletilmis_panel.crn"
write_crn_list(genler, yol)

print(f"\nYeni paket ({len(genler)} gen, {yol.stat().st_size / 1e6:.1f} MB):")
for gen in read_crn_list(yol):
    print(f"  chr{gen.features.first('chromosome'):<3} {gen}")

# Tek bir geni başkasıyla paylaşmak için paketten ayırmak yeterlidir.
paket["MLH1"].to_crn(CIKTI / "mlh1.crn")
print("\nÇıktılar:", CIKTI)
