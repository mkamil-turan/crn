"""Örnek 7: Metin dosyasındaki gen listesinden hatalara dayanıklı paket üret.

Senaryo: Bir meslektaşınız elle yazılmış bir gen listesi gönderdi. Listede yorum
satırları, tekrarlar, küçük harfler ve yanlış semboller olabilir. Bir hatalı
sembol bütün işi durdurmamalı.

Bu örnekte:
  1. Örnek liste dosyası yazılır ve temizlenerek okunur.
  2. Her sembol ayrı ayrı denenir; hata alanlar not edilir.
  3. Başarılı genler tek pakete yazılır.
  4. Her sembolün sonucu TSV raporuna kaydedilir.

Çıktılar: test/ncbi/cikti/ncbi_orn_7/
"""
import csv
from pathlib import Path

import requests

from crn import write_crn_list
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_7"
CIKTI.mkdir(parents=True, exist_ok=True)

liste_yolu = CIKTI / "gen_listesi.txt"
liste_yolu.write_text("""# Globin ve hemoliz çalışması
HBB
hba1
HBA2
HBB
# yazım hatası:
HBBX
G6PD
BU_GEN_YOK
""", encoding="utf-8")

# Sırayı koruyarak tekrarları at; dict.fromkeys ekleme sırasını saklar.
semboller = list(dict.fromkeys(
    satir.strip().upper() for satir in liste_yolu.read_text(encoding="utf-8").splitlines()
    if satir.strip() and not satir.startswith("#")))
print("Temizlenmiş liste:", ", ".join(semboller))

ncbi = Ncbi(verbose=False)
genler, rapor = [], []
for sembol in semboller:
    try:
        gen = ncbi.gene(sembol)
    except KeyError as hata:
        rapor.append({"sembol": sembol, "durum": "BULUNAMADI", "ayrinti": str(hata).strip("'\"")})
    except requests.RequestException as hata:
        rapor.append({"sembol": sembol, "durum": "AĞ HATASI", "ayrinti": type(hata).__name__})
    else:
        genler.append(gen)
        rapor.append({"sembol": sembol, "durum": "TAMAM", "ayrinti": str(gen)})

for satir in rapor:
    print(f"{satir['sembol']:<12}{satir['durum']:<12}{satir['ayrinti']}")

with open(CIKTI / "indirme_raporu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(rapor[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(rapor)

write_crn_list(genler, CIKTI / "globin_paketi.crn")
print(f"\n{len(genler)}/{len(semboller)} gen pakete yazıldı.")
print("Çıktılar:", CIKTI)
