"""Örnek 16: NcbiSession ile API anahtarlı çoklu gen indirme.

Senaryo: Kardiyomiyopati paneli için çok sayıda geni NCBI'dan indireceksiniz.
API anahtarı hız sınırını saniyede 3 istekten 10 isteğe çıkarır. Tek bir oturum
bütün genler için ortak kullanılır.

Bu örnekte:
  1. NcbiSession API anahtarı, zaman aşımı ve yeniden deneme ayarıyla açılır.
  2. Aynı oturum Ncbi nesnesine verilir.
  3. Genler sırayla indirilir; bir genin hatası diğerlerini durdurmaz.
  4. Her gen için süre ve sonuç raporlanır.
  5. Başarılı genler tek CRN paketine yazılır.
  6. with bloğu bitince oturumun kapandığı gösterilir.

Çıktılar: test/ncbi/cikti/ncbi_orn_16/
Not: API_KEY aşağıda bir yer tutucudur. Gerçek anahtarınızla değiştirin ve kodu
paylaşacaksanız anahtarı ortam değişkeninden okuyun. Anahtar hesabınızdan alınır:
https://www.ncbi.nlm.nih.gov/account/settings/
"""
import csv
import time
from pathlib import Path

import requests

from crn import read_crn_list, write_crn_list
from crn.fetch import Ncbi, NcbiSession

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_16"
CIKTI.mkdir(parents=True, exist_ok=True)

API_KEY = "this_is_an_api_key"
GENLER = ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "LMNA", "TPM1", "ACTC1", "MYL2", "MYL3", "PLN"]

genler, rapor = [], []

with NcbiSession(api_key=API_KEY, timeout=60, max_retries=3) as oturum:
    print("Oturum:", oturum)
    ncbi = Ncbi(session=oturum, verbose=False)

    toplam_bas = time.perf_counter()
    for sembol in GENLER:
        bas = time.perf_counter()
        try:
            gen = ncbi.gene(sembol)
        except KeyError as hata:
            durum, ayrinti = "BULUNAMADI", str(hata)
        except requests.HTTPError as hata:
            durum, ayrinti = f"HTTP {hata.response.status_code}", " ".join(hata.response.text.split())[:60]
            if "API key invalid" in hata.response.text:
                # Anahtar her istekte reddedilir; kalan genleri denemenin anlamı yok.
                rapor.append({"gen": sembol, "durum": durum, "sure_sn": 0, "ayrinti": ayrinti})
                print(f"{sembol:<8}{durum:<12}NCBI anahtarı reddetti, indirme durduruldu.")
                break
        except requests.RequestException as hata:
            durum, ayrinti = "AĞ HATASI", type(hata).__name__
        else:
            genler.append(gen)
            durum, ayrinti = "TAMAM", str(gen)
        sure = time.perf_counter() - bas
        rapor.append({"gen": sembol, "durum": durum, "sure_sn": round(sure, 1), "ayrinti": ayrinti})
        print(f"{sembol:<8}{durum:<12}{sure:>6.1f} sn  {ayrinti}")
    print(f"\nToplam süre: {time.perf_counter() - toplam_bas:.1f} sn")

print("Oturum açık mı:", oturum.is_open())
print("Oturum:", oturum)

with open(CIKTI / "indirme_raporu.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.DictWriter(dosya, fieldnames=list(rapor[0]), delimiter="\t")
    yazici.writeheader()
    yazici.writerows(rapor)

if not genler:
    raise SystemExit("Hiçbir gen indirilemedi: API anahtarını ve bağlantıyı denetleyin.")

yol = CIKTI / "kardiyomiyopati_paneli.crn"
write_crn_list(genler, yol)
print(f"\n{len(genler)}/{len(GENLER)} gen pakete yazıldı, {yol.stat().st_size / 1e6:.1f} MB")
print("Paket:", [g.name for g in read_crn_list(yol)])
print("Çıktılar:", CIKTI)
