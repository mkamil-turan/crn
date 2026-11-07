"""Örnek 9: Proteinde transmembran bölge adaylarını tara.

Senaryo: Bir akuaporinin zara kaç kez girip çıktığını kabaca görmek ve iki NPA
motifini bulmak istiyorsunuz. Yöntem Kyte & Doolittle (1982) hidropati grafiğidir.

Bu örnekte:
  1. Protein dizisi FASTA'dan okunur (önce örnek FASTA yazılır).
  2. sliding ile 19 amino asitlik pencereler çıkarılır; elemanlar PRT'dir.
  3. Pencere ortalaması 1.6'yı geçen bölgeler TM adayı olarak birleştirilir.
  4. Profil terminalde çizilir, NPA motifleri aranır.
  5. Sonuçlar TSV ve FASTA olarak kaydedilir.

Çıktılar: test/sequence/cikti/seq_orn_9/
Not: Hidropati grafiği kaba bir taramadır. Akuaporin 6 TM sarmalı içerir, ancak
bu yöntem yakın sarmalları tek bölge gösterebilir; kesin sonuç için DeepTMHMM
gibi araçlar ya da yapı verisi kullanın. Dizi örnek amaçlıdır; gerçek analizde
UniProt/NCBI kaydını alın.
"""
import csv
from pathlib import Path

from crn import PRT, write_fasta
from crn.tables import HYDROPATHY_SCALES

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_9"
CIKTI.mkdir(parents=True, exist_ok=True)

PENCERE = 19
ESIK = 1.6
OLCEK = HYDROPATHY_SCALES["Kyte-Doolittle"]

# 1. Girdi dosyası ve okuma.
girdi = CIKTI / "aqp1.fasta"
PRT("MASEFKKKLFWRAVVAEFLATTLFVFISIGSALGFKYPVGNNQTAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLG"
    "LLLSCQISIFRALMYIIAQCVGAIVATAILSGITSSLTGNSLGRNDLADGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDL"
    "GGSAPLAIGLSVALGHLLAIDYTGCGINPARSFGSAVITHNFSNHWIFWVGPFIGGALAVLIYDFILAPRSSDLTDRVKVWT"
    "SGQVEEYDLDADDINSRVEMKPK", name="AQP1_insan").to_fasta(girdi)
protein = PRT.from_fasta_file(girdi)[0]
print(protein)


# 2. Pencere hidropatisi. Standart dışı harfler (X, *) hesaba katılmaz.
def hidropati(pencere: PRT) -> float:
    degerler = [OLCEK[aa] for aa in pencere if aa in OLCEK]
    return sum(degerler) / len(degerler) if degerler else 0.0


pencereler = protein.sliding(PENCERE)
profil = [hidropati(p) for p in pencereler]
print(f"{len(pencereler)} pencere, sınıf: {type(pencereler[0]).__name__}, "
      f"en yüksek skor: {max(profil):.2f}")

# 3. Eşiği geçen ve birbiriyle örtüşen pencereleri tek bölgede birleştir.
bolgeler = []
yari = PENCERE // 2
for i, skor in enumerate(profil):
    if skor >= ESIK:
        if bolgeler and i < bolgeler[-1][1] + PENCERE:
            bolgeler[-1][1] = i
        else:
            bolgeler.append([i, i])

# 4. Terminal grafiği: pozitif taraf hidrofobik.
print(f"\nHidropati profili ({PENCERE} aa pencere, * eşik {ESIK} üstü)")
for i in range(0, len(profil), 10):
    skor = profil[i]
    cubuk = "#" * max(0, int((skor + 1) * 10))
    print(f"{i + yari + 1:>4} {skor:>6.2f} {cubuk}{' *' if skor >= ESIK else ''}")

satirlar = []
parcalar = []
for no, (ilk, son) in enumerate(bolgeler, start=1):
    bas, bit = ilk + 1, son + PENCERE     # 1 tabanlı, dıştaki pencere sınırları
    parca = protein[bas - 1:bit]
    parca.name = f"TM_adayi_{no}"
    parcalar.append(parca)
    tepe = max(profil[ilk:son + 1])
    satirlar.append({"aday": parca.name, "baslangic": bas, "bitis": bit,
                     "tepe_hidropati": round(tepe, 2), "dizi": parca.seq})

print(f"\nTM adayı: {len(satirlar)}")
for s in satirlar:
    print(f"  {s['aday']}: {s['baslangic']}-{s['bitis']}, tepe {s['tepe_hidropati']}")

motifler = [i + 1 for i in range(len(protein) - 2) if protein.seq[i:i + 3] == "NPA"]
print("NPA motif konumları:", motifler)

# 5. Kaydetme.
if parcalar:
    write_fasta(parcalar, CIKTI / "tm_adaylari.fasta")
    with open(CIKTI / "tm_adaylari.tsv", "w", newline="", encoding="utf-8") as dosya:
        yazici = csv.DictWriter(dosya, fieldnames=list(satirlar[0]), delimiter="\t")
        yazici.writeheader()
        yazici.writerows(satirlar)
print("\nÇıktılar:", CIKTI)
