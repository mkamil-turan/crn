"""Örnek 9: Genlerin promotör bölgelerini çek ve CpG içeriğini tara.

Senaryo: Farklı zincirlerdeki genlerin promotörlerinde metilasyon çalışacaksınız.
Her genin 5' ucunun 2 kb yukarısını doğru yönde almanız gerekiyor.

Bu örnekte:
  1. Gen koordinatından zincire uygun yukarı bölge hesaplanır.
  2. region() ile dizi 5'->3' yönünde alınır.
  3. Bölgenin genle kesintisiz birleştiği NCBI'dan alınan bir bağlantı
     parçasıyla doğrulanır.
  4. sliding ile 200 nt pencerelerde CpG adası ölçütü sayılır.
  5. Promotörler (DNA) tek CRN listesine ve FASTA'ya yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_9/
Not: Gen başlangıcı en uç transkripte göredir; asıl TSS doku ve izoforma göre değişir.
"""
from pathlib import Path

from crn import DNA, write_crn_list, write_fasta
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_9"
CIKTI.mkdir(parents=True, exist_ok=True)

GENLER = ["HBB", "TP53", "BRCA2", "PTEN", "MEFV"]
UZUNLUK = 2000
PENCERE = 200
ADIM = 10
BAGLANTI = 50


def cpg_olcutu(pencere: DNA) -> bool:
    """Gardiner-Garden & Frommer ölçütü: GC > %50, gözlenen/beklenen CpG > 0.6."""
    dizi = pencere.seq
    c, g = dizi.count("C"), dizi.count("G")
    if not c or not g:
        return False
    return (c + g) / len(dizi) > 0.5 and dizi.count("CG") * len(dizi) / (c * g) > 0.6


ncbi = Ncbi(verbose=False)
promotorler = []
print(f"{'Gen':<7}{'Z':<3}{'Bölge':<36}{'GC%':>6}{'CpG%':>7}  Bağlantı")
for sembol in GENLER:
    gen = ncbi.gene(sembol)
    erisim, bas, bit, zincir = gen.features.first("genomic_region")

    # (-) zincirde genin 5' ucu koordinatı büyük olan uçtur.
    if zincir == "+":
        yukari = (bas - UZUNLUK, bas - 1)
        baglanti = (bas - BAGLANTI, bas + BAGLANTI - 1)
    else:
        yukari = (bit + 1, bit + UZUNLUK)
        baglanti = (bit - BAGLANTI + 1, bit + BAGLANTI)

    promotor = ncbi.region(erisim, *yukari, zincir, name=f"{gen.name}_promotor")
    promotor.features.set("gene_name", gen.name)

    # promotörün sonu + genin başı, referanstaki kesintisiz diziyle aynı olmalı.
    beklenen = ncbi.region(erisim, *baglanti, zincir)
    birlesik = promotor.seq[-BAGLANTI:] + gen.dna.seq[:BAGLANTI]

    pencereler = promotor.sliding(PENCERE, step_size=ADIM)
    cpg = 100 * sum(cpg_olcutu(p) for p in pencereler) / len(pencereler)
    gc = 100 * (promotor.seq.count("G") + promotor.seq.count("C")) / len(promotor)
    promotor.features.set("cpg_pencere_yuzde", round(cpg, 1))
    promotorler.append(promotor)
    print(f"{gen.name:<7}{zincir:<3}{erisim + ':' + str(yukari[0]) + '-' + str(yukari[1]):<36}"
          f"{gc:>6.1f}{cpg:>7.1f}  {'doğru' if beklenen == birlesik else 'HATALI'}")

write_crn_list(promotorler, CIKTI / "promotorler.crn")
write_fasta(promotorler, CIKTI / "promotorler.fasta")
print("\nCpG%: ölçütü geçen pencerelerin yüzdesi")
print("Çıktılar:", CIKTI)
