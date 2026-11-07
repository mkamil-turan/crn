"""Örnek 6: PCR primeri tasarla ve bilgisayarda dene.

Senaryo: HBB kodlayan dizisinde beta-talasemiye yol açan c.118C>T konumunu
(CAG kodonu) kapsayan bir ürün çoğaltıp Sanger ile okumak istiyorsunuz.

Bu örnekte:
  1. sliding ile bütün 20 nt aday primerler çıkarılır.
  2. Adaylar GC, Tm, 3' ucu ve tekrar kurallarıyla elenir.
  3. Ters primerler, reverse_complement ile karşı zincirden aranır.
  4. Tm'leri yakın ve hedefi kapsayan en iyi çift seçilir.
  5. Çift dizide denenir (in silico PCR), ürün ve primerler kaydedilir.

Çıktılar: test/sequence/cikti/seq_orn_6/
Not: Tm burada basit bir formülle yaklaşık hesaplanır. Sipariş öncesinde
genom özgüllüğünü Primer-BLAST ile, Tm'yi tuz koşullarınızla doğrulayın.
"""
import csv
from pathlib import Path

from crn import DNA, write_fasta

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_6"
CIKTI.mkdir(parents=True, exist_ok=True)

PRIMER_BOYU = 20
GC_ARALIGI = (40.0, 60.0)
TM_ARALIGI = (50.0, 62.0)
URUN_ARALIGI = (150, 350)
EN_FAZLA_TM_FARKI = 2.0
HEDEF = (118, 120)  # CAG kodonunun 1 tabanlı konumları

hbb = DNA(
    "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTT"
    "GGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGAT"
    "CTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTT"
    "AGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAG"
    "CTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGC"
    "AAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCAC"
    "AAGTATCACTAA",
    name="HBB_CDS",
)


def gc_yuzde(dizi: str) -> float:
    return 100 * (dizi.count("G") + dizi.count("C")) / len(dizi)


def tm(dizi: str) -> float:
    """Yaklaşık erime sıcaklığı: Tm = 64.9 + 41 * (G+C - 16.4) / N."""
    return 64.9 + 41 * (dizi.count("G") + dizi.count("C") - 16.4) / len(dizi)


def uygun_mu(primer: DNA) -> bool:
    dizi = primer.seq
    if not GC_ARALIGI[0] <= gc_yuzde(dizi) <= GC_ARALIGI[1]:
        return False
    if not TM_ARALIGI[0] <= tm(dizi) <= TM_ARALIGI[1]:
        return False
    if dizi[-1] not in "GC":                       # 3' ucu G/C ile bitsin
        return False
    if any(baz * 4 in dizi for baz in "ACGT"):     # AAAA gibi tekrarlar kaymaya yol açar
        return False
    uc = primer[-4:].reverse_complement()          # 3' ucu primerin kendisine yapışmasın
    return uc.seq not in dizi[:-4]


# 1-3. Her pencerenin başlangıcı, listedeki sırasıdır.
ileri = [(i, p) for i, p in enumerate(hbb.sliding(PRIMER_BOYU)) if uygun_mu(p)]
karsi_zincir = hbb.reverse_complement()
geri = [(len(hbb) - i, p) for i, p in enumerate(karsi_zincir.sliding(PRIMER_BOYU)) if uygun_mu(p)]
print(f"Uygun ileri primer: {len(ileri)}, uygun geri primer: {len(geri)}")

# 4. Çift seçimi. geri listesindeki konum, ürünün dış sınırıdır.
ciftler = []
for bas, f in ileri:
    for son, r in geri:
        urun = son - bas
        if not URUN_ARALIGI[0] <= urun <= URUN_ARALIGI[1]:
            continue
        if not (bas + PRIMER_BOYU < HEDEF[0] and HEDEF[1] < son - PRIMER_BOYU):
            continue
        fark = abs(tm(f.seq) - tm(r.seq))
        if fark <= EN_FAZLA_TM_FARKI:
            ciftler.append((fark, urun, bas, son, f, r))
ciftler.sort(key=lambda c: (c[0], c[1]))

if not ciftler:
    raise SystemExit("Uygun çift yok: eşikleri gevşetin.")
print(f"Uygun çift: {len(ciftler)}\n")

fark, urun, bas, son, f, r = ciftler[0]
f.name, r.name = "HBB_F", "HBB_R"

# 5. In silico PCR: ileri primer dizide, geri primerin tümleyeni dizinin ileri yönünde bulunmalı.
assert f in hbb and r.reverse_complement() in hbb
amplikon = hbb[bas:son]
amplikon.name = "HBB_amplikon"

print(f"İleri primer : 5'-{f.seq}-3'  GC %{gc_yuzde(f.seq):.0f}, Tm {tm(f.seq):.1f}°C, konum {bas + 1}")
print(f"Geri primer  : 5'-{r.seq}-3'  GC %{gc_yuzde(r.seq):.0f}, Tm {tm(r.seq):.1f}°C, konum {son}")
print(f"Ürün         : {urun} bp (konum {bas + 1}-{son})")
print(f"Hedef kodon  : {hbb[HEDEF[0] - 1:HEDEF[1]].seq}, üründe konum {HEDEF[0] - bas}")

write_fasta([f, r, amplikon], CIKTI / "primerler_ve_urun.fasta")
with open(CIKTI / "primer_ciftleri.tsv", "w", newline="", encoding="utf-8") as dosya:
    yazici = csv.writer(dosya, delimiter="\t")
    yazici.writerow(["ileri", "geri", "urun_bp", "tm_ileri", "tm_geri", "tm_farki"])
    for fark_i, urun_i, _, _, f_i, r_i in ciftler[:20]:
        yazici.writerow([f_i.seq, r_i.seq, urun_i, round(tm(f_i.seq), 1),
                         round(tm(r_i.seq), 1), round(fark_i, 1)])
print("\nÇıktılar:", CIKTI)
