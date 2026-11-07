"""Örnek 12: Benzer gen ailesinde özgül qPCR primeri tasarla.

Senaryo: HBB ifadesini ölçmek istiyorsunuz ama HBD, HBG1, HBG2 ve HBE1 çok benzer.
Primerler bu genlerin transkriptlerine yapışmamalı.

Bu örnekte:
  1. Globin genleri indirilip tek pakete yazılır.
  2. HBB mRNA'sı RNA -> cDNA dönüşümüyle DNA'ya çevrilir.
  3. sliding ile aday primerler çıkarılır ve kurallarla süzülür.
  4. Her aday, primerin 3' ucu üzerinden diğer genlerin transkriptlerinde aranır.
  5. En iyi çift seçilir, in silico PCR yapılır ve sonuç kaydedilir.

Çıktılar: test/ncbi/cikti/ncbi_orn_12/
Not: Yalnızca bu beş gene karşı tam eşleşme aranır, yani bu bir ön elemedir.
Sipariş öncesinde Primer-BLAST ile tüm transkriptoma karşı denetleyin. Tm yaklaşıktır.
"""
from pathlib import Path

from crn import DNA, read_crn_list, write_crn_list, write_fasta
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_12"
CIKTI.mkdir(parents=True, exist_ok=True)

AILE = ["HBB", "HBD", "HBG1", "HBG2", "HBE1"]
BOY = 20
UC = 12                        # özgüllükte bakılan 3' uç uzunluğu
GC_ARALIGI = (40.0, 60.0)
TM_ARALIGI = (52.0, 60.0)
URUN_ARALIGI = (80, 200)

yol = CIKTI / "globin_ailesi.crn"
if not yol.is_file():
    ncbi = Ncbi(verbose=False)
    write_crn_list([ncbi.gene(s) for s in AILE], yol)
genler = {g.name: g for g in read_crn_list(yol)}


def sens_dna(rna) -> DNA:
    """mRNA ile aynı metni taşıyan 5'->3' DNA (cDNA'nın ikinci zinciri)."""
    return rna.reverse_transcription().complementary()


hedef = sens_dna(next(r for k, r in genler["HBB"].transcripts.items() if k.startswith("NM_")))
digerleri = [sens_dna(r) for ad, g in genler.items() if ad != "HBB" for r in g.transcripts.values()]
print(f"Hedef: {hedef}\nKarşılaştırılan transkript: {len(digerleri)}")


def gc(dizi: str) -> float:
    return 100 * (dizi.count("G") + dizi.count("C")) / len(dizi)


def tm(dizi: str) -> float:
    return 64.9 + 41 * (dizi.count("G") + dizi.count("C") - 16.4) / len(dizi)


def uygun(p: DNA) -> bool:
    d = p.seq
    return (GC_ARALIGI[0] <= gc(d) <= GC_ARALIGI[1] and TM_ARALIGI[0] <= tm(d) <= TM_ARALIGI[1]
            and d[-1] in "GC" and not any(b * 4 in d for b in "ACGT"))


def ozgul(p: DNA) -> bool:
    uc = p[-UC:]
    return not any(uc in t or uc.reverse_complement() in t for t in digerleri)


ileri_tum = [(i, p) for i, p in enumerate(hedef.sliding(BOY)) if uygun(p)]
geri_tum = [(len(hedef) - i, p) for i, p in enumerate(hedef.reverse_complement().sliding(BOY)) if uygun(p)]
ileri = [(i, p) for i, p in ileri_tum if ozgul(p)]
geri = [(i, p) for i, p in geri_tum if ozgul(p)]
print(f"Kurallara uyan -> özgül: ileri {len(ileri_tum)} -> {len(ileri)}, geri {len(geri_tum)} -> {len(geri)}")

ciftler = [(abs(tm(f.seq) - tm(r.seq)), son - bas, bas, son, f, r)
           for bas, f in ileri for son, r in geri if URUN_ARALIGI[0] <= son - bas <= URUN_ARALIGI[1]]
if not ciftler:
    raise SystemExit("Uygun çift yok: eşikleri gevşetin.")
fark, urun, bas, son, f, r = min(ciftler, key=lambda c: (c[0], c[1]))
f.name, r.name = "HBB_qPCR_F", "HBB_qPCR_R"
amplikon = hedef[bas:son]
amplikon.name = "HBB_qPCR_amplikon"
assert f in hedef and r.reverse_complement() in hedef

print(f"\nİleri: 5'-{f.seq}-3'  GC %{gc(f.seq):.0f}, Tm {tm(f.seq):.1f}°C")
print(f"Geri : 5'-{r.seq}-3'  GC %{gc(r.seq):.0f}, Tm {tm(r.seq):.1f}°C")
print(f"Ürün : {urun} bp, transkriptte {bas + 1}-{son}")

write_crn_list([f, r, amplikon], CIKTI / "hbb_qpcr.crn")
write_fasta([f, r, amplikon], CIKTI / "hbb_qpcr.fasta")
print("Çıktılar:", CIKTI)
