"""Grup 3: Hastalık mekanizmasına göre gen grupları ve kodon kullanımı.

Senaryo: Patojen varyantları çoğunlukla protein kesen genler (işlev kaybı,
haployetmezlik) ile çoğunlukla yanlış anlamlı olan genler (baskın negatif ya da
işlev kazanımı) arasında kodon kullanımı farkı olup olmadığını soruyorsunuz.
İşlev kaybına duyarlı genler yüksek ifadeli ve evrimsel olarak korunmuş olma
eğilimindedir, bu da çeviri seçilimiyle ilişkili bir kodon imzası bırakabilir.

Bu örnekte:
  1. Gen başına kesen ve yanlış anlamlı P/LP varyant sayıları genom ölçeğinde alınır.
  2. En az 20 varyantlı genlerden kesen payı ≥ %90 ve ≤ %10 olanlar iki grup olur;
     her gruptan en çok kanıtlı 80 gen seçilir.
  3. Gruplar gen listesiyle tanımlanıp tek CRN dosyasına paketlenir.
  4. Diferansiyel kodon kullanımı testi yapılır, Türkçe yorum yazılır.

Çıktılar: test/clinvar/cikti/clinvar_grup_3/
Not: Gruplar uzunluk ve GC içeriğine göre eşleştirilmemiştir; kesen varyant
toplayan genler ortalamada daha uzundur. Bulunan kodonlar GC farkını yansıtabilir.
"""
from pathlib import Path

from crn.cub import DifferentialCodonUsage
from crn.fetch import Clinvar, Gencode

CIKTI = Path(__file__).parent / "cikti" / "clinvar_grup_3"
CIKTI.mkdir(parents=True, exist_ok=True)

PLP = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 1}
BOYUT = 80
cv = Clinvar(verbose=False)
gc = Gencode(verbose=False)

kesen = {g["gene_symbol"]: g["variants"] for g in cv.genes(effect=["nonsense", "frameshift"], **PLP)}
yanlis = {g["gene_symbol"]: g["variants"] for g in cv.genes(effect="missense", **PLP)}
payi = {}
for gen in set(kesen) | set(yanlis):
    k, y = kesen.get(gen, 0), yanlis.get(gen, 0)
    if k + y >= 20:
        payi[gen] = (k / (k + y), k + y)

lof = sorted((g for g, (p, _) in payi.items() if p >= 0.9), key=lambda g: -payi[g][1])[:BOYUT]
mis = sorted((g for g, (p, _) in payi.items() if p <= 0.1), key=lambda g: -payi[g][1])[:BOYUT]
print(f"Kesen baskın: {len(lof)} gen, ör. {', '.join(lof[:8])}")
print(f"Yanlış anlamlı baskın: {len(mis)} gen, ör. {', '.join(mis[:8])}\n")

paket = CIKTI / "mekanizma.crn"
rapor = cv.pack({"islev_kaybi": {"gene": lof}, "yanlis_anlamli": {"gene": mis}},
                paket, source=gc, **PLP)
print(f"Paket {paket.stat().st_size / 1e6:.0f} MB, gruplar {rapor['groups']}, "
      f"atlanan {len(rapor['skipped'])}: {', '.join(sorted(rapor['skipped'])[:5])}")

dcu = DifferentialCodonUsage.from_crn(paket, "islev_kaybi", "yanlis_anlamli")
print(dcu)
anlamli = dcu.significant("family")
saglam = dcu.robust("family", n=2000)
print(f"Aile düzeyi Fisher q<0.05: {len(anlamli)} kodon, permütasyondan da geçen: {len(saglam)}\n")
for r in saglam[:8]:
    print(f"  {r['codon']} ({r['amino_acid']})  log2 OR {r['log2_or']:+.2f}  -> {r['favored']}")

Clinvar.write_tsv(saglam or anlamli, CIKTI / "kodonlar.tsv")
yorum = dcu.interpret(permutations=2000)
(CIKTI / "yorum.txt").write_text(yorum + "\n")
print("\n" + next(satir for satir in yorum.splitlines() if satir.startswith("Bileşim")))
print("Sağlam kodonların neredeyse hepsi A/U ile biten kodonu işlev kaybı grubunda gösteriyorsa\n"
      "fark büyük ölçüde GC3 farkıdır; ayrıştırmak için grupları GC3'e göre eşleştirin.")
print(f"\nYorum ve tablo: {CIKTI}")
