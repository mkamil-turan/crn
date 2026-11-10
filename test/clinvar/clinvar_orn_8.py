"""Örnek 8: Patojen eşanlamlı varyantlar ve kodon kullanımı.

Senaryo: Amino asidi değiştirmediği halde patojen sayılan "sessiz" varyantlarla
ilgileniyorsunuz. Bunların çoğu ekleme bölgesini bozar, ama bir kısmının etkisi
kodon kullanımı, mRNA yapısı ya da çeviri hızı üzerinden olabilir. Her varyantın
hangi kodonu hangisine çevirdiğini ve yeni kodonun daha nadir olup olmadığını
görmek, aday listesini daraltmanın ilk adımıdır.

Bu örnekte:
  1. P/LP, en az 1 yıldızlı eşanlamlı tek nükleotit varyantları alınır.
  2. Genler GENCODE ana transkriptiyle kurulur; ClinVar'daki referans baz CDS ile
     denetlenir (MANE transkriptleri aynı diziyi taşır).
  3. Eski ve yeni kodon, değişimin kodondaki ve CDS içindeki konumu, iki kodonun
     RSCU değeri hesaplanır.
  4. Yeni kodonun daha nadir olduğu varyantlar sayılır ve bu sonucun mutasyon
     yanlılığıyla nasıl karıştığı gösterilir; tablo TSV olarak yazılır.

Çıktılar: test/clinvar/cikti/clinvar_orn_8/eşanlamlı.tsv
Not: RSCU bu genlerin birleşik kodon kullanımından hesaplanır, genom ortalaması
değildir. Ekzon-intron sınırına yakın eşanlamlı varyantlar ekleme etkisi açısından
ayrıca SpliceAI gibi araçlarla değerlendirilmelidir. GENCODE verisi gerekir.
"""
import re
import warnings
from collections import Counter
from pathlib import Path

from crn.cub import CUB, Synonymous
from crn.fetch import Clinvar, Gencode

CIKTI = Path(__file__).parent / "cikti" / "clinvar_orn_8"
CIKTI.mkdir(parents=True, exist_ok=True)

cv = Clinvar(verbose=False)
gc = Gencode(verbose=False)
SNV = re.compile(r"c\.(\d+)([ACGT])>([ACGT])")


def sembol(satir):
    """Varyantın genini verir; adlandırılmamış LOC lokusları (ABCA4;LOC126805793) atlanır."""
    adli = [s for s in satir["gene_symbol"].split(";") if not re.fullmatch(r"LOC\d+", s)]
    return adli[0] if adli else ""


varyantlar = [r for r in cv.variants(effect="synonymous", min_stars=1,
                                     significance=["pathogenic", "likely_pathogenic"])
              if SNV.fullmatch(r["hgvs_c"] or "")]
semboller = sorted({sembol(r) for r in varyantlar} - {""})
print(f"{len(varyantlar)} eşanlamlı SNV, {len(semboller)} gen")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")              # kurulamayan genler aşağıda sayılır
    genler = cv.to_genes(gene=semboller, source=gc)
adlar = {g.features.first("gene_name"): g.name for g in genler}
cub = CUB(genler)
rscu = Synonymous(cub).rscu()                    # bu genlerin birleşik kodon kullanımı
print(f"GENCODE ile kurulan gen: {len(genler)}, kurulamayan: {len(semboller) - len(genler)}\n")

tablo, uyumsuz = [], 0
for r in varyantlar:
    if sembol(r) not in adlar:
        continue
    cds = cub.cds(adlar[sembol(r)]).replace("U", "T")
    konum, ref, alt = SNV.fullmatch(r["hgvs_c"]).groups()
    konum = int(konum)
    if konum > len(cds) or cds[konum - 1] != ref:
        uyumsuz += 1                             # ClinVar transkripti ana transkript değil
        continue
    bas = (konum - 1) // 3 * 3
    eski = cds[bas:bas + 3]
    yeni = eski[:(konum - 1) % 3] + alt + eski[(konum - 1) % 3 + 1:]
    tablo.append({"gene": sembol(r), "hgvs_c": r["hgvs_c"], "hgvs_p": r["hgvs_p"],
                  "codon_position": (konum - 1) % 3 + 1, "cds_fraction": round(konum / len(cds), 3),
                  "old_codon": eski, "new_codon": yeni,
                  "old_rscu": round(rscu[eski.replace("T", "U")], 3),
                  "new_rscu": round(rscu[yeni.replace("T", "U")], 3),
                  "stars": r["stars"], "conditions": r["conditions"]})

print(f"CDS ile eşleşen: {len(tablo)}, referans bazı tutmayan: {uyumsuz}")
print("Kodondaki konum:", dict(sorted(Counter(t["codon_position"] for t in tablo).items())))
nadir = [t for t in tablo if t["new_rscu"] < t["old_rscu"]]
print(f"Yeni kodonu daha nadir olan: {len(nadir)} / {len(tablo)}")
gc_at = [x for x in tablo if x["codon_position"] == 3 and x["old_codon"][2] in "GC"
         and x["new_codon"][2] in "AT"]
print(f"Üçüncü konumda G/C -> A/T: {len(gc_at)} / {len(tablo)}; bunların "
      f"{sum(x['new_rscu'] < x['old_rscu'] for x in gc_at)} tanesinde yeni kodon daha nadir.")
print("İnsanda G/C ile biten kodonlar daha sık kullanıldığından ve G>A geçişi en sık\n"
      "mutasyon olduğundan 'yeni kodon daha nadir' sonucu mutasyon yanlılığından beklenen\n"
      "yöndür; tek başına çeviri etkisine kanıt sayılmaz.")
print("En sık kodon değişimleri:",
      ", ".join(f"{a}>{b} ({n})" for (a, b), n in
                Counter((t["old_codon"], t["new_codon"]) for t in tablo).most_common(6)))

print("\nRSCU'su en çok düşen beş varyant:")
for t in sorted(tablo, key=lambda t: t["new_rscu"] - t["old_rscu"])[:5]:
    print(f"  {t['gene']:<7} {t['hgvs_c']:<12} p.{t['hgvs_p']:<9} {t['old_codon']}>{t['new_codon']}  "
          f"RSCU {t['old_rscu']:.2f} -> {t['new_rscu']:.2f}")

cv.write_tsv(tablo, CIKTI / "eşanlamlı.tsv")
print(f"\nYazıldı: {CIKTI / 'eşanlamlı.tsv'}")
