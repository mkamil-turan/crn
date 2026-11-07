"""Örnek 1: Yüksek ve düşük ifadeli genlerin diferansiyel kodon kullanımı.

Senaryo: Ribozomal proteinler gibi çok yüksek ifade edilen genler, koku reseptörü
gibi düşük ifadeli genlerden farklı eş anlamlı kodonlar mı seçiyor?

Bu örnekte:
  1. İki grubun genleri NCBI'dan alınır, her biri tek izoforma indirilir.
  2. Grup adı her genin features'ına yazılır.
  3. İki grup tek CRN dosyasına, tek gen listesi olarak paketlenir.
  4. Test dosyadan başlatılır; küresel, aile ve alt aile düzeyleri hesaplanır.
  5. Tablolar TSV'ye yazılır, yanardağ grafiği çizilir, yorum dosyaya kaydedilir.

Çıktılar: test/cub/cikti/diferansiyel_orn_1/
Not: İlk çalıştırma, genler sys/ içinde yoksa internet ister. 12'şer genlik gruplar
yöntemi göstermek içindir; yayın için daha büyük ve dengeli gruplar kullanın.
"""
import warnings
from pathlib import Path

from crn import write_crn_list
from crn.cub import DifferentialCodonUsage
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "diferansiyel_orn_1"
CIKTI.mkdir(parents=True, exist_ok=True)

GRUPLAR = {
    "yuksek": ["RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2", "EEF1A1", "RPL5", "RPL10",
               "RPS4X", "RPL13", "RPS8"],
    "dusuk": ["OR1A1", "OR2W1", "OR5A1", "OR7A5", "OR1D2", "OR2J2", "OR4D1", "OR6A2",
              "OR10A2", "OR51E2", "OR2T4", "OR1E1"],
}

ncbi = Ncbi(verbose=False)
genler = []
for grup, semboller in GRUPLAR.items():
    for sembol in semboller:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")      # izoform ayıklama bildirimleri
            izoform = next(g for g in ncbi.gene(sembol).isoforms()
                           if g.features.first("transcript_id").startswith("NM_"))
        izoform.features.set("group", grup)
        genler.append(izoform)

paket = CIKTI / "gruplar.crn"
write_crn_list(genler, paket)
print(f"Paket: {len(genler)} gen, {paket.stat().st_size / 1e6:.1f} MB\n")

dcu = DifferentialCodonUsage.from_crn(paket, "yuksek", "dusuk")
print(dcu)

saglam = {s["codon"] for s in dcu.robust("family")}
print(f"\nAile düzeyinde Fisher'ın bulduğu kodonlar (* gen düzeyinde permütasyondan da geçti):")
print(f"{'Kodon':<7}{'AA':<4}{'yuksek':>8}{'dusuk':>8}{'pay A':>8}{'pay B':>8}{'OR':>8}{'log2OR':>8}{'q':>11}  Lehine")
for s in dcu.significant("family")[:12]:
    print(f"{s['codon']:<7}{s['amino_acid']:<4}{s['a']:>8}{s['c']:>8}{s['share_a']:>8.2f}{s['share_b']:>8.2f}"
          f"{s['odds_ratio']:>8.2f}{s['log2_or']:>8.2f}{s['q']:>11.2e}  {s['favored']}"
          f"{' *' if s['codon'] in saglam else ''}")

print("\nDüzeylere göre anlamlı kodon sayısı:")
for duzey in ("global", "family", "subfamily"):
    satirlar = dcu.table(duzey, file_name=CIKTI / f"{duzey}.tsv")
    print(f"  {duzey:<11}{sum(s['significant'] for s in satirlar):>3} / {len(satirlar)}")

yorum = dcu.interpret()
(CIKTI / "yorum.txt").write_text(yorum + "\n", encoding="utf-8")
print("\nYorum\n-----")
print(yorum)

dcu.plot_volcano("family", file_name=CIKTI / "yanardag_aile.png")
print("\nÇıktılar:", CIKTI)
