"""Grup 5: Uçtan uca çalışma tasarımı; iki hastalık grubu, kontrol ve yinelenebilirlik.

Senaryo: Nörogelişimsel bozukluk genlerini ve mitokondriyal hastalık genlerini,
ClinVar'da çok sayıda varyantı olduğu halde hiç patojen varyantı bildirilmemiş
genlerden oluşan bir kontrol grubuyla karşılaştıracaksınız. Makalenin yöntem
bölümünde grupların tam tanımını, veri sürümlerini ve dosya özetlerini vermeniz
gerekecek.

Bu örnekte:
  1. Hastalık grupları kimlik ya da ad parçası, P/LP ve yıldız filtresiyle tanımlanır.
  2. Kontrol, en az 100 varyantı olan ama P/LP varyantı olmayan genlerden sabit
     tohumla seçilir; aynı tohum aynı grubu verir. Kontrolün kendi filtresi vardır.
  3. Üç grup tek pakete yazılır; ortak ve atlanan genler raporlanır.
  4. Her grup kontrolle ayrı ayrı test edilir, sonuçlar tabloya dökülür.
  5. Grup tanımları, ClinVar ve GENCODE sürümleri, arşivin SHA-256 özeti ve tohum
     JSON dosyasına yazılır.

Çıktılar: test/clinvar/cikti/clinvar_grup_5/
Not: "Patojen varyantı bildirilmemiş" gen hastalıkla ilgisiz demek değildir; yalnızca
ClinVar'da bu kanıtın bulunmadığını söyler. Kontrol gen uzunluğu ve GC'ye göre
eşleştirilmemiştir; bulunan kodonlar GC farkını da yansıtabilir. mtDNA genleri pakete
girmez; onları ayrıca incelemek için CUB(..., genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL)
kullanılır.
"""
import json
import random
import warnings
from collections import Counter
from pathlib import Path

from crn import read_crn_list
from crn.cub import CubWarning, DifferentialCodonUsage
from crn.fetch import Clinvar, Gencode

CIKTI = Path(__file__).parent / "cikti" / "clinvar_grup_5"
CIKTI.mkdir(parents=True, exist_ok=True)

TOHUM = 7
KONTROL_BOYUTU = 150
PLP = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 2}
cv = Clinvar(verbose=False)
gc = Gencode(verbose=False)

hepsi = cv.genes()                                         # tek genli varyantlar, tüm sınıflar
adaylar = sorted(g["gene_symbol"] for g in hepsi
                 if g["variants"] >= 100 and g["pathogenic"] == 0 and g["ncbi_gene_id"])
kontrol = random.Random(TOHUM).sample(adaylar, KONTROL_BOYUTU)
print(f"Kontrol havuzu: {len(adaylar)} gen; seçilen {KONTROL_BOYUTU} (tohum {TOHUM})")

GRUPLAR = {
    "norogelisimsel": {"condition": ["MONDO:0700092", "intellectual disability"], **PLP},
    "mitokondriyal": {"condition": "mitochondrial", **PLP},   # adında "mitochondrial" geçen hastalıklar
    "kontrol": {"gene": kontrol},
}
for ad, tanim in GRUPLAR.items():
    print(f"  {ad:<15} filtreden geçen {len(cv.genes(**tanim))} gen")

paket = CIKTI / "calisma.crn"
with warnings.catch_warnings():
    warnings.simplefilter("ignore")                        # atlananlar raporda yazılıdır
    rapor = cv.pack(GRUPLAR, paket, source=gc)
print(f"\nPaket {paket.stat().st_size / 1e6:.0f} MB; gruplar {rapor['groups']}")
print(f"Çakıştığı için çıkarılan {len(rapor['overlap'])}, atlanan {len(rapor['skipped'])} gen")
mt = sorted(g for g in rapor["skipped"] if g.startswith("MT-"))
print(f"Atlananların {len(mt)} tanesi mitokondriyal DNA geni ({', '.join(mt[:4])} ...); "
      "tRNA genleri protein kodlamaz,\nprotein genleri omurgalı mitokondriyal koduyla çevrilir. "
      "Mitokondriyal grup yalnızca çekirdek genlerinden oluşur.")
for gen in sorted(set(rapor["skipped"]) - set(mt)):
    print(f"  atlanan {gen}: {rapor['skipped'][gen]}")

genler = read_crn_list(paket)
print("Dosyadan:", dict(Counter(g.features.first("group") for g in genler)))

sonuc = []
for grup in ("norogelisimsel", "mitokondriyal"):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CubWarning)        # üçüncü grup bu çiftte kullanılmaz
        dcu = DifferentialCodonUsage.from_crn(paket, grup, "kontrol")
    saglam = dcu.robust("family", n=2000)
    sonuc.append({"group": grup, "genes": len(dcu.a), "control": len(dcu.b),
                  "fisher_q05": len(dcu.significant("family")), "robust": len(saglam),
                  "robust_codons": ",".join(r["codon"] for r in saglam)})
    print(f"\n{grup} / kontrol: Fisher {sonuc[-1]['fisher_q05']} kodon, sağlam {len(saglam)}: "
          f"{sonuc[-1]['robust_codons'] or '-'}")
    (CIKTI / f"{grup}_yorum.txt").write_text(dcu.interpret(permutations=2000) + "\n")
Clinvar.write_tsv(sonuc, CIKTI / "sonuc.tsv")

manifest = json.loads((cv.release_dir / "manifest.json").read_text())
ayarlar = {
    "clinvar": {"release": cv.release, "assembly": cv.assembly, "sha256": manifest["sha256"],
                "source": manifest["source"]},
    "gencode": gc.release, "seed": TOHUM, "control_size": KONTROL_BOYUTU,
    "control_pool": "genes with >=100 single-gene ClinVar variants and no P/LP variant",
    "groups": {ad: {k: v for k, v in t.items() if k != "gene"} or {"gene": "kontrol listesi"}
               for ad, t in GRUPLAR.items()},
    "pack_report": {"groups": rapor["groups"], "overlap": rapor["overlap"],
                    "skipped": sorted(rapor["skipped"])},
}
(CIKTI / "ayarlar.json").write_text(json.dumps(ayarlar, ensure_ascii=False, indent=2) + "\n")
print(f"\nSonuçlar, yorumlar ve ayarlar.json: {CIKTI}")
