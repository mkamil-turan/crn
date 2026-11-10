"""Grup 4: İki hastalık grubu, tek CRN dosyası, diferansiyel kodon testi.

Senaryo: Kalıtsal kanser yatkınlığı genleri ile kardiyomiyopati genlerini kodon
kullanımı açısından karşılaştıracaksınız. Grupları ClinVar hastalık kimlikleriyle
kurup tek komutla pakete yazacak, ortak genlerin nasıl ele alındığını ve hangi
genlerin atlandığını rapordan okuyacaksınız.

Bu örnekte:
  1. Gruplar MONDO kimlikleriyle, bütün gruplara ortak P/LP ve 2 yıldız filtresiyle tanımlanır.
  2. İki grupta birden olan genler çıkarılır (overlap="exclude").
  3. Paket yeniden okunur; her genin features'ındaki ClinVar özeti denetlenir.
  4. Dosya DifferentialCodonUsage.from_crn ile açılır ve test yapılır.
  5. Aynı paket komut satırından da açılır (crn dco ... --file).

Çıktılar: test/clinvar/cikti/clinvar_grup_4/
Not: Bu iki grup biyolojik olarak eşleştirilmiş değildir; örnek iş akışını gösterir.
Kardiyomiyopati grubunda kodonların yarısı TTN'den gelir; yorum bunu uyarı olarak
yazar. Böyle bir gen sonucu tek başına sürükleyebilir.
Sağlam kodon çıkmazsa crn dco dosya yazmaz ve "sonuç yok" bildirir.
"""
import subprocess
import sys
from collections import Counter
from pathlib import Path

from crn import read_crn_list
from crn.cub import DifferentialCodonUsage
from crn.fetch import Clinvar

CIKTI = Path(__file__).parent / "cikti" / "clinvar_grup_4"
CIKTI.mkdir(parents=True, exist_ok=True)

GRUPLAR = {
    "kanser": "MONDO:0015356",        # hereditary cancer-predisposing syndrome
    "kardiyo": "MONDO:0004994",       # cardiomyopathy
}
FILTRE = {"significance": ["pathogenic", "likely_pathogenic"], "min_stars": 2}
cv = Clinvar(verbose=False)

for ad, kimlik in GRUPLAR.items():
    print(f"{ad}: {kimlik}, filtreden geçen {len(cv.genes(kimlik, **FILTRE))} gen")

paket = CIKTI / "kanser_kardiyo.crn"
rapor = cv.pack(GRUPLAR, paket, **FILTRE)
print(f"\nClinVar {rapor['release']}, paket {paket.stat().st_size / 1e6:.0f} MB")
print("Grup büyüklükleri:", rapor["groups"])
print("Ortak olduğu için çıkarılan:", ", ".join(rapor["overlap"]) or "yok")
for gen, neden in sorted(rapor["skipped"].items()):
    print(f"  atlanan {gen}: {neden}")

genler = read_crn_list(paket)
print("\nDosyadaki gruplar:", dict(Counter(g.features.first("group") for g in genler)))
ornek = max(genler, key=lambda g: g.features.first("clinvar_pathogenic"))
print(f"En çok P/LP varyantlı gen: {ornek.features.first('gene_name')} "
      f"({ornek.features.first('clinvar_pathogenic')} P/LP, ≥2 yıldız)")
print("  hastalıklar:", "; ".join(ornek.features.get("clinvar_condition")[:3]))

dcu = DifferentialCodonUsage.from_crn(paket, "kanser", "kardiyo")
print(f"\n{dcu}")
print(dcu.interpret(permutations=2000))

print("\nAynı paket komut satırında:")
komut = [sys.executable, "-m", "crn.cli", "dco", "kanser", "kardiyo", "--file", str(paket),
         "--permutations", "2000", "-o", str(CIKTI / "robust.tsv")]
print("  " + " ".join(["crn"] + komut[3:]))
sonuc = subprocess.run(komut, capture_output=True, text=True)
print(f"  çıkış kodu {sonuc.returncode}; {(CIKTI / 'robust.tsv').exists() and 'robust.tsv yazıldı' or sonuc.stderr.strip()}")
