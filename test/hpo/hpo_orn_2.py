"""Örnek 2: Semptomdan güvenilir gen listesi çıkarmak.

Senaryo: Ataksi ile ilişkili genleri arıyorsunuz, ama ham liste binin üzerinde.
Filtreler eklendikçe listenin nasıl daraldığını görmek, eşiği bilinçli seçmenizi
sağlar. Her adım tek bir kriter ekler.

Bu örnekte:
  1. Terim alt terimleriyle sorgulanır, sonra yalnızca kendisiyle karşılaştırılır.
  2. Sıklık, ilişki tipi, kanıt ve kalıtım filtreleri sırayla eklenir.
  3. Son liste TSV olarak kaydedilir.

Çıktılar: test/hpo/cikti/hpo_orn_2/
Not: HPO'da gen, fenotipe hastalık üzerinden bağlıdır. Liste "bu semptomun
görüldüğü hastalıklarla ilişkili genler" olarak okunmalıdır.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_orn_2"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
TERIM = "Ataxia"

adimlar = [
    ("Yalnızca terimin kendisi", {"descendants": False}),
    ("Alt terimler dahil", {}),
    ("+ sıklık en az %30", {"min_frequency": 0.3}),
    ("+ Mendel tipi ilişki", {"min_frequency": 0.3, "association": "MENDELIAN"}),
    ("+ yayın ya da uzman kanıtı", {"min_frequency": 0.3, "association": "MENDELIAN",
                                    "evidence": ["PCS", "TAS"]}),
    ("+ otozomal dominant", {"min_frequency": 0.3, "association": "MENDELIAN",
                             "evidence": ["PCS", "TAS"],
                             "inheritance": "Autosomal dominant inheritance"}),
]

print(f"{TERIM} için gen sayısı:")
for ad, filtreler in adimlar:
    print(f"  {len(hpo.symbols(TERIM, **filtreler)):>5}  {ad}")

_, son_filtre = adimlar[-1]
genler = hpo.genes(TERIM, **son_filtre)
Hpo.write_tsv(genler, CIKTI / "ataksi_ad.tsv")

print("\nEn çok hastalıkla ilişkili genler:")
for satir in genler[:8]:
    print(f"  {satir['gene_symbol']:<10} {satir['diseases']} hastalık, {satir['terms']} terim, "
          f"en yüksek sıklık {satir['max_frequency']:.2f}")
print(f"\nYazıldı: {CIKTI / 'ataksi_ad.tsv'}")
