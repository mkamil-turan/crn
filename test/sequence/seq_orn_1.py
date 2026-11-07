"""Örnek 1: İlk dizini tanımla, not ekle, kaydet ve geri oku.

Senaryo: Laboratuvarda çalıştığınız insan HBB (beta-globin) kodlayan dizisini
bilgisayarda tutmak istiyorsunuz. Dizinin yanında gen adı, organizma ve deney
notlarını da saklamak gerekiyor.

Bu örnekte:
  1. DNA nesnesi oluşturulur; küçük harfli giriş otomatik düzeltilir.
  2. Hatalı karakter içeren dizinin nasıl reddedildiği gösterilir.
  3. Bilgiler features içine yazılır.
  4. Dizi hem FASTA hem CRN olarak kaydedilir ve geri okunur.
  5. İki biçimin farkı karşılaştırılır: FASTA notları kaybeder, CRN korur.
  6. Bozulan bir CRN dosyasının okunurken yakalandığı gösterilir.

"""
from pathlib import Path

from crn import DNA, Sequence

CIKTI = Path(__file__).parent / "cikti" / "seq_orn_1"
CIKTI.mkdir(parents=True, exist_ok=True)

# 1. Dizi tanımlama: küçük harf yazılsa da büyük harfe çevrilir.
hbb = DNA(
    "atggtgcatctgactcctgaggagaagtctgccgttactgccctgtggggcaaggtgaacgtggatgaagtt"
    "ggtggtgaggccctgggcaggctgctggtggtctacccttggacccagaggttctttgagtcctttggggat"
    "ctgtccactcctgatgctgttatgggcaaccctaaggtgaaggctcatggcaagaaagtgctcggtgccttt"
    "agtgatggcctggctcacctggacaacctcaagggcacctttgccacactgagtgagctgcactgtgacaag"
    "ctgcacgtggatcctgagaacttcaggctcctgggcaacgtgctggtctgtgtgctggcccatcactttggc"
    "aaagaattcaccccaccagtgcaggctgcctatcagaaagtggtggctggtgtggctaatgccctggcccac"
    "aagtatcactaa",
    name="HBB_CDS",
)
print("Dizi       :", hbb)
print("Uzunluk    :", len(hbb), "nt")
print("İlk kodon  :", hbb[:3].seq)
print("EcoRI var mı (GAATTC):", "GAATTC" in hbb)

# 2. Alfabe denetimi: yazım hataları analize sızmadan yakalanır.
try:
    DNA("ATGXXGCA", name="hatali")
except ValueError as hata:
    print("Reddedildi :", hata)

# 3. Notlar. add aynı anahtara değer ekler, set tek değerle değiştirir.
hbb.features.set("gene", "HBB")
hbb.features.set("organism", "Homo sapiens")
hbb.features.set("cds", (1, len(hbb)))
hbb.features.add("note", "Orak hücre çalışması için referans")
hbb.features.add("note", "Primer tasarımında kullanıldı")

# 4. Kaydetme.
fasta_yolu = CIKTI / "hbb.fasta"
crn_yolu = CIKTI / "hbb.crn"
hbb.to_fasta(fasta_yolu)
hbb.to_crn(crn_yolu)

fastadan = DNA.from_fasta_file(fasta_yolu)[0]  # FASTA birden çok kayıt tutabilir, liste döner
crnden = Sequence.from_crn(crn_yolu)           # sınıf dosyadan anlaşılır

# 5. Karşılaştırma.
print("\nDizi aynı mı     FASTA:", fastadan == hbb, "| CRN:", crnden == hbb)
print("Sınıf            FASTA:", type(fastadan).__name__, "| CRN:", type(crnden).__name__)
print("Gen notu         FASTA:", fastadan.features.first("gene"),
      "| CRN:", crnden.features.first("gene"))
print("Tüm notlar (CRN) :", crnden.features.get("note"))

# 6. Bozuk dosya. Tek bayt değişse bile sağlama toplamı (CRC32) bunu yakalar.
bozuk_yolu = CIKTI / "hbb_bozuk.crn"
veri = bytearray(crn_yolu.read_bytes())
veri[len(veri) // 2] ^= 0xFF
bozuk_yolu.write_bytes(veri)
try:
    Sequence.from_crn(bozuk_yolu)
except ValueError as hata:
    print("\nBozuk dosya      :", hata)

print("\nKaydedilen dosyalar:", CIKTI)
