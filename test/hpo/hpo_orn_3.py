"""Örnek 3: Hastalık ve gen tarafından sorgulamak.

Senaryo: Elinizde bir hastalık ya da aday gen var. Hastalığın genlerini, genin
hangi hastalık ve fenotiplerle ilişkili olduğunu öğrenmek istiyorsunuz.

Bu örnekte:
  1. Adında "parkinson" geçen hastalıklar kalıtım ve başlangıç bilgisiyle listelenir.
  2. Bir hastalığın genleri çıkarılır.
  3. Adayı SCN1A olan bir olguda genin fenotipleri ve hastalıkları incelenir.
  4. Hastalıkta olmadığı yazılı kayıtların sonuçtan çıkarıldığı gösterilir.

Çıktı: yalnızca ekran.
"""
from crn.fetch import Hpo

hpo = Hpo(verbose=False)

print("Parkinson hastalıkları (Mendel tipi, en çok gen önce):")
for satir in hpo.diseases("parkinson", association="MENDELIAN")[:6]:
    print(f"  {satir['disease_id']:<13} {satir['genes']:>3} gen  {satir['disease_name'][:42]:<44}"
          f"{satir['onset'] or '-'}")

print("\nOMIM:168600 genleri:", ", ".join(hpo.symbols(disease="OMIM:168600")))

print("\nSCN1A, en çok hastalıkta geçen fenotipler:")
for satir in hpo.phenotypes("SCN1A")[:6]:
    print(f"  {satir['hpo_id']}  {satir['diseases']} hastalık  {satir['hpo_name']}")

print("\nSCN1A hastalıkları:")
for satir in hpo.diseases(gene="SCN1A")[:6]:
    print(f"  {satir['disease_id']:<13} {satir['disease_name']}")

print("\nFEM1B geninde nöbet: tek kayıt 0/5, yani hastaların hiçbirinde görülmemiş.")
print("  Sorguda FEM1B var mı:", "FEM1B" in hpo.symbols("Seizure"))
ozet = hpo.summary()
print(f"  Yüklemede çıkarılan: {ozet['removed_not']} NOT satırı, "
      f"{ozet['removed_zero_frequency']} sıfır sıklık satırı")
