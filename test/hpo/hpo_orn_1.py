"""Örnek 1: Bir semptomu ontolojide bulmak.

Senaryo: Klinikte "nöbet geçiriyor" denen bir hastalık grubuyla çalışacaksınız.
HPO terimleri İngilizcedir ve hiyerarşiktir. Önce doğru terimi bulmanız, sonra
sorgunun hangi alt terimleri kapsayacağını görmeniz gerekir.

Bu örnekte:
  1. HPO verisi ilk sorguda indirilir (tqdm ilerleme çubuğuyla).
  2. Metinle terim aranır; halk dilindeki eş anlamlılar da bulunur.
  3. Tek terimin tanımı, üst ve alt terimleri okunur.
  4. Eski kimliklerin güncel terime çevrildiği gösterilir.

Çıktı: yalnızca ekran.
Not: İlk çalıştırma yaklaşık 66 MB indirir, sonrakiler internet istemez.
"""
from crn.fetch import Hpo

hpo = Hpo()

ozet = hpo.summary()
print(f"Sürüm {ozet['release']}: {ozet['terms']} terim, {ozet['genes']} gen, {ozet['diseases']} hastalık\n")

print("'seizure' araması, en çok gen içeren 6 terim:")
for satir in hpo.search("seizure", limit=6):
    print(f"  {satir['id']}  {satir['genes']:>5} gen  {satir['name']}")

print("\nHalk dilindeki ifade de bulunur ('fainting'):")
for satir in hpo.search("fainting", limit=3):
    print(f"  {satir['id']}  {satir['name']}  (eşleşen: {satir['matched']})")

nobet = hpo.term("Seizure")
print(f"\n{nobet['id']} {nobet['name']}")
print(f"  Tanım: {nobet['definition'][:100]}...")
print(f"  Doğrudan alt terim: {len(nobet['children'])}, tüm alt terimler: {len(hpo.descendants('Seizure'))}")
print("  Üst terimler:")
for kimlik in hpo.ancestors("Seizure"):
    print(f"    {kimlik} {hpo.term(kimlik)['name']}")

eski = "HP:0001275"
print(f"\nEski kimlik {eski} -> {hpo.term(eski)['id']} {hpo.term(eski)['name']}")
