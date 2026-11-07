"""Grup 7: Başlangıç yaşına göre gruplar.

Senaryo: Ataksi çocuklukta da yetişkinlikte de başlayabilir. Erken başlangıçlı ve
geç başlangıçlı ataksi genlerini karşılaştıracaksınız. Yaş, hastalığın özelliğidir;
bu yüzden filtre hastalık düzeyinde uygulanır.

Bu örnekte:
  1. Başlangıç terimlerinin alt terimleri gösterilir.
  2. Erken ve geç gruplar kurulur.
  3. İki gruba birden giren genler ayrılır, kesin gruplar yazılır.
  4. Başlangıç bilgisi olmayan genlerin oranı raporlanır.

Çıktılar: test/hpo/cikti/hpo_grup_7/
Not: Çok hastalıkta başlangıç bilgisi yoktur. Bu genler iki gruba da girmez, yani
gruplar tüm ataksi genlerini temsil etmez.
"""
from pathlib import Path

from crn.fetch import Hpo

CIKTI = Path(__file__).parent / "cikti" / "hpo_grup_7"
CIKTI.mkdir(parents=True, exist_ok=True)

hpo = Hpo(verbose=False)
TERIM = "Ataxia"
KALITE = {"association": "MENDELIAN"}
ERKEN = ["Congenital onset", "Infantile onset", "Childhood onset"]
GEC = ["Adult onset"]

print("Geç başlangıç terimi şunları kapsar:")
for kimlik in hpo.descendants("Adult onset"):
    print(f"  {kimlik} {hpo.term(kimlik)['name']}")

tum = set(hpo.symbols(TERIM, **KALITE))
erken = set(hpo.symbols(TERIM, onset=ERKEN, **KALITE))
gec = set(hpo.symbols(TERIM, onset=GEC, **KALITE))
iki_yas = erken & gec
bilinmeyen = tum - erken - gec

print(f"\nTüm ataksi genleri: {len(tum)}")
print(f"  Erken başlangıç: {len(erken)}")
print(f"  Geç başlangıç:   {len(gec)}")
print(f"  İkisinde birden: {len(iki_yas)}")
print(f"  Başlangıç bilgisi yok: {len(bilinmeyen)} (%{100 * len(bilinmeyen) / len(tum):.0f})")

(CIKTI / "erken.txt").write_text("\n".join(sorted(erken - gec)) + "\n")
(CIKTI / "gec.txt").write_text("\n".join(sorted(gec - erken)) + "\n")
(CIKTI / "iki_yas.txt").write_text("\n".join(sorted(iki_yas)) + "\n")
print(f"\nKesin gruplar: erken {len(erken - gec)}, geç {len(gec - erken)}")
print("Yalnızca geç başlangıçlı:", ", ".join(sorted(gec - erken)[:10]))
