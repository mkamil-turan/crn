"""Örnek 3: Periferik sinir ve retina; CMT ile retinitis pigmentosa.

Klinik arka plan: Charcot-Marie-Tooth (CMT) en sık kalıtsal periferik nöropatidir; uzun
aksonlar ve Schwann hücreleri bozulur. Retinitis pigmentosa (RP) fotoreseptörlerin
ilerleyici kaybıyla körlüğe gider. İki hastalık da tek tip hücrenin yavaş kaybıdır ama
dokuları çok farklıdır.

Soru: Dokuya özgü kalıtsal dejenerasyonlar farklı kodon profili taşır mı? Fotoreseptör
yoğun protein üretir (dış segment her gün yenilenir), bu yüzden RP genlerinde çeviriye
uyumun farklı olması beklenebilir.

Bu örnekte:
  1. CMT ve RP ayrı ayrı kendi eşleştirilmiş kontrolüyle karşılaştırılır.
  2. İki hastalık doğrudan karşılaştırılır. Bu eşleştirme içermediği için önce
     uzunluk ve GC farkına bakılır, kodon farkı buna göre okunur.

Çıktılar: test/dco/cikti/dco_3_cmt, dco_3_rp, dco_3_dogrudan
Not: Kontrol havuzunda göz notu bulunan genler yoktur, yani RP için de geçerlidir.
"""
import re

from dco_ortak import ILISKI, analiz, gen_adi, gen_al, hazirla

k = hazirla()
filtre = {"association": ILISKI}
cmt = gen_al(k, disease=re.compile("charcot-marie-tooth", re.I), **filtre)
rp = gen_al(k, disease=re.compile("retinitis pigmentosa", re.I), **filtre)

# Bir gen iki hastalıkta birden geçiyorsa doğrudan karşılaştırmadan çıkarılır.
ortak = {gen_adi(g) for g in cmt} & {gen_adi(g) for g in rp}
print(f"CMT {len(cmt)} gen, RP {len(rp)} gen, ortak: {sorted(ortak) or 'yok'}\n")

analiz(k, "CMT", cmt, "dco_3_cmt")
analiz(k, "RP", rp, "dco_3_rp")
analiz(k, "CMT", [g for g in cmt if gen_adi(g) not in ortak], "dco_3_dogrudan",
       karsi=[g for g in rp if gen_adi(g) not in ortak], kontrol_adi="RP")
