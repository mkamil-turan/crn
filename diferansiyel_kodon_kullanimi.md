# Diferansiyel kodon kullanımı (crn.cub, J bölümü)

## Ne işe yarar

Diferansiyel kodon kullanımı testi, iki gen grubunu alıp her eş anlamlı kodon
için şu soruyu yanıtlar: bu kodon, bir grupta ötekine göre anlamlı biçimde daha
çok mu yoksa daha az mı seçiliyor? Buradaki ölçü çıplak sayı değildir. Uzun ya
da kalabalık bir grup her kodonu daha çok içerir; o yüzden kodonun mutlak geçiş
sayısı değil, eş anlamlı rakipleri arasındaki payı karşılaştırılır. Aynı amino
asidi kodlayan kodonlar arasında grupların tercihinin nerede ayrıştığını
gösterir.

Tipik kullanım, bir ilgi grubunu bir arka plana karşı koymaktır: hastalık
genlerine karşı kontroller, yüksek ifadeli genlere karşı düşük ifadeli genler,
bir dokuda zenginleşen genlere karşı genom geneli, ya da bir virüsün genlerine
karşı konak genleri. Çıktı, hangi kodonların bu iki durum arasında ayrım
yarattığını kodon kodon söyler.

crn'de bunun bir ucu zaten var: iki `CUB` nesnesi birleştirilebiliyor ve
optimal kodon seçiminde ki-kare kullanılabiliyor (`method="chi2"`). Ama o,
yalnızca hangi kodonların tercih edildiğine karar verip bırakan bir elemedir;
her kodon için yön, etki büyüklüğü ve çoklu teste göre düzeltilmiş anlamlılık
vermez. Diferansiyel test bu boşluğu doldurur.


## Nasıl kullanılır

Kullanıcı iki gen seti hazırlar; her biri bir `CUB` nesnesidir. Örneğin yüksek
ifadeli ribozomal genler bir set, koku reseptörü genleri öteki set olur.
Test bu iki seti ve bir karşılaştırma düzeyini alır, bir tablo döndürür. Her
satır bir kodondur ve kodonun amino asidini, iki gruptaki sayımlarını, olasılık
oranını, ham p değerini ve düzeltilmiş q değerini taşır.

Karşılaştırma üç düzeyde yapılabilir ve düzey seçimi sonucun anlamını belirler.

Küresel düzeyde kodon, diziideki tüm öteki kodonlara karşı konur. Bu, kodonun
genel bolluğunu ölçer ama iki grubun amino asit bileşimi farklıysa karışır;
çünkü bir amino asit bir grupta daha sıksa, onun tüm kodonları da sırf bu
yüzden yüksek çıkar.

Aile düzeyi esas ilgi çeken düzeydir. Kodon, yalnızca aynı amino asidi kodlayan
öteki eş anlamlı kodonlara karşı konur; örneğin GCC yalnızca GCU, GCA ve GCG
ile yarışır. Böylece amino asit bileşiminden arınılır ve geriye tam da sinonim
tercih kalır: aynı Ala için hangi kodon seçiliyor.

Alt aile düzeyi aileyi daha da böler. Eş anlamlı kodonlar ilk iki bazı ortak
olanlara ayrılır; Löysin'de CUN grubuyla UUR grubu ayrı değerlendirilir. Bu,
wobble konumundaki ince tercihleri yakalamak içindir ve çoğu zaman aile
düzeyinin bir ayrıntılandırmasıdır.

Pratik akış şudur: iki set kurulur, uygun düzeyle test çağrılır, dönen
tablodan q değeri eşiğin altındaki satırlar alınır ve olasılık oranının
büyüklüğüne göre sıralanır. Sonuç bir yanardağ grafiğiyle de özetlenebilir;
yatay eksende etki yönü, dikey eksende anlamlılık.


## Ne anlama gelir

Her kodon için bir 2x2 sayım tablosu kurulur. Grup A ve grup B satırlar,
ilgilenilen kodon ile karşılaştırma havuzundaki geri kalan sütunlardır.

                 kodon      geri kalan
    Grup A         a            b
    Grup B         c            d

Geri kalanın ne olduğu, seçilen düzeye göre değişir: küresel düzeyde tüm öteki
kodonlar, aile düzeyinde aynı amino asidin öteki kodonları, alt aile düzeyinde
aynı alt ailenin öteki kodonları.

Bu tablodan üç şey okunur.

Birincisi olasılık oranı, etkinin yönü ve büyüklüğü. Grup A'da kodonun geri
kalana oranının, grup B'deki aynı orana bölümüdür:

    OR = (a/b) / (c/d) = ad / bc

LibreOffice Math söz dizimiyle:

    OR = {a/b} over {c/d} = {a cdot d} over {b cdot c}

Bir hücre sıfır olduğunda oran tanımsız kalmasın diye her hücreye 0,5 eklenir
(Haldane-Anscombe düzeltmesi) ve genellikle ikili logaritması raporlanır:

    OR = ((a+0.5)(d+0.5)) / ((b+0.5)(c+0.5))

LibreOffice Math söz dizimiyle:

    OR = {(a+0.5) cdot (d+0.5)} over {(b+0.5) cdot (c+0.5)}

İkinci okunan değer p değeridir; Fisher kesin testiyle hesaplanır. Satır ve
sütun toplamları sabit tutulduğunda, gözlenen ve ondan daha uç tabloların
hipergeometrik olasılığı toplanır. Tek bir tablonun olasılığı:

    P = C(a+b, a) · C(c+d, c) / C(N, a+c)

burada N = a+b+c+d. LibreOffice Math söz dizimiyle:

    P = {C(a+b,a) cdot C(c+d,c)} over {C(N, a+c)}

Fisher kesin testi, ki-karenin aksine küçük sayılarda da geçerlidir; nadir
kodonlarda ve küçük gen gruplarında ki-karenin yaklaşıklığı bozulduğu için
tercih edilir.

Üçüncü değer q değeridir, yani çoklu karşılaştırmaya göre düzeltilmiş
anlamlılık. Bir düzeyde elliyi aşkın kodon ayrı ayrı test edildiğinden,
yalnızca şans eseri bazıları düşük p verir. Benjamini-Hochberg düzeltmesi
yanlış keşif oranını denetler; p değerleri küçükten büyüğe sıralanıp her sıraya
oranlı bir eşik uygulanır:

    q_(i) = min_{k ≥ i} [ (m / k) · p_(k) ]

m o düzeydeki test sayısıdır. LibreOffice Math söz dizimiyle:

    q_i = min from {k >= i} left ( {m} over {k} cdot p_k right )


## Sonuç nasıl yorumlanır

Yorum, yönü olasılık oranından, güveni q değerinden alır. İkisini birlikte
okumak gerekir; tek başına düşük q anlamlılığı söyler ama etkinin ne kadar
büyük olduğunu söylemez, tek başına büyük olasılık oranı da küçük sayılardan
gelen gürültü olabilir.

q değeri eşiğin (çoğunlukla 0,05) altındaysa kodon o iki grup arasında anlamlı
farklı kullanılıyordur. Olasılık oranı birden büyükse (log2 OR pozitif) kodon
grup A'da göreli daha çok seçiliyor, birden küçükse grup B'de daha çok
seçiliyor demektir. Örneğin yüksek ifadeli genleri düşük ifadeli genlere karşı
koyduğunuzda, bir amino asidin tercih edilen kodonunun grup A lehine yüksek
oranla ve düşük q ile çıkması, o kodonun çeviri için optimal kodon olduğunun
işaretidir; aynı amino asidin öteki kodonu ise ters yönde, oran bir altında
çıkar. Böyle bir çift, sinonim tercihin yönünü doğrudan gösterir.

Düzeyler arası okuma da anlamlıdır. Aile düzeyinde anlamlı çıkan ama küresel
düzeyde kaybolan bir kodon, farkın saf sinonim tercihten geldiğini; küresel
düzeyde güçlü ama aile düzeyinde zayıflayan bir kodon ise farkın büyük ölçüde
amino asit bileşiminden kaynaklandığını düşündürür. Alt aile düzeyinde beliren
bir ayrım, tercihin wobble konumuna indiğini gösterir.

Birkaç uyarı yorumu güvenli tutar. Her iki grupta da neredeyse hiç geçmeyen bir
kodonun olasılık oranı Haldane düzeltmesiyle hesaplanabilir ama güveni düşüktür;
bu satırlar q yüksek çıkar ve zaten elenir. Olasılık oranı büyük görünüp q
yüksekse, etki küçük örneklemden gelen bir yanılsamadır. Test bir ilişkiyi
gösterir, nedeni değil: bir kodonun bir grupta zenginleşmesi çeviri seçilimi,
mutasyon yanlılığı ya da GC bileşimi farkı gibi birkaç nedenden olabilir; hangisi
olduğu, GC3s, ENC ya da tAI gibi öteki crn indeksleriyle birlikte okunarak
yorumlanır.
