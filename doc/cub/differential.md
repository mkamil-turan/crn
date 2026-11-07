# `crn.cub.differential`: J. Diferansiyel kodon kullanımı

Kaynak: [`crn/cub/differential.py`](../../crn/cub/differential.py)

İki gen grubunu kodon kodon karşılaştırır: bir kodon, karşılaştırma havuzu içinde bir grupta ötekine göre anlamlı biçimde daha çok mu seçiliyor? Her kodon için bir 2×2 tablo kurulur.

```text
             kodon      geri kalan
Grup A         a            b
Grup B         c            d
```

Yöntemin gerekçesi ve ayrıntılı yorumu [`diferansiyel_kodon_kullanimi.md`](../../diferansiyel_kodon_kullanimi.md) dosyasındadır. Uçtan uca örnekler `test/dco/` altında durur.

<!-- toc -->
- [1. Karşılaştırma düzeyi](#1-karşılaştırma-düzeyi)
- [2. Raporlanan değerler](#2-raporlanan-değerler)
- [3. Neden permütasyon gerekir](#3-neden-permütasyon-gerekir)
- [4. Kurulum](#4-kurulum): [`DifferentialCodonUsage`](#differentialcodonusage) · [`from_crn`](#differentialcodonusagefrom_crn)
- [5. Sonuç tabloları](#5-sonuç-tabloları): [`table`](#differentialcodonusagetable) · [`significant`](#differentialcodonusagesignificant) · [`permutation`](#differentialcodonusagepermutation) · [`robust`](#differentialcodonusagerobust) · [`compare_levels`](#differentialcodonusagecompare_levels) · [`interpret`](#differentialcodonusageinterpret) · [`plot_volcano`](#differentialcodonusageplot_volcano)
- [6. Modül düzeyinde fonksiyonlar](#6-modül-düzeyinde-fonksiyonlar): [`fisher_log10_p`](#fisher_log10_p) · [`benjamini_hochberg_log10`](#benjamini_hochberg_log10)
- [7. İç yardımcılar](#7-iç-yardımcılar)
<!-- /toc -->

## 1. Karşılaştırma düzeyi

"Geri kalan" düzeye bağlıdır ve düzey seçimi sonucun anlamını belirler. `LEVELS` kabul edilen değerlerdir.

| Düzey | Geri kalan | Ne ölçer |
|---|---|---|
| `global` | Bütün öteki anlamlı kodonlar | Genel bolluk; amino asit bileşimi farkıyla karışır |
| `family` | Aynı amino asidin öteki kodonları | Saf sinonim tercih; esas düzey budur |
| `subfamily` | Aynı amino asit ve aynı ilk iki baz | Wobble konumundaki tercih |

Stoplar hiçbir düzeye girmez, tek kodonlu amino asitler aile ve alt aile düzeylerinde test edilmez.

## 2. Raporlanan değerler

Olasılık oranı Haldane-Anscombe düzeltmesiyle hesaplanır, `OR = ((a+0.5)(d+0.5)) / ((b+0.5)(c+0.5))`. `log2_or` pozitifse kodon A grubunda, negatifse B grubunda göreli daha çok seçilir. p iki yönlü Fisher kesin testidir (`scipy.stats.fisher_exact`). q o düzeydeki testler üzerinde Benjamini-Hochberg düzeltmesidir.

Genom ölçeğinde p değerleri çift duyarlıklı sayının altına inebilir. `fisher_exact` yaklaşık 1e-295'in altında kuyruklardan birini kaybeder, daha aşağıda 0 döndürür. Bu yüzden p `_SCIPY_MIN_P` (1e-200) altına indiğinde aynı toplam log uzayında yeniden hesaplanır: `p` ve `q` 0'a yuvarlansa bile `log10_p` ve `log10_q` sonludur.

## 3. Neden permütasyon gerekir

Fisher testi sayımları grup içinde havuzlar ve her kodonu bağımsız gözlem sayar. Oysa genler kodon kullanımında birbirinden zaten farklıdır. Rastgele bölünmüş gen gruplarında 59 kodondan 53'üne kadarı anlamlı çıkabilir, yani tek başına Fisher sonucu bir bulgu değildir. `permutation` genlerin grup etiketlerini karıştırır ve log2 OR'un gerçek değerini rastgele gruplarla karşılaştırır. Aynı rastgele bölünmelerde permütasyon hiç kodon bulmaz, gerçek fark eklenen pozitif kontrolde 59 kodonun hepsini bulur. Fisher'da da permütasyonda da geçen kodonlar `robust` işaretlenir, `interpret` yalnızca bunları yorumlar.

Permütasyon genler arası çeşitliliği hesaba katar, ama iki grubun uzunluk ya da bölgesel GC farkını gidermez. Kontrol grubu bu özelliklere göre eşleştirilmelidir; `test/dco/dco_ortak.py` bunun bir örneğidir.

## 4. Kurulum

### DifferentialCodonUsage

`DifferentialCodonUsage(group_a, group_b, name_a='A', name_b='B', genetic_code=GeneticCode.STANDARD)`

`group_a` ve `group_b` `CUB`, `Gene` listesi ya da `write_crn_list` ile yazılmış dosya yoludur. `name_a` ve `name_b` tablo ve yorumda kullanılacak adlardır. `genetic_code` grup `CUB` olarak verilmediğinde kullanılır. Grup bu tiplerden biri değilse `TypeError`; genetik kodlar farklıysa ya da bir gen iki grupta birden varsa `ValueError` verir. Genler `CUB` kuralına uymalıdır.

```python
from crn import write_crn_list
from crn.fetch import Gencode
from crn.cub import CUB, DifferentialCodonUsage

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

ribozomal = kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2", "RPL5", "RPS8"])
koku = kur(["OR1A1", "OR2A2", "OR5A1", "OR10A2", "OR4D1", "OR6A2", "OR1E1", "OR2B2"])
dcu = DifferentialCodonUsage(CUB(ribozomal), koku, name_a="ribozom", name_b="koku")
print(dcu)  # -> DifferentialCodonUsage: ribozom (8 gen) / koku (8 gen)
```

### DifferentialCodonUsage.from_crn

`sınıf metodu: DifferentialCodonUsage.from_crn(file_name, group_a, group_b, key='group', genetic_code=GeneticCode.STANDARD)`

Tek CRN listesindeki genleri features'taki grup bilgisine göre ayırır. `file_name` gen listesidir, `group_a` ve `group_b` grup adlarıdır, `key` grup adının tutulduğu anahtardır, `genetic_code` kullanılacak koddur. Dosyada `Gene` olmayan kayıt varsa `TypeError`, gruplardan biri boşsa dosyadaki grupları sayan `ValueError` verir.

```python
for gen in ribozomal:
    gen.features.set("grup", "ribozom")
for gen in koku:
    gen.features.set("grup", "koku")
write_crn_list(ribozomal + koku, "gruplar.crn")

dcu = DifferentialCodonUsage.from_crn("gruplar.crn", "ribozom", "koku", key="grup")
print(dcu)  # -> DifferentialCodonUsage: ribozom (8 gen) / koku (8 gen)
```

## 5. Sonuç tabloları

Bu sayfadaki permütasyon örnekleri kısa sürsün diye `n=200` kullanır. Bu değer aile düzeyinin 59 testi için yetersizdir, `perm_q` eşiğe ulaşamaz; gerçek analizde varsayılan 5000 bırakılmalıdır.

### DifferentialCodonUsage.table

`DifferentialCodonUsage.table(level='family', q=0.05, file_name=None)`

Düzeydeki her kodon için bir satır verir. Anahtarlar: `level, codon, amino_acid, pool, a, b, c, d, share_a, share_b, odds_ratio, log2_or, p, log10_p, q, log10_q, significant, favored`. `share_*` kodonun havuzdaki payıdır, `favored` anlamlı satırda kodonu daha çok seçen grubun adıdır. `level` düzey, `q` eşik, `file_name` verilirse TSV çıktısıdır. Geçersiz değer `ValueError` verir.

```python
satirlar = dcu.table(level="family", q=0.05, file_name="tablo.tsv")
print(len(satirlar), satirlar[0]["codon"], satirlar[0]["amino_acid"])  # -> 59 UUU F
```

### DifferentialCodonUsage.significant

`DifferentialCodonUsage.significant(level='family', q=0.05, min_log2_or=0.0)`

Eşiği geçen satırları |log2 OR|'a göre büyükten küçüğe verir. `level` ve `q` önceki metotla aynıdır, `min_log2_or` satırın kalması için gereken en küçük |log2 OR|'dur.

```python
print(len(dcu.significant("family")), len(dcu.significant(level="global", q=0.01, min_log2_or=1.0)))  # -> 4 17
```

### DifferentialCodonUsage.permutation

`DifferentialCodonUsage.permutation(level='family', q=0.05, n=PERMUTATIONS, seed=0, file_name=None)`

Gen düzeyinde permütasyon testi, `table` satırlarına `perm_p`, `perm_q` ve `robust` ekler. Genler grup büyüklükleri korunarak yeniden dağıtılır ve her kodon için log2 OR yeniden hesaplanır. `perm_p`, rastgele grupların en az gerçek fark kadar uç sonuç verme oranıdır, `(uç sonuç + 1) / (n + 1)`. `robust` True ise kodon iki testte de `q` eşiğini geçmiştir.

`level` düzey, `n` karıştırma sayısı (varsayılan `PERMUTATIONS`, 5000), `seed` tohumdur; sonuç tohum için tekrar üretilebilir ve (level, n, seed) başına saklanır. `file_name` verilirse TSV yazılır. En küçük p `1 / (n + 1)` olduğundan `n`, düzeydeki test sayısının `1 / q` katından az olmamalıdır: aile düzeyinde 59 test ve q 0,05 için en az 1180 karıştırma gerekir. Geçersiz değer `ValueError` verir.

```python
satirlar = dcu.permutation("family", q=0.05, n=200, seed=1, file_name="perm.tsv")
print(sorted(set(satirlar[0]) - set(dcu.table("family")[0])))  # -> ['perm_p', 'perm_q', 'robust']
```

### DifferentialCodonUsage.robust

`DifferentialCodonUsage.robust(level='family', q=0.05, min_log2_or=0.0, n=PERMUTATIONS, seed=0)`

Fisher ve permütasyonda birlikte anlamlı satırları |log2 OR|'a göre verir. Raporlanması gereken liste budur. `level`, `q` ve `min_log2_or` `significant`, `n` ve `seed` `permutation` ile aynıdır.

```python
print(isinstance(dcu.robust("family", q=0.05, min_log2_or=0.0, n=200, seed=1), list))  # -> True
```

### DifferentialCodonUsage.compare_levels

`DifferentialCodonUsage.compare_levels(q=0.05, permutations=0, seed=0)`

Eş anlamlı kodonların düzeyler arası okumasını verir. Satır `codon, amino_acid, log2_or_global, log2_or_family, log2_or_subfamily, reading` içerir. `reading` şu değerlerden biridir: `sinonim` (ailede anlamlı, küreselde değil; fark saf sinonim tercih), `bilesim` (küreselde anlamlı, ailede değil ya da ailede etkisi yarıdan aza iniyor; fark büyük ölçüde amino asit bileşiminden), `karma` (iki düzeyde de anlamlı, etki korunuyor), `wobble` (yalnızca alt ailede anlamlı), boş (hiçbir düzeyde anlamlı değil).

`q` eşiktir. `permutations` 0 ise "anlamlı" Fisher'a göre okunur, büyükse o kadar karıştırmayla `robust` satırlar kullanılır. Fisher gen çeşitliliğini aşırı saydığı için ikincisi önerilir. `seed` permütasyon tohumudur.

```python
okuma = dcu.compare_levels(q=0.05, permutations=200, seed=1)
print(sorted(okuma[0]))  # -> ['amino_acid', 'codon', 'log2_or_family', 'log2_or_global', 'log2_or_subfamily', 'reading']
```

### DifferentialCodonUsage.interpret

`DifferentialCodonUsage.interpret(q=0.05, top=5, small_effect=0.25, min_count=10, permutations=PERMUTATIONS, seed=0)`

Sonuçların Türkçe yorumunu metin olarak verir. Yalnızca permütasyondan geçen kodonlara dayanır: Fisher'ın bulup permütasyonun doğrulamadığı kodonlar sayı olarak anılır, biyolojik yorum yapılmaz. Ardından düzeyler arası okuma, GC bağlamı ve güvenlik uyarıları gelir.

`q` eşiktir. `top` her grup için anılacak kodon sayısıdır. `small_effect` altındaki |log2 OR| "etkisi küçük" sayılır (0,25, OR yaklaşık 1,19). `min_count` iki gruptaki toplam geçişi bundan az olan kodonları düşük sayımlı sayar. `permutations` karıştırma sayısı, `seed` tohumdur. Karıştırma sayısı ya da grup büyüklüğü anlamlılığa ulaşmaya yetmiyorsa metin bunu söyler.

```python
yorum = dcu.interpret(q=0.05, top=3, small_effect=0.25, min_count=10, permutations=200, seed=1)
print(yorum.splitlines()[0][:60])  # -> ribozom (8 gen, 3,147 kodon) ile koku (8 gen, 2,553 kodon) k
```

### DifferentialCodonUsage.plot_volcano

`DifferentialCodonUsage.plot_volcano(level='family', q=0.05, labels=True, ax=None, file_name=None)`

Yanardağ grafiği: x log2 OR, y −log10 Fisher q, anlamlılar grup rengindedir. Grafik Fisher sonucunu gösterir, gen düzeyinde doğrulama için `robust` kullanın. `level` düzey, `q` yatay kılavuz çizgisidir. `labels=True` anlamlı kodonların adlarını yazar.

```python
ax = dcu.plot_volcano(level="family", q=0.05, labels=True, file_name="yanardag.png")
```

## 6. Modül düzeyinde fonksiyonlar

### fisher_log10_p

`fisher_log10_p(a, b, c, d)`

[[`a`, `b`], [`c`, `d`]] tablosunun iki yönlü Fisher testini log10(p) olarak verir. p `_SCIPY_MIN_P` altındaysa sonuç log uzayında yeniden bulunur, böylece değer her zaman sonludur. scipy ilk çağrıda yüklenir.

```python
from crn.cub.differential import fisher_log10_p, benjamini_hochberg_log10

print(round(fisher_log10_p(a=8, b=2, c=1, d=5), 4))  # -> -1.4564
print(round(fisher_log10_p(2000000, 1000000, 1000000, 2000000), 1))  # -> -147575.1
```

### benjamini_hochberg_log10

`benjamini_hochberg_log10(log10_p)`

`log10_p` değerlerinden Benjamini-Hochberg log10(q) değerlerini aynı sırayla verir.

```python
print([round(x, 3) for x in benjamini_hochberg_log10([-3.0, -1.0, -0.5])])  # -> [-2.523, -0.824, -0.5]
```

## 7. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_LOG10` | ln 10; log dönüşümü için |
| `_LARGE_EFFECT` | `interpret`'in "büyük görünüyor ama anlamlı değil" uyarısının |log2 OR| eşiği: 1,0 |
| `_DOMINANT_GENE` | Tek genin grup kodonlarındaki payı bunu aşarsa uyarı verilir: 0,25 |
| `_MIN_THIRD_BASE` | GC3 ayrışması yorumu için gereken en az anlamlı kodon: 3 |
| `_GC3S_GAP` | GC3s farkı bunu aşarsa GC uyarısı eklenir: 0,05 |
| `_log_comb` | Binom katsayısının logaritması |
| `_log_add` | log(e^x + e^y); çok küçük değerlerde taşmaz |
| `_fisher_log10_tails` | Fisher p'sini tümüyle log uzayında hesaplar; en olası tablodan uzaklaşarak ardışık oranlarla toplar, hesap olasılık küçülünce biter |
| `_power10` | 10^x, alt taşmada 0 |
| `DifferentialCodonUsage._as_cub` | Grup girdisini `CUB`'a çevirir |
| `DifferentialCodonUsage._pools` | Düzeyin test birimlerini (havuz adı, kodonlar) verir |
| `DifferentialCodonUsage._compute` | Düzeyin tablosunu hesaplar ve saklar |
| `DifferentialCodonUsage._gene_counts` | Gen × kodon sayım matrisi (numpy) |
| `DifferentialCodonUsage._permute` | Permütasyon döngüsü; `perm_p` ve `perm_q` üretir ve saklar |
| `DifferentialCodonUsage._smallest_p` | Ulaşılabilecek en küçük `perm_p`; küçük grupta ayrık bölünme sayısıyla sınırlıdır, eşit gruplarda iki etiketin yer değiştirmesi aynı sonucu verdiği için sayı yarıya iner |
| `DifferentialCodonUsage._opposing_pairs` | Aynı amino asitte zıt yönde anlamlı en güçlü kodon çiftleri |
| `DifferentialCodonUsage._third_base` | İki grup lehine anlamlı kodonlarda 3. konumu G/C olanların sayısı |
| `DifferentialCodonUsage._dominant_gene` | Grubun en uzun geni ve kodon payı |
| `DifferentialCodonUsage.__repr__` | `DifferentialCodonUsage: ribozom (8 gen) / koku (8 gen)` |
