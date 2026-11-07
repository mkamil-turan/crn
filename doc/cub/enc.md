# `crn.cub.enc`: F. Etkin kodon sayısı ve grafikler

Kaynak: [`crn/cub/enc.py`](../../crn/cub/enc.py)

Etkin kodon sayısı (ENC, Nc; Wright 1990) yanlılığın referans gerektirmeyen ölçüsüdür: 20 her amino asit için tek kodon, 61 hiç yanlılık yok demektir. Her amino asit için F (homozigotluk) hesaplanır, aynı dejenerelik sınıfındaki (aile büyüklüğü 2, 3, 4, 6) amino asitlerin F'leri ortalanır ve Nc = (tek kodonlu amino asit sayısı) + Σ_k (sınıftaki amino asit sayısı) / F̄_k bulunur. En az iki kez kullanılmamış (n < 2) ya da F'si 0 çıkan amino asitler hesaba girmez. Bir sınıfta hiç F kalmazsa: 3'lü sınıf (standart kodda yalnızca Ile) için F3 = (F2 + F4) / 2 alınır, öbür sınıflarda ENC hesaplanamaz ve `nan` döner. Sonuç amino asit kodlayan kodon sayısını aşarsa o değere indirilir. Kısa genlerde `nan` görmek bu yüzden normaldir.

Nötrallik ve ENC-GC3s grafikleri genler arası karşılaştırmadır, her gen bir noktadır. Nötrallik grafiğinde eğim 1'e yakınsa mutasyon baskısı, 0'a yakınsa seçilim baskın kabul edilir. PR2 varsayılan olarak dört katlı dejenere kodonlarda hesaplanır, (0,5; 0,5) merkezi A=U ve G=C dengesidir.

`gene`, `ax`, `file_name`, `labels` ve `groups` parametreleri [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`EffectiveCodons`](#effectivecodons) · [`homozygosity`](#effectivecodonshomozygosity) · [`enc`](#effectivecodonsenc) · [`expected_enc`](#effectivecodonsexpected_enc) · [`enc_ratio`](#effectivecodonsenc_ratio) · [`enc_points`](#effectivecodonsenc_points) · [`neutrality`](#effectivecodonsneutrality) · [`pr2`](#effectivecodonspr2) · [`pr2_points`](#effectivecodonspr2_points)
- [2. Grafikler](#2-grafikler): [`plot_enc`](#effectivecodonsplot_enc) · [`plot_enc_ratio`](#effectivecodonsplot_enc_ratio) · [`plot_neutrality`](#effectivecodonsplot_neutrality) · [`plot_pr2`](#effectivecodonsplot_pr2)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sınıf

### EffectiveCodons

`EffectiveCodons(cub)`

`cub` bir `CUB` olmalıdır, değilse `TypeError` verilir. `EffectiveCodons.composition` içteki `Composition` nesnesidir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, EffectiveCodons

gc = Gencode(verbose=False)
cub = CUB([gc.canonical(s).isoforms()[0] for s in ["TP53", "HBB", "ACTB", "MEFV", "RPL3", "OR1A1"]])
ad = cub.names[0]
eff = EffectiveCodons(cub)
```

### EffectiveCodons.homozygosity

`EffectiveCodons.homozygosity(gene=None)`

Madde 23. Her eş anlamlı amino asidin F değeri, (nΣp² − 1) / (n − 1). n < 2 olan amino asitler dönmez.

```python
print(round(eff.homozygosity(ad)["L"], 3))  # -> 0.274
```

### EffectiveCodons.enc

`EffectiveCodons.enc(gene=None)`

Madde 23. Etkin kodon sayısı Nc. Bir dejenerelik sınıfı için F bulunamazsa `nan` döner.

```python
print(round(eff.enc(ad), 2), round(eff.enc(), 2))  # -> 51.23 45.48
```

### EffectiveCodons.expected_enc

`statik: EffectiveCodons.expected_enc(gc3s)`

Madde 24. Yanlılık yalnızca GC3s'ten gelseydi beklenecek Nc: 2 + s + 29 / (s² + (1 − s)²). `gc3s` s değeridir. Statik metottur.

```python
print(round(EffectiveCodons.expected_enc(gc3s=0.5), 2))  # -> 60.5
```

### EffectiveCodons.enc_ratio

`EffectiveCodons.enc_ratio(gene=None)`

Madde 25. (ENC_beklenen − ENC_gözlenen) / ENC_beklenen. 0'dan uzaklaştıkça yanlılığın GC3s dışı bir nedeni vardır.

```python
print(round(eff.enc_ratio(ad), 3))  # -> 0.116
```

### EffectiveCodons.enc_points

`EffectiveCodons.enc_points()`

Madde 24-25. Her gen için GC3s, gözlenen ve beklenen ENC ve ENC oranı.

```python
print(sorted(eff.enc_points()[0]))  # -> ['enc', 'enc_ratio', 'expected_enc', 'gc3s', 'gene']
```

### EffectiveCodons.neutrality

`EffectiveCodons.neutrality()`

Madde 26. Her gen için (GC3, GC12) ve GC12 = eğim·GC3 + kesişim. `points`, `slope`, `intercept`, `r2` döner; iki genden azsa doğru `nan` olur.

```python
n = eff.neutrality()
print(round(n["slope"], 3), round(n["r2"], 3), len(n["points"]))  # -> 0.039 0.004 6
```

### EffectiveCodons.pr2

`EffectiveCodons.pr2(gene=None, codons=None)`

Madde 27. PR2 konumu: x = G3 / (G3 + C3), y = A3 / (A3 + U3); payda 0 ise 0,5. `codons` üçüncü konumu sayılacak kodonlardır, verilmezse dört katlı dejenere kodonlar kullanılır.

```python
print({k: round(v, 3) for k, v in eff.pr2(ad).items()})  # -> {'x': 0.393, 'y': 0.425}
print(eff.pr2(ad, codons=["GCG"])["x"])  # -> 1.0
```

### EffectiveCodons.pr2_points

`EffectiveCodons.pr2_points(codons=None)`

Madde 27. Her gen için (x, y). `codons` `pr2` ile aynıdır.

```python
print(sorted(eff.pr2_points()[0]), len(eff.pr2_points(codons=["GCG", "GCC"])))  # -> ['gene', 'x', 'y'] 6
```

## 2. Grafikler

### EffectiveCodons.plot_enc

`EffectiveCodons.plot_enc(ax=None, labels=False, groups=None, file_name=None)`

ENC-GC3s grafiği: genler ve Wright'ın beklenen ENC eğrisi. x 0-1, y `ENC_AXIS_LIMITS` (20-61) aralığındadır, eğri `ENC_CURVE_STEPS` noktayla çizilir. Eğrinin altında kalan genlerde yanlılık yalnızca GC3s ile açıklanamaz. ENC'si `nan` olan genler çizilmez.

```python
ax = eff.plot_enc(labels=True, file_name="enc.png")
```

### EffectiveCodons.plot_enc_ratio

`EffectiveCodons.plot_enc_ratio(ax=None, bin_width=ENC_RATIO_BIN_WIDTH, groups=None, file_name=None)`

ENC oranının frekans dağılımı. `bin_width` çubuk genişliğidir, varsayılanı `ENC_RATIO_BIN_WIDTH` (0,05).

```python
ax = eff.plot_enc_ratio(bin_width=0.1, file_name="oran.png")
```

### EffectiveCodons.plot_neutrality

`EffectiveCodons.plot_neutrality(ax=None, labels=False, groups=None, file_name=None)`

Nötrallik grafiği: GC12 ~ GC3, regresyon doğrusu, denklem, R² ve y = x köşegeni. Tek gende yalnızca nokta çizilir.

```python
ax = eff.plot_neutrality(file_name="notr.png")
```

### EffectiveCodons.plot_pr2

`EffectiveCodons.plot_pr2(codons=None, ax=None, labels=False, groups=None, file_name=None)`

PR2 grafiği, merkezde 0,5 çizgileriyle. `codons` `pr2` ile aynıdır.

```python
ax = eff.plot_pr2(codons=None, file_name="pr2.png")
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_finish` | Noktaları çizer, istenirse sembolleri yazar, kaydeder |
| `_linear_regression` | En küçük kareler doğrusu ve r²; iki noktadan azsa `nan` |
| `EffectiveCodons.__repr__` | `EffectiveCodons(CUB: ...)` |
