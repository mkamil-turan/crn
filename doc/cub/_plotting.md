# `crn.cub._plotting`: ortak grafik yardımcıları

Kaynak: [`crn/cub/_plotting.py`](../../crn/cub/_plotting.py)

Bütün `plot_*` metotları aynı biçimi buradan alır, yani renkler ve stiller tek yerden değişir. matplotlib yalnızca `axes` çağrılınca yüklenir, bu yüzden grafik çizmeyen analizler onu hiç içe aktarmaz.

<!-- toc -->
- [1. `groups` parametresi](#1-groups-parametresi)
- [2. Biçim sabitleri](#2-biçim-sabitleri)
- [3. Fonksiyonlar](#3-fonksiyonlar): [`axes`](#axes) · [`save`](#save) · [`gene_symbols`](#gene_symbols) · [`gene_groups`](#gene_groups) · [`group_colors`](#group_colors) · [`group_legend`](#group_legend) · [`scatter`](#scatter) · [`histogram`](#histogram)
<!-- /toc -->

## 1. `groups` parametresi

Gen başına nokta ya da değer çizen grafikler `groups` alır. Her grup ayrı renkte çizilir, gösterge grafiğin sağına konur. `Groups` bu parametrenin tipidir ve üç biçim kabul eder: metin (genin features'ındaki anahtar, ör. `"organism"`), sözlük (gen adı → grup) ya da fonksiyon (`Gene` alır, grup adı döner). Grubu bulunamayan gen `ValueError` verir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, EffectiveCodons
from crn.cub import _plotting as p

gc = Gencode(verbose=False)
cub = CUB([gc.canonical(s).isoforms()[0] for s in ["TP53", "HBB", "ACTB", "RPL3"]])

gruplar = {ad: ("ribozom" if ad.startswith("RPL") else "diger") for ad in cub.names}
ax = EffectiveCodons(cub).plot_enc(groups=gruplar, labels=True, file_name="enc.png")
ax = EffectiveCodons(cub).plot_pr2(groups=lambda gen: gen.features.first("gene_type"), file_name="pr2.png")
```

## 2. Biçim sabitleri

| Ad | Değer |
|---|---|
| `PLOT_POINT_STYLE` | Nokta stili (`scatter` argümanları) |
| `PLOT_CURVE_STYLE` | Eğri ve regresyon doğrusu stili |
| `PLOT_GUIDE_STYLE` | Kesikli kılavuz çizgileri |
| `PLOT_LABEL_SIZE` | Eksen etiketi yazı boyutu: 7 |
| `PLOT_GENE_LABEL_SIZE` | Gen sembolü yazı boyutu: 4 |
| `PLOT_LEGEND_SIZE` | Gösterge yazı boyutu: 7 |
| `PLOT_GROUP_COLORS` | Grup renkleri; ilk görünüş sırasıyla atanır |

## 3. Fonksiyonlar

### axes

`axes(ax, figsize=(5, 4.5))`

`ax` verilmişse onu, yoksa `figsize` boyutunda yeni şeklin `Axes`'ini verir.

```python
ax = p.axes(None, figsize=(4, 3))
print(type(ax).__name__)  # -> Axes
```

### save

`save(ax, file_name)`

`file_name` verilmişse şekli 300 dpi kaydeder ve `ax` nesnesini döndürür.

```python
print(p.save(ax, None) is ax)  # -> True
```

### gene_symbols

`gene_symbols(cub, names)`

`names` genlerinin sembollerini verir: features'taki `gene_name`, yoksa genin adı. `cub` kaynak settir.

```python
print(p.gene_symbols(cub, cub.names[:2]))  # -> ['TP53', 'HBB']
```

### gene_groups

`gene_groups(cub, names, groups)`

`names` genlerinin grubunu verir, `groups` `None` ise `None`. `cub` kaynak settir. Grup bulunamazsa `ValueError` verir.

```python
print(p.gene_groups(cub, cub.names[:1], gruplar))  # -> ['diger']
```

### group_colors

`group_colors(groups)`

`groups` listesindeki gruplara ilk görünüş sırasıyla renk atar.

```python
print(list(p.group_colors(["a", "b", "a"])))  # -> ['a', 'b']
```

### group_legend

`group_legend(ax, handles, title=None)`

Grup göstergesini `ax` grafiğinin sağına koyar, var olan göstergeyi korur. `handles` gösterge nesneleri, `title` başlıktır.

```python
ax.scatter([1], [1], label="a")
p.group_legend(ax, ax.get_legend_handles_labels()[0], title="grup")
```

### scatter

`scatter(ax, cub, points, groups=None, labels=False)`

Gen noktalarını (`points`: x, y, gen adı) `ax` üzerine çizer. `cub` sembol ve grup bilgisi için gereklidir, `groups` renkleri, `labels` sembolleri belirler.

```python
ax = p.axes(None)
p.scatter(ax, cub, [(0.1, 0.2, cub.names[0])], groups=None, labels=True)
```

### histogram

`histogram(ax, cub, values, bins, groups=None, **kwargs)`

(gen adı, değer) çiftlerini (`values`) histogram çizer. `ax` hedef eksen, `cub` kaynak set, `bins` çubuk tanımıdır. `groups` verilirse çubuklar renkli ve yığılmış olur. `**kwargs` `Axes.hist`'e gider.

```python
ax = p.axes(None)
p.histogram(ax, cub, [(ad, 0.1) for ad in cub.names], bins=3, groups=None)
```
