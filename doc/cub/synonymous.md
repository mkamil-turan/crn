# `crn.cub.synonymous`: C. Sinonim kodon tercihi

Kaynak: [`crn/cub/synonymous.py`](../../crn/cub/synonymous.py)

Bir amino asidin kodonları arasındaki tercihi ölçer. RSCU'da 1 tarafsız kullanımı, 1'in üstü tercihi gösterir; bir ailenin RSCU değerlerinin toplamı aile büyüklüğüne eşittir. Hiç kullanılmayan amino asidin kodonlarında kodon oranı ve RSCU 0'dır. Tek kodonlu amino asitlerin (Met, Trp) RSCU'su kullanıldıklarında hep 1'dir, bu yüzden tercih listelerine girmezler.

Fop için optimal kodonlar bir referans gen setinden gelir, genellikle yüksek ifade edilen genlerden kurulan ikinci bir `CUB`. Referans verilmezse `cub`'ın kendisi kullanılır ve Fop genin tercihini kendisiyle ölçmüş olur, bu tek gende anlamsızdır.

`gene`, `ax` ve `file_name` parametreleri [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sabitler](#1-sabitler)
- [2. Sınıf](#2-sınıf): [`Synonymous`](#synonymous) · [`codon_ratio`](#synonymouscodon_ratio) · [`rscu`](#synonymousrscu) · [`preferred_codons`](#synonymouspreferred_codons) · [`underused_codons`](#synonymousunderused_codons) · [`optimal_codons`](#synonymousoptimal_codons) · [`fop`](#synonymousfop)
- [3. Grafikler](#3-grafikler): [`plot_rscu`](#synonymousplot_rscu) · [`plot_rscu_heatmap`](#synonymousplot_rscu_heatmap)
- [4. Modül düzeyinde denetimler](#4-modül-düzeyinde-denetimler): [`check_choice`](#check_choice) · [`check_cub`](#check_cub)
- [5. İç yardımcılar](#5-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler

`PREFERRED_RSCU` (1,6) ve `UNDERUSED_RSCU` (0,6) tercih eşikleridir. `OPTIMAL_ALPHA` (0,01) χ² yöntemindeki anlamlılık eşiğidir. `OPTIMAL_METHODS` (`rscu`, `chi2`) ve `FOP_DENOMINATORS` (`synonymous`, `optimal`) kabul edilen değerlerdir.

## 2. Sınıf

### Synonymous

`Synonymous(cub, reference=None, background=None)`

`cub` ölçülecek gen setidir. `reference` optimal kodonların belirleneceği set, `background` `method="chi2"` iken referansla karşılaştırılan settir (ör. düşük ifadeli genler); ikisi de verilmezse `cub` kullanılır. Bir `CUB` verilmezse `TypeError`, genetik kodlar farklıysa `ValueError` verir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, Synonymous

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "HBB", "ACTB"]))
yuksek = CUB(kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2"]))
dusuk = CUB(kur(["OR1A1", "OR2A2", "OR5A1", "OR10A2", "OR4D1", "OR6A2"]))
ad = cub.names[0]
syn = Synonymous(cub, reference=yuksek, background=dusuk)
```

### Synonymous.codon_ratio

`Synonymous.codon_ratio(gene=None)`

Madde 14. Kodonun kendi amino asidi içindeki payı, N_kodon / N_amino asit. Bir ailenin oranlarının toplamı 1'dir.

```python
print(round(syn.codon_ratio(ad)["GCC"], 3))  # -> 0.5
```

### Synonymous.rscu

`Synonymous.rscu(gene=None)`

Madde 15. N_kodon / (N_amino asit / aile büyüklüğü), yani kodon oranı × aile büyüklüğü.

```python
print(round(syn.rscu(ad)["GCC"], 3))  # -> 2.0
```

### Synonymous.preferred_codons

`Synonymous.preferred_codons(gene=None, threshold=PREFERRED_RSCU)`

Madde 16. RSCU'su `threshold`'u aşan eş anlamlı kodonlar.

```python
print(syn.preferred_codons(ad)[:5], len(syn.preferred_codons(ad, threshold=2.0)))  # -> ['UCC', 'CUG', 'CGC', 'AUC', 'ACC'] 4
```

### Synonymous.underused_codons

`Synonymous.underused_codons(gene=None, threshold=UNDERUSED_RSCU)`

Madde 16. RSCU'su `threshold`'un altındaki eş anlamlı kodonlar. Hiç kullanılmayan amino asitlerin kodonları dahil edilmez.

```python
print(len(syn.underused_codons(ad)), len(syn.underused_codons(ad, threshold=0.3)))  # -> 13 5
```

### Synonymous.optimal_codons

`Synonymous.optimal_codons(method='rscu', threshold=PREFERRED_RSCU, alpha=OPTIMAL_ALPHA)`

Fop ve CBI için optimal kodonları `CODONS` sırasıyla verir, sonuç önbelleğe alınır. `method="rscu"` referansın RSCU'su `threshold`'u aşan kodonlardır. `method="chi2"` (Ikemura; CodonW) her kodon için referans ve background'da [kodon, ailenin öteki kodonları] sayılarından 2×2 χ² kurar; payı referansta yüksek ve p < `alpha` olanları alır. `method` geçersizse ya da `"chi2"` için reference ile background aynı setse `ValueError` verir.

```python
print(len(syn.optimal_codons()), len(syn.optimal_codons(method="chi2", alpha=0.001)))  # -> 9 1
```

### Synonymous.fop

`Synonymous.fop(gene=None, optimal_codons=None, method='rscu', denominator='synonymous', threshold=PREFERRED_RSCU, alpha=OPTIMAL_ALPHA)`

Madde 17. Optimal kodon sayısı / payda (0-1; payda 0 ise 0). `optimal_codons` listeyi doğrudan verir, verilmezse `method`, `threshold` ve `alpha` ile hesaplanır. `denominator="synonymous"` bütün eş anlamlı kodonları, `"optimal"` (CodonW) yalnızca en az bir optimal kodonu olan amino asitlerin kodonlarını payda alır. Geçersiz değer `ValueError` verir.

```python
print(round(syn.fop(ad), 3), round(syn.fop(ad, method="chi2", denominator="optimal"), 3))  # -> 0.244 0.345
print(round(syn.fop(ad, optimal_codons=["GCC"]), 4))  # -> 0.0318
```

## 3. Grafikler

### Synonymous.plot_rscu

`Synonymous.plot_rscu(gene=None, ax=None, file_name=None)`

RSCU çubuk grafiği: eş anlamlı kodonlar amino aside göre gruplanır, ardışık gruplar iki renkte boyanır, 1 çizgisi tarafsız kullanımı gösterir.

```python
ax = syn.plot_rscu(gene=ad, file_name="rscu.png")
print(type(ax).__name__)  # -> Axes
```

### Synonymous.plot_rscu_heatmap

`Synonymous.plot_rscu_heatmap(ax=None, genes=None, cmap='YlGnBu', groups=None, file_name=None)`

RSCU ısı haritası: satırlar genler, sütunlar eş anlamlı kodonlar. `genes` verilmezse bütün genler, `cmap` matplotlib renk haritasıdır. `groups` verilirse satırlar gruplara göre sıralanır, sembol renkleri grubu gösterir ve gruplar arasına ayraç çizilir.

```python
ax = syn.plot_rscu_heatmap(genes=cub.names[:2], cmap="viridis", file_name="isi.png")
```

## 4. Modül düzeyinde denetimler

### check_choice

`check_choice(name, value, choices)`

`value` `choices` içinde değilse `name` adıyla `ValueError` verir. Paketin başka modülleri de seçenek denetiminde bunu kullanır.

```python
from crn.cub.synonymous import check_choice, check_cub

check_choice("method", "rscu", ("rscu", "chi2"))
```

### check_cub

`check_cub(name, value, cub=None)`

`value` bir `CUB` değilse `name` adıyla `TypeError`, genetik kodu `cub`'ınkinden farklıysa `ValueError` verir.

```python
check_cub("reference", yuksek, cub=cub)
```

## 5. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_chi2_p` | 2×2 tablonun Pearson χ² p değeri (1 sd); boş satır ya da sütunda 1 |
| `Synonymous._find_optimal` | `optimal_codons`'ın hesap gövdesi; sonucu önbellek tutar |
| `Synonymous._synonymous_layout` | Eş anlamlı kodonları amino asit sırasıyla gruplar; grafikler ve COA bu sırayı kullanır |
| `Synonymous.__repr__` | Referans ve karşılaştırma setini gösterir |
