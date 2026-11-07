# `crn.cub.composition`: B. Nükleotit bileşimi

Kaynak: [`crn/cub/composition.py`](../../crn/cub/composition.py)

Kodon sayımlarından nükleotit bileşimini hesaplar. Bütün ölçümler L ile aynı kodonlar üzerindendir: stoplar ve belirsiz bazlı kodonlar sayılmaz. GC% yüzdedir, ötekiler 0-1 arası orandır. Diziler RNA olduğu için T yerine U yazılır. GC3 bütün kodonların üçüncü konumuna, GC3s yalnızca eş anlamlı kodonların üçüncü konumuna bakar; Met ve Trp'nin üçüncü bazı seçilemediği için kodon tercihi analizinde doğru ölçü GC3s'tir.

`gene` parametresi [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`Composition`](#composition) · [`nucleotide_counts`](#compositionnucleotide_counts) · [`nucleotide_freqs`](#compositionnucleotide_freqs)
- [2. B bölümü: ölçümler](#2-b-bölümü-ölçümler): [`gc_percent`](#compositiongc_percent) · [`gc1`](#compositiongc1) · [`gc2`](#compositiongc2) · [`gc12`](#compositiongc12) · [`gc3`](#compositiongc3) · [`gc3s`](#compositiongc3s) · [`base3s`](#compositionbase3s) · [`at3s`](#compositionat3s) · [`summary`](#compositionsummary)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sınıf

### Composition

`Composition(cub)`

`cub` bir `CUB` olmalıdır, değilse `TypeError` verilir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, Composition

gc = Gencode(verbose=False)
cub = CUB([gc.canonical(s).isoforms()[0] for s in ["TP53", "HBB", "ACTB"]])
ad = cub.names[0]
comp = Composition(cub)
print(comp)  # -> Composition(CUB: 3 gen, L=915)
```

### Composition.nucleotide_counts

`Composition.nucleotide_counts(gene=None, position=None, codons=CODONS, bases=RNA_BASES)`

`codons` içindeki amino asit kodlayan kodonlarda baz sayılarını verir. `position` `None` ise üç konum birden, 1, 2 ya da 3 ise o konum; başka değer `ValueError` verir. `bases` sayılacak bazlardır.

```python
print(comp.nucleotide_counts(ad, position=3))  # -> {'U': 86, 'C': 131, 'A': 61, 'G': 115}
print(comp.nucleotide_counts(ad, position=1, codons=["AUG"], bases="A"))  # -> {'A': 12}
```

### Composition.nucleotide_freqs

`Composition.nucleotide_freqs(gene=None, position=None, codons=CODONS, bases=RNA_BASES)`

Baz oranlarını verir, toplamları 1'dir. `position`, `codons` ve `bases` `nucleotide_counts` ile aynıdır.

```python
print(round(comp.nucleotide_freqs(ad, position=3)["G"], 3))  # -> 0.293
print(round(comp.nucleotide_freqs(ad, position=None, codons=["GCC", "GCG"], bases="GC")["C"], 3))  # -> 0.641
```

## 2. B bölümü: ölçümler

### Composition.gc_percent

`Composition.gc_percent(gene=None)`

Madde 6. GC%: (G + C) / toplam nükleotit × 100.

```python
print(round(comp.gc_percent(ad), 2))  # -> 56.91
```

### Composition.gc1

`Composition.gc1(gene=None)`

Madde 7. Kodonların 1. konumundaki G + C oranı.

```python
print(round(comp.gc1(ad), 3))  # -> 0.593
```

### Composition.gc2

`Composition.gc2(gene=None)`

Madde 8. Kodonların 2. konumundaki G + C oranı.

```python
print(round(comp.gc2(ad), 3))  # -> 0.489
```

### Composition.gc12

`Composition.gc12(gene=None)`

Madde 9. (GC1 + GC2) / 2. Nötrallik grafiğinin y eksenidir.

```python
print(round(comp.gc12(ad), 3))  # -> 0.541
```

### Composition.gc3

`Composition.gc3(gene=None)`

Madde 10. Kodonların 3. konumundaki G + C oranı.

```python
print(round(comp.gc3(ad), 3))  # -> 0.626
```

### Composition.gc3s

`Composition.gc3s(gene=None)`

Madde 11. Eş anlamlı kodonlarda 3. konum G + C oranı.

```python
print(round(comp.gc3s(ad), 3))  # -> 0.61
```

### Composition.base3s

`Composition.base3s(gene=None)`

Madde 12. A3s, U3s, G3s, C3s: eş anlamlı kodonlarda 3. konumdaki baz oranı. Toplamı 1'dir ve `G3s + C3s = GC3s`.

```python
b = comp.base3s(ad)
print(round(b["G"] + b["C"], 3))  # -> 0.61
```

### Composition.at3s

`Composition.at3s(gene=None)`

Madde 13. AT3s: 1 − GC3s.

```python
print(round(comp.at3s(ad), 3))  # -> 0.39
```

### Composition.summary

`Composition.summary(gene=None)`

6-13 arasındaki bütün ölçümleri tek sözlükte verir.

```python
print(list(comp.summary(ad)))  # -> ['GC%', 'GC1', 'GC2', 'GC12', 'GC3', 'GC3s', 'A3s', 'U3s', 'G3s', 'C3s', 'AT3s']
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_ratio` | Payda 0 ise 0 veren bölme |
| `Composition._gc` | Verilen konum ve kodonlarda G + C oranı |
| `Composition.__repr__` | `Composition(CUB: ...)` |
