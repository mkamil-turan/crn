# `crn.cub.ordination`: uyum çözümlemesi (COA)

Kaynak: [`crn/cub/ordination.py`](../../crn/cub/ordination.py)

Genlerin RSCU değerlerinden uyum çözümlemesiyle genleri birkaç ana eksene yerleştirir. Literatürde kodon kullanımının başlıca değişim eksenini (çoğu zaman GC3 ya da ifade düzeyi) bulmak için kullanılır, codonW'nin COA'sına karşılık gelir.

Gen × eş anlamlı kodon RSCU tablosu X için:

```text
P = X / ΣX,   satır kütlesi r_i = Σ_j P_ij,   sütun kütlesi c_j = Σ_i P_ij
S_ij = (P_ij − r_i c_j) / √(r_i c_j)            standartlaştırılmış artık
```

Özçözüm C = SᵀS matrisinden alınır, küçük bir saf Python Jacobi algoritmasıyla yapılır. Eksen k'nin açıkladığı atalet oranı λ_k / Σλ'dır, genlerin asıl koordinatları F = D_r^{-1/2} (S V). Hiçbir gende kullanılmayan kodonlar (sütun kütlesi 0) çıkarılır. Gerekli eksen sayısı min(n − 1, p − 1) ile sınırlıdır.

`ax`, `file_name`, `labels` ve `groups` parametreleri [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`CorrespondenceAnalysis`](#correspondenceanalysis) · [`inertia`](#correspondenceanalysisinertia) · [`coordinates`](#correspondenceanalysiscoordinates) · [`plot_coa`](#correspondenceanalysisplot_coa)
- [2. İç yardımcılar](#2-iç-yardımcılar)
<!-- /toc -->

## 1. Sınıf

### CorrespondenceAnalysis

`CorrespondenceAnalysis(cub)`

`cub` en az iki gen içermelidir. Bir `CUB` değilse `TypeError`; ikiden az gen varsa ya da RSCU tablosu boşsa `ValueError` verir. Kurulum sırasında çözümleme yapılır, sonuçlar `genes`, `codons` ve `eigenvalues` özniteliklerinde durur.

```python
from crn.fetch import Gencode
from crn.cub import CUB, CorrespondenceAnalysis

gc = Gencode(verbose=False)
cub = CUB([gc.canonical(s).isoforms()[0] for s in ["TP53", "HBB", "ACTB", "RPL3", "OR1A1", "MEFV"]])
coa = CorrespondenceAnalysis(cub)
print(len(coa.genes), len(coa.codons), len(coa.eigenvalues))  # -> 6 59 5
```

### CorrespondenceAnalysis.inertia

`CorrespondenceAnalysis.inertia()`

Her eksenin açıkladığı atalet oranını azalan sırayla verir, toplamı 1'dir.

```python
print(round(coa.inertia()[0], 3))  # -> 0.46
```

### CorrespondenceAnalysis.coordinates

`CorrespondenceAnalysis.coordinates(n_axes=2)`

Her gen için ilk `n_axes` eksendeki asıl koordinatı verir.

```python
print(sorted(coa.coordinates(n_axes=3)[0]))  # -> ['axis1', 'axis2', 'axis3', 'gene']
```

### CorrespondenceAnalysis.plot_coa

`CorrespondenceAnalysis.plot_coa(x_axis=1, y_axis=2, ax=None, labels=False, groups=None, file_name=None)`

COA saçılım grafiği, atalet oranları eksen adlarında yazılır. `x_axis` ve `y_axis` çizilecek eksenlerdir, 1 tabanlıdır.

```python
ax = coa.plot_coa(x_axis=1, y_axis=2, labels=True, file_name="coa.png")
```

## 2. İç yardımcılar

| Ad | Görevi |
|---|---|
| `JACOBI_MAX_SWEEPS` | Jacobi döngüsünün üst sınırı: 100 |
| `JACOBI_TOLERANCE` | Yakınsama eşiği: 1e-12 |
| `_jacobi_eigen` | Simetrik matrisin özdeğer ve özvektörleri, özdeğerler azalan sırada |
| `CorrespondenceAnalysis._fit` | Artık matrisini kurar, özçözümü alır, koordinatları hesaplar |
| `CorrespondenceAnalysis.__repr__` | Gen, kodon ve eksen sayısı |
