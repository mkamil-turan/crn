# `crn.cub.adaptation`: D. Referans setine dayalı indeksler

Kaynak: [`crn/cub/adaptation.py`](../../crn/cub/adaptation.py)

Bir genin kodonlarını referans gen setine göre değerlendirir. Referans, karşılaştırma seti ve optimal kodon seçimi [`Synonymous`](synonymous.md) ile aynıdır.

w (Sharp ve Li, 1987) yalnızca eş anlamlı kodonlar için tanımlıdır: bir kodonun referanstaki sayısının ailesinin en büyük sayısına oranıdır. Tek kodonlu amino asitler hep w = 1 alacakları için CAI'ye girmez, stoplar da girmez; bu yüzden CAI'deki L genin eş anlamlı kodon sayısıdır. Referansta hiç geçmeyen kodonun w'si 0 olur ve ln 0 tanımsızdır, bu nedenle böyle kodonların sayısı `absent_count` kabul edilir. Referansta hiç kullanılmayan bir amino asidin bütün kodonları böylece w = 1 alır.

CBI (Bennetzen ve Hall, 1982) yalnızca en az bir optimal kodonu olan amino asitlerin kodonlarını sayar: (N_opt − N_rastgele) / (N_toplam − N_rastgele). N_rastgele, kodonlar rastgele seçilseydi beklenecek optimal kodon sayısıdır. 1 yalnızca optimal kodon, 0 rastgele kullanım demektir; optimal kodonlar rastgeleden az kullanılırsa negatif olur.

`gene`, `ax`, `file_name` ve `groups` parametreleri [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`Adaptation`](#adaptation) · [`relative_adaptiveness`](#adaptationrelative_adaptiveness) · [`cai`](#adaptationcai) · [`cbi`](#adaptationcbi) · [`cai_points`](#adaptationcai_points) · [`plot_cai`](#adaptationplot_cai)
<!-- /toc -->

## 1. Sınıf

`ABSENT_COUNT` (0,5) `absent_count` parametresinin varsayılanıdır.

### Adaptation

`Adaptation(cub, reference=None, background=None)`

`cub` ölçülecek set, `reference` w'lerin ve optimal kodonların belirleneceği set, `background` `method="chi2"` için karşılaştırma setidir. Verilmeyen değerler `cub` olur. Hatalar `Synonymous` ile aynıdır. `Adaptation.synonymous` içteki `Synonymous` nesnesidir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, Adaptation

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "HBB", "ACTB", "MEFV"]))
yuksek = CUB(kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2"]))
dusuk = CUB(kur(["OR1A1", "OR2A2", "OR5A1", "OR10A2"]))
ad = cub.names[0]
ada = Adaptation(cub, reference=yuksek, background=dusuk)
```

### Adaptation.relative_adaptiveness

`Adaptation.relative_adaptiveness(absent_count=ABSENT_COUNT)`

Madde 18. w: referanstaki kodon sayısı / ailesinin en büyük sayısı. Yalnızca eş anlamlı kodonlar için döner, sonuç önbelleğe alınır. `absent_count` referansta hiç geçmeyen kodonun yerine kullanılacak sayıdır; 0 ya da negatifse `ValueError` verir.

```python
w = ada.relative_adaptiveness()
print(len(w), max(w.values()))  # -> 59 1.0
print(round(ada.relative_adaptiveness(absent_count=0.1)["GCC"], 3))  # -> 1.0
```

### Adaptation.cai

`Adaptation.cai(gene=None, w=None, absent_count=ABSENT_COUNT)`

Madde 19. Genin eş anlamlı kodonlarının w değerlerinin geometrik ortalaması, 0-1 arası. `w` kodon sözlüğünü doğrudan verir, verilmezse `relative_adaptiveness(absent_count)` kullanılır. Genin eş anlamlı kodonu yoksa 0 döner. `w` içinde 0 ya da negatif değer varsa `ValueError` verir.

```python
print(round(ada.cai(ad), 3), round(ada.cai(ad, w=w), 3), round(ada.cai(ad, absent_count=0.1), 3))  # -> 0.677 0.677 0.677
```

### Adaptation.cbi

`Adaptation.cbi(gene=None, optimal_codons=None, method='rscu', threshold=PREFERRED_RSCU, alpha=OPTIMAL_ALPHA)`

Madde 20. (N_opt − N_rastgele) / (N_toplam − N_rastgele), payda 0 ise 0. `optimal_codons` listeyi doğrudan verir, verilmezse `Synonymous.optimal_codons(method, threshold, alpha)` kullanılır.

```python
print(round(ada.cbi(ad), 3), round(ada.cbi(ad, method="chi2", alpha=0.01), 3), round(ada.cbi(ad, threshold=1.3), 3))  # -> 0.279 0.363 0.281
print(round(ada.cbi(ad, optimal_codons=["GCC", "CUG"]), 3))  # -> 0.351
```

### Adaptation.cai_points

`Adaptation.cai_points()`

Her gen için CAI'yi `[{"gene", "cai"}, ...]` biçiminde verir.

```python
print(sorted(ada.cai_points()[0]))  # -> ['cai', 'gene']
```

### Adaptation.plot_cai

`Adaptation.plot_cai(ax=None, bins=20, groups=None, file_name=None)`

CAI frekans dağılımı. Yüksek CAI referans genlere benzer kodon kullanımı gösterir, kesikli çizgi ortalamadır. `bins` çubuk sayısıdır. `groups` verilirse çubuklar gruplara göre renkli ve yığılmış çizilir.

```python
ax = ada.plot_cai(bins=5, file_name="cai.png")
```

`Adaptation.__repr__` içteki `Synonymous` gösterimini sarar.
