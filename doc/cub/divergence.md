# `crn.cub.divergence`: G. Entropi ve sapma tabanlı ölçümler

Kaynak: [`crn/cub/divergence.py`](../../crn/cub/divergence.py)

Kodon kullanımının düzensizliğini ve beklenen kullanımdan sapmasını ölçer.

SCUO (Wan ve ark. 2004): eş anlamlı her amino asit için kodonların Shannon entropisi H ve en büyük entropi Hmax = log2(aile büyüklüğü) hesaplanır, (Hmax − H) / Hmax amino asidin gendeki payıyla ağırlıklandırılıp toplanır. 0 yanlılık yok, 1 her amino asit için tek kodon demektir.

CDC ve RCBS'de beklenen kodon oranı konum baz oranlarının çarpımıdır, f1(x)·f2(y)·f3(z). Stoplar iki tarafta da yoktur, beklenen oranlar amino asit kodlayan kodonlar üzerinden toplamı 1 olacak biçimde normalize edilir. CDC'de Zhang ve ark. (2012) arka planı her konumun GC ve pürin içeriğiyle modeller ve bootstrap ile sınar; burada beklenen bileşim doğrudan konum baz oranlarından kurulur ve bootstrap yapılmaz, yani değerler özgün yazılımla birebir aynı olmayabilir. RCBS (Roymondal ve ark. 2009) (1 + d)'nin geometrik ortalamasıdır.

RCDI (Mueller ve ark. 2006) genin kodon oranının konaktakine oranını genin kodonlarıyla ağırlıklandırır. Yalnızca eş anlamlı kodonlar sayılır. 1 konakla aynı kullanım, büyük değer konağa uyumsuzluk demektir.

`gene` parametresi [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`Divergence`](#divergence) · [`expected_codon_freqs`](#divergenceexpected_codon_freqs) · [`scuo`](#divergencescuo) · [`cdc`](#divergencecdc) · [`rcdi`](#divergencercdi) · [`rcbs`](#divergencercbs)
- [2. İç yardımcılar](#2-iç-yardımcılar)
<!-- /toc -->

## 1. Sınıf

### Divergence

`Divergence(cub, reference=None)`

`cub` ölçülecek set, `reference` RCDI'deki konak setidir, verilmezse `cub`. Hatalar `Synonymous` ile aynıdır. `Divergence.synonymous` ve `Divergence.composition` içteki nesnelerdir.

```python
from crn.fetch import Gencode
from crn.cub import CUB, Divergence

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "HBB"]))
konak = CUB(kur(["ACTB", "GAPDH", "TUBB", "RPL3", "RPS6", "EEF2", "MEFV", "EEF1A1"]))
ad = cub.names[0]
div = Divergence(cub, reference=konak)
```

### Divergence.expected_codon_freqs

`Divergence.expected_codon_freqs(gene=None)`

Konum baz oranlarından beklenen kodon oranlarını verir.

```python
print(round(sum(div.expected_codon_freqs(ad).values()), 6))  # -> 1.0
```

### Divergence.scuo

`Divergence.scuo(gene=None)`

Madde 28. Amino asit payıyla ağırlıklı normalize entropi farkı.

```python
print(round(div.scuo(ad), 4))  # -> 0.1175
```

### Divergence.cdc

`Divergence.cdc(gene=None)`

Madde 29. Gözlenen ve beklenen kodon oranları arasındaki kosinüs uzaklığı (0-1).

```python
print(round(div.cdc(ad), 4))  # -> 0.1118
```

### Divergence.rcdi

`Divergence.rcdi(gene=None)`

Madde 30. Σ (oran_gen / oran_ref) × N_kodon / N, eş anlamlı kodonlarda. Referansta hiç geçmeyen kodonlar çıkarılır ve `CubWarning` ile bildirilir. Sayılabilecek kodon yoksa 0 döner.

```python
print(round(div.rcdi(ad), 4))  # -> 1.1295
```

### Divergence.rcbs

`Divergence.rcbs(gene=None)`

Madde 31. Genin kodonları üzerinden (1 + d)'nin geometrik ortalaması − 1. 0 yanlılık yok demektir.

```python
print(round(div.rcbs(ad), 4))  # -> 0.1731
```

## 2. İç yardımcılar

| Ad | Görevi |
|---|---|
| `Divergence._reference_ratio` | Referansın kodon oranları; bir kez hesaplanır |
| `Divergence.__repr__` | `Divergence(Synonymous(...))` |
