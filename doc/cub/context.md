# `crn.cub.context`: H. tRNA ve kodon bağlamı

Kaynak: [`crn/cub/context.py`](../../crn/cub/context.py)

Kodonları tRNA havuzuna ve komşu kodonlarına göre değerlendirir.

tAI (dos Reis ve ark. 2004): tGCN, kodonu okuyabilen antikodonun tRNA gen kopya sayısıdır. `trna_counts` ile verilir ve GtRNAdb FASTA dosyasından `read_trna_counts` ile okunabilir. Kodonun 3. bazını okuyan antikodonlar ve wobble cezaları [`TAI_WOBBLE_PAIRS` ve `TAI_S_VALUES`](../tables.md#7-tai-sabitleri) içindedir. `kingdom="prokaryote"` ise AUA kodonu lizidinli tRNA-Ile CAU'dan da okunur. Met kodonları başlatıcı tRNA ile karışmasın diye, stoplar da doğal olarak hesaba girmez. W'si 0 olan kodonların w'si öteki w'lerin geometrik ortalaması alınır.

Dinükleotit ρ (Karlin ve Burge 1995): `frame` `None` ise kodlayan dizinin tamamı kullanılır. `"12"` ve `"23"` kodon içi konum çiftleri, `"31"` kodon sınırını aşan çifttir (bir kodonun 3. ve sonrakinin 1. bazı); bu durumlarda f(X) ve f(Y) o konumların baz oranlarıdır. ρ < 0,78 az, > 1,23 çok temsil edilmiş kabul edilir.

CPS ve CPB (Coleman ve ark. 2008): sayımlar referans gen setinden alınır, genellikle genomun bütün ORF'leri. A ve B komşu kodonlar, X ve Y onların amino asitleridir. Pozitif CPS çiftin beklenenden çok kullanıldığını gösterir, CPB genin kodon çiftlerinin CPS ortalamasıdır.

`gene` parametresi [ortak kurala](README.md#1-ortak-kurallar) uyar.

<!-- toc -->
- [1. tRNA sayımı](#1-trna-sayımı): [`read_trna_counts`](#read_trna_counts)
- [2. Sınıf](#2-sınıf): [`CodonContext`](#codoncontext) · [`trna_weights`](#codoncontexttrna_weights) · [`tai`](#codoncontexttai) · [`dinucleotide_counts`](#codoncontextdinucleotide_counts) · [`rho`](#codoncontextrho) · [`codon_pair_scores`](#codoncontextcodon_pair_scores) · [`cpb`](#codoncontextcpb)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. tRNA sayımı

`TAI_KINGDOMS` kabul edilen değerlerdir: `eukaryote`, `prokaryote`.

### read_trna_counts

`read_trna_counts(file_name)`

`file_name` GtRNAdb FASTA dosyasından antikodon başına tRNA gen kopya sayısını verir. Antikodonlar 5'→3' RNA alfabesiyle döner. Belirsiz, baskılayıcı, selenosistein ve yalancı gen kayıtları (`_TRNA_SKIPPED`) sayılmaz.

```python
from crn.fetch import Gencode
from crn.fetch._store import DEFAULT_SYS_DIR
from crn.cub import CUB, CodonContext, read_trna_counts

trna = read_trna_counts(DEFAULT_SYS_DIR / "gtrnadb" / "hg38-tRNAs.fa")
print(len(trna), trna["UUC"])  # -> 47 8
```

## 2. Sınıf

### CodonContext

`CodonContext(cub, reference=None)`

`cub` ölçülecek set, `reference` CPS'nin hesaplandığı settir, verilmezse `cub`. Bir `CUB` verilmezse `TypeError` verir. `CodonContext.composition` içteki `Composition` nesnesidir.

```python
gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "HBB"]))
referans = CUB(kur(["ACTB", "GAPDH", "TUBB", "RPL3", "EEF2", "MEFV"]))
ad = cub.names[0]
ctx = CodonContext(cub, reference=referans)
```

### CodonContext.trna_weights

`CodonContext.trna_weights(trna_counts, s_values=TAI_S_VALUES, kingdom='eukaryote')`

Madde 32. Her kodonun tAI ağırlığı: w = W / W_en büyük, W = Σ (1 − s)·tGCN. `trna_counts` antikodon → kopya sayısıdır (5'→3', T ya da U). `s_values` wobble cezaları, `kingdom` alemdir. Met ve stop dışındaki kodonlar döner. `kingdom` geçersizse ya da hiçbir kodon için tRNA yoksa `ValueError` verir.

```python
w = ctx.trna_weights(trna, kingdom="eukaryote")
print(len(w), max(w.values()))  # -> 60 1.0
```

### CodonContext.tai

`CodonContext.tai(gene=None, trna_counts=None, s_values=TAI_S_VALUES, kingdom='eukaryote')`

Madde 32. Genin kodonlarının tAI ağırlıklarının geometrik ortalaması. `trna_counts` zorunludur, verilmezse `ValueError`. `s_values` ve `kingdom` `trna_weights` ile aynıdır.

```python
print(round(ctx.tai(ad, trna_counts=trna), 4))  # -> 0.327
```

### CodonContext.dinucleotide_counts

`CodonContext.dinucleotide_counts(gene=None, frame=None, bases=RNA_BASES)`

Madde 33. Kodlayan dizideki dinükleotit sayıları, stoplar hariç. `frame` `None`, `"12"`, `"23"` ya da `"31"`, başka değer `ValueError` verir. `bases` kullanılacak bazlardır.

```python
print(ctx.dinucleotide_counts(ad, frame="31", bases="CG")["CG"])  # -> 13
```

### CodonContext.rho

`CodonContext.rho(gene=None, frame=None, bases=RNA_BASES)`

Madde 33. Göreli bolluk ρ(XY) = f(XY) / (f(X)·f(Y)); 1 beklenen düzeydir. Bazlardan biri hiç geçmiyorsa `nan`. `frame` ve `bases` önceki metotla aynıdır.

```python
print(round(ctx.rho(ad)["CG"], 3), round(ctx.rho(ad, frame="31", bases="ACGU")["CG"], 3))  # -> 0.444 0.34
```

### CodonContext.codon_pair_scores

`CodonContext.codon_pair_scores()`

Madde 34. Referansta geçen her kodon çiftinin CPS'si: ln [N(AB) / (N(A)·N(B) / (N(X)·N(Y)) × N(XY))]. Sonuç önbelleğe alınır.

```python
print(len(ctx.codon_pair_scores()) > 1000)  # -> True
```

### CodonContext.cpb

`CodonContext.cpb(gene=None)`

Madde 34. Genin kodon çiftlerinin CPS ortalaması. Referansta geçmeyen çiftler çıkarılır ve `CubWarning` ile bildirilir, hesaba girecek çift yoksa `nan` döner. Küçük referansla çok çift düşer, bu yüzden gerçek analizde bütün genom kullanın.

```python
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cpb = ctx.cpb(ad)
print(round(cpb, 4))  # -> 0.3499
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `DINUCLEOTIDE_FRAMES` | `frame` değerini iki konuma çevirir |
| `_TRNA_HEADER_PATTERNS` | GtRNAdb başlıklarından amino asit ve antikodon çıkaran iki kalıp |
| `_anticodon` | Kodonun 1. ve 2. bazını eşleyen, 34. bazı verilen antikodon |
| `CodonContext._pair_scores` | CPS hesabının gövdesi; bir kez hesaplanır |
| `CodonContext.__repr__` | Ölçülen ve referans seti gösterir |
