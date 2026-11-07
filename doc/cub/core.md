# `crn.cub.core`: A. Ham sayımlar

Kaynak: [`crn/cub/core.py`](../../crn/cub/core.py)

`CUB` genleri hesaba hazırlar ve A bölümünün ham sayımlarını verir; öbür bölümler bir `CUB` alıp bu sayımları kullanır.

Girdi bir `Gene` listesidir (ya da `.crn` yolları). Her gende tam olarak bir genomik DNA, bir transkript ve bir protein bulunmalıdır; çok izoformlu genler `Gene.isoforms()` ile ayrılır. Kodonlar yalnızca okuma çerçevesinde anlamlı olduğundan, transkriptteki çevirisi proteinle birebir aynı ORF bulunur ve hesap bu CDS üzerinden yapılır. Verilen genler değiştirilmez; `cub.genes` transkripti CDS'e indirgenmiş kopyaları verir.

Stop kodonları `codon_counts` içinde sayılır ama L'ye, amino asit sayılarına ve frekanslara girmez. Belirsiz baz içeren kodonlar hiçbir sayıma girmez.

<!-- toc -->
- [1. Uyarı sınıfı](#1-uyarı-sınıfı): [`CubWarning`](#cubwarning)
- [2. Kurulum](#2-kurulum): [`CUB`](#cub) · [`genes`](#cubgenes) · [`names`](#cubnames) · [`genetic_code`](#cubgenetic_code)
- [3. Dizi erişimi](#3-dizi-erişimi): [`cds`](#cubcds) · [`codons`](#cubcodons)
- [4. Kodon kümeleri](#4-kodon-kümeleri): [`sense_codons`](#cubsense_codons) · [`codon_families`](#cubcodon_families) · [`synonymous_codons`](#cubsynonymous_codons) · [`fourfold_codons`](#cubfourfold_codons)
- [5. A bölümü: ölçümler](#5-a-bölümü-ölçümler): [`codon_counts`](#cubcodon_counts) · [`total_codons`](#cubtotal_codons) · [`amino_acid_counts`](#cubamino_acid_counts) · [`codon_freqs`](#cubcodon_freqs) · [`codon_per_thousand`](#cubcodon_per_thousand)
- [6. Yapı bilgileri](#6-yapı-bilgileri): [`termini`](#cubtermini) · [`chain_count`](#cubchain_count) · [`codon_pair_counts`](#cubcodon_pair_counts) · [`table`](#cubtable)
- [7. İç yardımcılar](#7-iç-yardımcılar)
<!-- /toc -->

## 1. Uyarı sınıfı

### CubWarning

Bir verinin hesap dışında bırakıldığını bildirir, `UserWarning` alt sınıfıdır. `Divergence.rcdi` referansta olmayan kodonları, `CodonContext.cpb` referansta olmayan çiftleri böyle bildirir. Susturmak için `warnings.filterwarnings("ignore", category=CubWarning)`.

## 2. Kurulum

### CUB

`CUB(genes, genetic_code=GeneticCode.STANDARD)`

`genes` `Gene` nesneleri ya da `.crn` gen dosyası yolları listesidir. `genetic_code` ORF arama ve kodon çevirisi için kullanılır. `genes` liste değilse ya da elemanı geçersizse `TypeError`; boşsa, iki gen aynı adlıysa, bir gen "bir DNA, bir transkript, bir protein" kuralına uymuyorsa ya da proteinle örtüşen ORF yoksa `ValueError` verilir. Hata mesajı sorunlu genleri tek tek listeler.

İki `CUB` `+` ile birleştirilebilir (`CUB.__add__`). Sayımlar yeniden yapılmaz; genetik kodlar farklıysa ya da gen adları çakışıyorsa `ValueError` verir.

```python
from crn.fetch import Gencode
from crn.cub import CUB
from crn import GeneticCode

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

a = CUB(kur(["TP53", "HBB"]), genetic_code=GeneticCode.STANDARD)
b = CUB(kur(["ACTB"]))
cub = a + b
print(cub, len(cub))  # -> CUB: 3 gen, L=915 3
```

Genler `.crn` olarak saklanırsa ikinci çalıştırmada GENCODE'a gerek kalmaz.

```python
hbb = kur(["HBB"])[0]
hbb.to_crn("hbb.crn")
print(CUB(["hbb.crn"]).names)  # -> ['HBB_ENST00000335295.4']
```

### CUB.genes

`CUB.genes  (özellik)`

Hesaba giren genleri verir: transkriptleri CDS'e indirgenmiş kopyalar.

```python
print(cub.genes[0])  # -> TP53_ENST00000269305.9: DNA <25768>, 1 transkript, 1 protein
```

### CUB.names

`CUB.names  (özellik)`

Gen adlarını ekleme sırasıyla verir.

```python
ad = cub.names[0]
print(ad)  # -> TP53_ENST00000269305.9
```

### CUB.genetic_code

`CUB.genetic_code  (özellik)`

Kullanılan genetik kodu verir.

```python
print(cub.genetic_code)  # -> GeneticCode.STANDARD
```

## 3. Dizi erişimi

### CUB.cds

`CUB.cds(gene)`

`gene` geninin hesaba giren CDS dizisini verir (RNA, 5'→3', stop dahil).

```python
print(cub.cds(ad)[:9], cub.cds(ad)[-3:])  # -> AUGGAGGAG UGA
```

### CUB.codons

`CUB.codons(gene=None)`

Kodon listesini okuma sırasıyla verir, belirsiz bazlılar dahil. `gene` `None` ise genler art arda gelir.

```python
print(cub.codons(ad)[:3])  # -> ['AUG', 'GAG', 'GAG']
```

## 4. Kodon kümeleri

### CUB.sense_codons

`CUB.sense_codons(codons=CODONS)`

`codons` içinden amino asit kodlayanları verir.

```python
print(len(cub.sense_codons()), cub.sense_codons(["AUG", "UAA"]))  # -> 61 ['AUG']
```

### CUB.codon_families

`CUB.codon_families(amino_acids=AMINO_ACIDS)`

`amino_acids` için amino asitten kodon ailesine giden sözlüğü verir.

```python
print(cub.codon_families(amino_acids="F"))  # -> {'F': ['UUU', 'UUC']}
```

### CUB.synonymous_codons

`CUB.synonymous_codons(codons=CODONS)`

`codons` içinden eş anlamlı olanları verir; tek kodonlu amino asitler ve stop hariç.

```python
print(cub.synonymous_codons(["AUG", "GCC", "UAA"]))  # -> ['GCC']
```

### CUB.fourfold_codons

`CUB.fourfold_codons(codons=CODONS)`

`codons` içinden dört katlı dejenere olanları verir.

```python
print(cub.fourfold_codons(["GCC", "AAA"]))  # -> ['GCC']
```

## 5. A bölümü: ölçümler

### CUB.codon_counts

`CUB.codon_counts(gene=None, codons=CODONS)`

Madde 1. `codons` içindeki her kodonun kaç kez geçtiğini verir, stoplar dahil.

```python
print(cub.codon_counts(ad, codons=["AUG", "UGA"]))  # -> {'AUG': 12, 'UGA': 1}
```

### CUB.total_codons

`CUB.total_codons(gene=None)`

Madde 2. Toplam kodon sayısı L, stoplar hariç.

```python
print(cub.total_codons(ad), cub.total_codons())  # -> 393 915
```

### CUB.amino_acid_counts

`CUB.amino_acid_counts(gene=None, amino_acids=AMINO_ACIDS)`

Madde 3. Aynı amino asidi kodlayan kodonların toplamı, stop hariç. `amino_acids` sayılacak harflerdir.

```python
print(cub.amino_acid_counts(ad, amino_acids="MW"))  # -> {'M': 12, 'W': 4}
```

### CUB.codon_freqs

`CUB.codon_freqs(gene=None)`

Madde 4. Kodon frekansı N_kodon / L; yalnızca amino asit kodlayan kodonlar, toplamı 1.

```python
print(round(cub.codon_freqs(ad)["AUG"], 4))  # -> 0.0305
```

### CUB.codon_per_thousand

`CUB.codon_per_thousand(gene=None)`

Madde 5. Binde frekans: kodon frekansı × 1000. Kodon kullanım tabloları genellikle bu birimi kullanır.

```python
print(round(cub.codon_per_thousand(ad)["AUG"], 2))  # -> 30.53
```

## 6. Yapı bilgileri

### CUB.termini

`CUB.termini(gene=None)`

Her proteinin (ilk, son) kalıntısını verir; son kalıntı stop'tan öncekidir.

```python
print(cub.termini(ad))  # -> [('M', 'D')]
```

### CUB.chain_count

`CUB.chain_count(gene=None)`

Hesaba giren protein zinciri sayısı: `gene` `None` ise gen sayısı, değilse 1.

```python
print(cub.chain_count(), cub.chain_count(ad))  # -> 3 1
```

### CUB.codon_pair_counts

`CUB.codon_pair_counts(gene=None)`

Komşu kodon çiftlerinin sayılarını 6 harflik anahtarla verir. Stop ya da belirsiz baz içeren çiftler sayılmaz, sonuç önbelleğe alınır.

```python
ciftler = cub.codon_pair_counts(ad)
print(sum(ciftler.values()), len(next(iter(ciftler))))  # -> 392 6
```

### CUB.table

`CUB.table(**kwargs)`

Her gen için bütün skaler ölçümleri tek tabloda verir. `**kwargs` [`summary_table`](report.md#summary_table) parametreleridir.

```python
satir = cub.table(protein=False)[0]
print(satir["gene"] == ad, round(satir["GC3s"], 3))  # -> True 0.61
```

## 7. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_load_gene` | Liste elemanını `Gene`'e çevirir: `Gene` ise aynen, yol ise `Gene.from_crn` |
| `_gene_problems` | Genin "bir DNA, bir transkript, bir protein" kuralına uymayan yanlarını listeler |
| `CUB._finish` | Havuz sayımını ve kodon çifti önbelleğini kurar |
| `CUB._fit_to_protein` | Transkripti proteinle örtüşen ORF'ye indirgenmiş gen kopyası üretir; `RNA.match_orf` kullanır |
| `CUB._split` | Diziyi üçerli kodonlara böler |
| `CUB._check` | Gen adını doğrular |
| `CUB._selected` | `None` ise bütün adlar, değilse yalnızca o ad |
| `CUB._counter` | Havuz ya da genin `Counter` nesnesi |
| `CUB.__len__` | Gen sayısı |
| `CUB.__repr__` | `CUB: 3 gen, L=915` |
