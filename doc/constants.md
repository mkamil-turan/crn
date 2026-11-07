# `crn.constants`: alfabeler, kodonlar ve genetik kodlar

Kaynak: [`crn/constants.py`](../crn/constants.py)

Dizi sınıflarının kabul ettiği harfler, tümleme tabloları ve genetik kodlar burada durur. Protein ve tRNA ölçümlerinin sayısal tabloları ayrı bir modüldedir: [`crn.tables`](tables.md).

<!-- toc -->
- [1. Sabitler](#1-sabitler)
- [2. Genetik kodlar](#2-genetik-kodlar): [`GeneticCode`](#geneticcode) · [`ncbi_id`](#geneticcodencbi_id) · [`table`](#geneticcodetable) · [`starts`](#geneticcodestarts) · [`stops`](#geneticcodestops) · [`sense_codons`](#geneticcodesense_codons) · [`families`](#geneticcodefamilies) · [`synonymous_codons`](#geneticcodesynonymous_codons) · [`fourfold_codons`](#geneticcodefourfold_codons)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler

Bütün sabitler `from crn.constants import ...` ile alınır.

| Ad | Değer | Kullanıldığı yer |
|---|---|---|
| `DNA_ALPHABET` | `ACGT`, IUPAC belirsizlik kodları (`RYKMSWBDHVN`) ve `-` | `DNA` alfabe denetimi |
| `RNA_ALPHABET` | `ACGU`, aynı belirsizlik kodları ve `-` | `RNA` alfabe denetimi |
| `PRT_ALPHABET` | 20 amino asit, `B Z J X U O`, `*` ve `-` | `PRT` alfabe denetimi |
| `GAP` | `"-"` | Boşluk simgesi; `seq_5to3` ve `find_orfs` bunu atlar |
| `STOP_SYMBOL` | `"*"` | Çeviride stop kodonunun simgesi |
| `DNA_BASES` | `"TCAG"` | Kodon tablolarının geleneksel baz sırası |
| `RNA_BASES` | `"UCAG"` | Aynı sıra, RNA alfabesiyle |
| `CODONS` | 64 kodonluk demet | Bütün kodon sözlüklerinin sırası (`UUU`, `UUC`, ... `GGG`) |
| `AMINO_ACIDS` | `"ACDEFGHIKLMNPQRSTVWY"` | Amino asit sayımlarının sırası |
| `DNA_COMPLEMENT` | `str.translate` tablosu | DNA tümleme; belirsiz kodları da tümler (`R` ↔ `Y`) |
| `RNA_COMPLEMENT` | `str.translate` tablosu | RNA tümleme |
| `DNA_TO_RNA` | `T` → `U` | Transkripsiyon |
| `RNA_TO_DNA` | `U` → `T` | `U` içeren metni DNA alfabesine çevirir; ters transkripsiyon bunu kullanır |
| `TRANSCRIPT_PRIORITY` | `("NM_", "XM_")` | `Gene.isoforms` aynı CDS'li kayıtlardan hangisini tutacağına bu sırayla karar verir |

```python
from crn.constants import CODONS, DNA_COMPLEMENT, DNA_TO_RNA, TRANSCRIPT_PRIORITY

print(len(CODONS), CODONS[0], CODONS[-1])  # -> 64 UUU GGG
print("ATGRN".translate(DNA_COMPLEMENT))  # -> TACYN
print("ATG".translate(DNA_TO_RNA))  # -> AUG
print(TRANSCRIPT_PRIORITY)  # -> ('NM_', 'XM_')
```

## 2. Genetik kodlar

### GeneticCode

Çeviride, ORF aramasında ve kodon sayımında kullanılacak genetik kodu seçen numaralandırmadır. İki üyesi vardır: `GeneticCode.STANDARD` (NCBI 1) ve `GeneticCode.VERTEBRATE_MITOCHONDRIAL` (NCBI 2). Her üyenin değeri NCBI tablo numarası, ad, 64 harflik amino asit dizisi ve başlangıç kodonlarından oluşan bir demettir. Kodonlar RNA alfabesiyle yazılır, stoplar `*` ile gösterilir. Aşağıdaki türetilmiş kümeler ilk kullanımda bir kez hesaplanıp saklanır.

```python
from crn import GeneticCode

kod = GeneticCode.VERTEBRATE_MITOCHONDRIAL
print(kod)  # -> GeneticCode.VERTEBRATE_MITOCHONDRIAL
print(kod.table["UGA"], GeneticCode.STANDARD.table["UGA"])  # -> W *
```

### GeneticCode.ncbi_id

`GeneticCode.ncbi_id  (özellik)`

NCBI genetik kod tablosu numarasını verir.

```python
print(GeneticCode.STANDARD.ncbi_id)  # -> 1
```

### GeneticCode.table

`GeneticCode.table  (özellik)`

Kodondan amino aside giden sözlüktür. Stop kodonlarının değeri `*` olur.

```python
print(GeneticCode.STANDARD.table["AUG"])  # -> M
```

### GeneticCode.starts

`GeneticCode.starts  (özellik)`

Kodun tanıdığı bütün başlangıç kodonlarını `frozenset` olarak verir. `find_orfs(alternative_starts=True)` bu kümeyi kullanır.

```python
print(sorted(GeneticCode.STANDARD.starts))  # -> ['AUG', 'CUG', 'GUG', 'UUG']
```

### GeneticCode.stops

`GeneticCode.stops  (özellik)`

Stop kodonlarının kümesidir. Mitokondriyal kodda `UGA` stop değildir, `AGA` ve `AGG` stoptur.

```python
print(sorted(GeneticCode.VERTEBRATE_MITOCHONDRIAL.stops))  # -> ['AGA', 'AGG', 'UAA', 'UAG']
```

### GeneticCode.sense_codons

`GeneticCode.sense_codons  (özellik)`

Amino asit kodlayan kodonları `CODONS` sırasıyla demet olarak verir.

```python
print(len(GeneticCode.STANDARD.sense_codons))  # -> 61
```

### GeneticCode.families

`GeneticCode.families  (özellik)`

Amino asitten onu kodlayan kodonlara giden sözlüktür. Eş anlamlı kodon ailesi tanımı bütün `crn.cub` hesaplarında buradan gelir.

```python
print(GeneticCode.STANDARD.families["L"])  # -> ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG')
```

### GeneticCode.synonymous_codons

`GeneticCode.synonymous_codons  (özellik)`

Amino asidi birden çok kodonla kodlanan kodonların kümesidir. Standart kodda Met, Trp ve stoplar bu kümenin dışında kalır, yani 59 kodon içerir.

```python
print(len(GeneticCode.STANDARD.synonymous_codons))  # -> 59
```

### GeneticCode.fourfold_codons

`GeneticCode.fourfold_codons  (özellik)`

Üçüncü konumu tümüyle eş anlamlı `xyN` kutularındaki kodonları verir. PR2 analizi bu kümeyi kullanır çünkü bu konumda baz değişimi amino asidi değiştirmez.

```python
print(len(GeneticCode.STANDARD.fourfold_codons))  # -> 32
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_codon_table` | NCBI biçimindeki 64 harflik diziyi `CODONS` sırasıyla eşleyip kodon sözlüğüne çevirir. |
| `GeneticCode.__repr__` | `<GeneticCode.STANDARD: NCBI 1, Standard>` biçiminde yazar. |
