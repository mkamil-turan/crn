# crn_io_type

`crn/io/crn_io_type.py` modülü dizilerle çalışırken kullanılan üç enum'u tanımlar:

| Enum            | Ne anlatır                                          | Nerede tutulur                     |
|-----------------|-----------------------------------------------------|------------------------------------|
| `Molecule`      | Dizinin molekül tipi (DNA, RNA, protein, genel)     | Dizi nesnesinin `_type` alanı      |
| `Strand`        | Dizinin yazılış yönü (5'→3' / 3'→5')                | DNA'nın `strand` alanı, features   |
| `GenomicStrand` | RNA'nın transkript edildiği genomik zincir (+/−)    | Yalnızca RNA features'ı            |

Üyelerin değerleri birer `tuple`'dır. Karşılaştırmada `is` kullanın (`m is Molecule.DNA`).

Örnek betik: [tests/io/test_type.py](../tests/io/test_type.py)

## İçe aktarma

```python
from crn.io import Molecule, Strand, GenomicStrand
```

## `Molecule`

| Üye            | `value[0]`   | `value[1]` |
|----------------|--------------|------------|
| `Molecule.DNA` | `"Dna"`      | `"DNA"`    |
| `Molecule.RNA` | `"Rna"`      | `"RNA"`    |
| `Molecule.PRT` | `"Protein"`  | `"PRT"`    |
| `Molecule.SEQ` | `"Sequence"` | `"SEQ"`    |

Her dizi sınıfının bir `_type` alanı vardır (`Sequence` → `SEQ`, `DNA` → `DNA`). Bu alan features'ta değil, nesnenin kendisinde (root) tutulur ve CRN dosyasına da root alanı olarak yazılır.

```python
from crn.io import DNA, Molecule

DNA("ATG")._type is Molecule.DNA   # True
Molecule["PRT"]                    # Molecule.PRT
```

## `Strand`

Dizinin hangi uçtan hangi uca yazıldığını gösterir.

| Üye                 | `value[0]` (başlangıç) | `value[1]` (bitiş) | `value[2]`      |
|---------------------|------------------------|--------------------|-----------------|
| `Strand.FivePrime`  | `"5'"`                 | `"3'"`             | `"Five prime"`  |
| `Strand.ThreePrime` | `"3'"`                 | `"5'"`             | `"Three prime"` |

```python
start, end, label = Strand.FivePrime.value
f"{start}-ATGGCTAAG-{end}"      # "5'-ATGGCTAAG-3'"
```

`Strand` bir `+`/`−` işareti taşımaz: 5'→3' yazılmış bir dizi genomun (−) zincirinden de gelmiş olabilir. Genomik zincir için `GenomicStrand` kullanılır.

## `GenomicStrand`

RNA'nın referans genomda hangi zincirden transkript edildiğini gösterir. GFF/GTF/BED dosyalarının strand sütunuyla aynı notasyonu kullanır.

> **Yalnızca RNA özelliğidir; DNA dizileri için kullanılmaz.**

| Üye                        | `value[0]` | `value[1]`     | Anlamı                              |
|----------------------------|------------|----------------|-------------------------------------|
| `GenomicStrand.Plus`       | `"+"`      | `"Plus"`       | Referans genomda yazılan zincir     |
| `GenomicStrand.Minus`      | `"-"`      | `"Minus"`      | Referans zincirin tümleyeni         |
| `GenomicStrand.Unstranded` | `"."`      | `"Unstranded"` | Zincir bilgisi anlamsız             |
| `GenomicStrand.Unknown`    | `"?"`      | `"Unknown"`    | Zincir anlamlı ama bilinmiyor       |

Örnek: MEFV geni kromozom 16'nın (−) zincirindedir. mRNA'sı 5'→3' yazılır:

```python
mrna.features.set("strand", Strand.FivePrime)
mrna.features.set("genomic_strand", GenomicStrand.Minus)
```

GFF satırından okumak için:

```python
def genomic_strand_from_gff(sign: str) -> GenomicStrand:
    for g in GenomicStrand:
        if g.value[0] == sign:
            return g
    raise ValueError(f"Geçersiz GFF strand değeri: {sign}")
```

## `Maps` (crn_io_const)

Tümleyici çeviri tablolarını tutar. Varsayılan tablo IUPAC kodlarını da çevirir:

| Kod | A | C | G | T | R | Y | K | M | S | W | B | V | D | H | N | - |
|-----|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Tümleyen | T | G | C | A | Y | R | M | K | S | W | V | B | H | D | N | - |

```python
from crn.io import Maps
"ATGCRY".translate(Maps.DEFAULT_COMPLEMENT_MAP_FOR_DNA.value)   # 'TACGYR'
```

## Notlar

- Değerler `tuple` olduğundan `Strand("5'")` gibi tek alanla arama çalışmaz; yukarıdaki gibi döngüyle arayın.
- CRN dosyasına `Molecule`, `Strand`, `GenomicStrand` ve `Maps` değerleri enum adı + üye adı olarak yazılır ve aynen geri okunur.
