# crn_io_features

`crn/io/crn_io_features.py` modülü, dizilere ait her türlü ek bilgiyi tutan **`Features`** sınıfını tanımlar.

Örnek betik: [tests/io/test_features.py](../tests/io/test_features.py)

## Amaç ve tasarım ilkesi

`Features`, `Sequence`, `DNA`, `RNA` ve `PRT` sınıflarına bilgi eklemek için kullanılan **ana omurgadır**.

- Bir dizinin nereden indirildiği, hangi gene ait olduğu, lokus tipi, versiyonu, strand bilgisi gibi **her türlü ek bilgi** `Features` içine yazılır.
- Dizi sınıflarına bu bilgiler için **ayrı property eklenmez**. Yeni bir bilgi türü gerektiğinde sınıfı değiştirmek yerine yeni bir anahtar kullanılır.
- İstisnalar (sınıfa doğrudan property olarak eklenmesi gereken bilgiler) tartışılarak belirlenir. Şu ana kadar kararlaştırılan root alanları: `_type` (tüm diziler), `strand` ve `maps` (DNA).

```text
Sequence / DNA / RNA / PRT
 ├── dizinin kendisi
 └── features : Features
       ├── "source"     -> ["NCBI"]
       ├── "gene"       -> ["BRCA1", "BRCA2"]
       ├── "version"    -> ["1.2"]
       └── "gene_id"    -> [672]
```

## İçe aktarma

```python
from crn.io import Features
# veya
from crn.io.crn_io_features import Features
```

## Veri modeli

- Her **anahtar** boş olmayan bir `str` olmalıdır.
- Her anahtar bir **değer listesi** tutar. Aynı anahtara birden fazla değer eklenebilir.
- Değerler herhangi bir Python nesnesi olabilir (`str`, `int`, `tuple`, `Strand` enum'u vb.). CRN dosyasına yazılabilecek tipler: `None`, `bool`, `int`, `float`, `str`, `bytes`, `Molecule`/`Strand`/`GenomicStrand`/`Maps` ve bunlardan oluşan `tuple`/`list`.
- Anahtarlar **eklenme sırasını korur**.
- Okuma metotları (`get`, `items`, `to_dict`, `copy`) listelerin **kopyasını** döndürür. Dönen liste değiştirilse bile depo etkilenmez.

## Metot özeti

| Metot / işlem               | Açıklama                                                         | Dönüş                              |
|-----------------------------|------------------------------------------------------------------|------------------------------------|
| `Features()`                | Boş bir depo oluşturur.                                          | `Features`                         |
| `add(key, value)`           | Değeri anahtarın listesinin sonuna ekler.                        | `None`                             |
| `set(key, value)`           | Anahtarın tüm değerlerini silip yerine tek değer koyar.          | `None`                             |
| `get(key, default=None)`    | Değer listesinin kopyasını döndürür.                             | `List[Any]` veya `default`         |
| `first(key, default=None)`  | Anahtarın ilk değerini döndürür.                                 | `Any`                              |
| `keys()`                    | Anahtarları sırasıyla döndürür.                                  | `List[str]`                        |
| `items()`                   | `(anahtar, değer listesi)` çiftleri üzerinde dolaşır.            | `Iterator[Tuple[str, List[Any]]]`  |
| `to_dict()`                 | Tüm özellikleri sözlük olarak döndürür.                          | `Dict[str, List[Any]]`             |
| `copy()`                    | Bağımsız (derin) bir kopya döndürür.                             | `Features`                         |
| `remove(key)`               | Anahtarı tüm değerleriyle siler.                                 | `bool`                             |
| `clear()`                   | Tüm özellikleri siler.                                           | `None`                             |
| `key in f`                  | Anahtarın var olup olmadığını kontrol eder.                      | `bool`                             |
| `len(f)`                    | Anahtar sayısı.                                                  | `int`                              |
| `for key in f`              | Anahtarlar üzerinde dolaşır.                                     | —                                  |
| `repr(f)`                   | `Keys:source,gene` biçiminde gösterim.                           | `str`                              |

## Ayrıntılı kullanım

### Oluşturma

```python
f = Features()
len(f)      # 0
f           # Keys:
```

### `add` — çok değerli bilgi ekleme

Aynı anahtara birden fazla değer eklemek için kullanılır. Anahtar yoksa oluşturulur.

```python
f.add("gene", "BRCA1")
f.add("gene", "BRCA2")
f.get("gene")    # ['BRCA1', 'BRCA2']
```

### `set` — tek değerli bilgi yazma

Versiyon, strand, kaynak gibi **tek değer** taşıması gereken bilgiler için kullanılır. Önceki değerlerin hepsi silinir.

```python
f.add("version", "1.0")
f.add("version", "1.1")
f.set("version", "2.0")
f.get("version")   # ['2.0']
```

> **Kural:** Tek değerli bilgilerde `add` yerine `set` kullanın. `add` kullanılırsa aynı anahtarda farkında olmadan birden fazla değer birikebilir.

### `get` — değer listesini okuma

Anahtar varsa değer listesi, yoksa `default` (varsayılan `None`) döner.

```python
f.get("gene")      # ['BRCA1', 'BRCA2']
f.get("yok")       # None
f.get("yok", [])   # []

for v in f.get("gene", []):
    print(v)
```

Dönen liste bir kopyadır:

```python
genes = f.get("gene")
genes.append("TP53")
f.get("gene")      # ['BRCA1', 'BRCA2']  -> depo değişmedi
```

Depoya değer eklemek için her zaman `add` veya `set` kullanılmalıdır.

### `first` — tek değeri okuma

Anahtarın ilk değerini döndürür. Anahtar yoksa `default` döner. `set` ile yazılan bilgileri okumanın en kısa yoludur.

```python
f.set("source", "NCBI")
f.first("source")              # 'NCBI'
f.first("organism")            # None
f.first("organism", "unknown") # 'unknown'
```

### `keys`, `items`, döngü

```python
f.keys()                 # ['gene', 'version', 'source']

for key in f:
    print(key)

for key, values in f.items():
    print(key, values)
```

### Varlık kontrolü ve sayım

```python
"gene" in f      # True
"locus" in f     # False
len(f)           # 3
```

### `to_dict` — dışa aktarma

Tüm bilgileri kopya bir sözlük olarak verir. JSON'a yazma, loglama veya başka bir formata dönüştürme için uygundur.

```python
import json

f.to_dict()
# {'gene': ['BRCA1', 'BRCA2'], 'version': ['2.0'], 'source': ['NCBI']}

json.dumps(f.to_dict())
```

> Değerler arasında `Strand` gibi enum'lar varsa `json.dumps` doğrudan çalışmaz; önce `f"{type(v).__name__}.{v.name}"` gibi bir dönüşüm yapılmalıdır. Enum'ları kayıpsız saklamak için `to_crn` kullanın.

### `copy` — bağımsız kopya

Dizi dilimlenirken, `complementary()` ya da `DNA.from_sequence()` çağrılırken features bu metotla kopyalanır. Kopyada yapılan değişiklik orijinali etkilemez.

```python
k = f.copy()
k.add("gene", "TP53")
f.get("gene")      # ['BRCA1', 'BRCA2']
```

### `remove` ve `clear` — silme

```python
f.remove("version")   # True  -> silindi
f.remove("version")   # False -> zaten yoktu

f.clear()
len(f)                # 0
```

## Hata durumları

`add` ve `set` anahtarı doğrular:

```python
f.add(1, "x")     # TypeError: Anahtar str olmalı, int verildi
f.set("", "x")    # ValueError: Anahtar boş olamaz
```

`get`, `first`, `remove` ve `in` geçersiz anahtarlarda hata vermez; sadece "bulunamadı" sonucunu döndürür (`None`, `default`, `False`).

## `crn_io_type` ile birlikte kullanım

Tip ve yön bilgileri `str` yerine [`crn_io_type`](crn_io_type.md) içindeki enum'larla saklanmalıdır:

```python
from crn.io import Features, Molecule, Strand, GenomicStrand

f = Features()
f.set("molecule", Molecule.RNA)
f.set("strand", Strand.FivePrime)
f.set("genomic_strand", GenomicStrand.Minus)   # yalnızca RNA

f.first("strand") is Strand.FivePrime      # True
f.first("genomic_strand").value[0]          # '-'
```

## Dizi sınıflarında kullanım

Her `Sequence` (ve `DNA`) nesnesinin bir `features` alanı vardır:

```python
from crn.io import DNA

d = DNA("ATGGATTTATCTGCTCTTCGC", name="BRCA1_part")
d.features.set("source", "NCBI")
d.features.set("accession", "NM_007294")
d.features.set("version", "4")
d.features.add("gene", "BRCA1")

d.features.first("source")   # 'NCBI'
d[0:9].features.get("gene")  # ['BRCA1']  -> dilim kendi kopyasını taşır
```

Yeni bir bilgi gerektiğinde (örneğin lokus tipi) sınıfa dokunulmaz:

```python
d.features.set("locus_type", "protein-coding")
```

`features` alanına yalnızca bir `Features` nesnesi atanabilir; başka bir tip `TypeError` verir.

`from_fasta_file` her diziye `fasta_header` ve `source_file` anahtarlarını kendisi ekler.

## Önerilen anahtar adları

Aynı bilginin farklı yerlerde farklı adlarla (`"Source"`, `"src"`, `"source"`) yazılmasını önlemek için ortak anahtar adları kullanılmalıdır. Başlangıç önerisi:

| Anahtar        | Değer tipi        | Tek / çok | Açıklama                                   |
|----------------|-------------------|-----------|--------------------------------------------|
| `source`       | `str`             | tek       | Dizinin indirildiği yer (NCBI, Ensembl …)  |
| `accession`    | `str`             | tek       | Kaynaktaki kayıt numarası                  |
| `version`      | `str`             | tek       | Kayıt versiyonu                            |
| `strand`       | `Strand`          | tek       | Yazılış yönü (DNA'da root alanıdır)        |
| `genomic_strand` | `GenomicStrand` | tek       | Genomik zincir (yalnızca RNA)              |
| `fasta_header` | `str`             | tek       | FASTA başlığının tamamı (otomatik)         |
| `source_file`  | `str`             | tek       | Okunan dosyanın yolu (otomatik)            |
| `gene`         | `str`             | çok       | İlgili gen adları                          |
| `locus_type`   | `str`             | tek       | Lokus tipi                                 |
| `gene_type`    | `str`             | tek       | Gen tipi                                   |

Bu liste ihtiyaç oldukça genişletilecektir. Anahtar adları küçük harf ve `snake_case` olmalıdır.
