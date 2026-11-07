# `crn.fileio`: FASTA ve CRN dosyaları

Kaynak: [`crn/fileio.py`](../crn/fileio.py)

İki biçim vardır. FASTA başka programlarla alışveriş içindir ve kayıplıdır: tip, yön ve features yazılmaz. CRN paketin kendi ikili biçimidir ve nesneyi eksiksiz saklar. `pickle` kullanılmaz, bu yüzden güvenilmeyen bir CRN dosyasını açmak kod çalıştırmaz; her dosyada CRC32 denetimi vardır.

<!-- toc -->
- [1. FASTA](#1-fasta): [`read_fasta`](#read_fasta) · [`write_fasta`](#write_fasta)
- [2. CRN yazma ve okuma](#2-crn-yazma-ve-okuma): [`write_crn`](#write_crn) · [`read_crn`](#read_crn) · [`write_gene_crn`](#write_gene_crn) · [`read_gene_crn`](#read_gene_crn) · [`write_crn_list`](#write_crn_list) · [`read_crn_list`](#read_crn_list)
- [3. CRN dosya biçimi](#3-crn-dosya-biçimi)
- [4. İç yardımcılar](#4-iç-yardımcılar)
<!-- /toc -->

## 1. FASTA

### read_fasta

`read_fasta(file_name, cls=Sequence)`

`file_name` FASTA dosyasındaki dizileri `cls` tipinde liste olarak okur. `.gz` dosyalar gzip ile açılır. Dizinin adı başlığın ilk kelimesidir. Boş satırlar ve `;` yorumları atlanır, çok satırlı diziler birleştirilir. features'a `fasta_header` (başlığın tamamı) ve `source_file` yazılır. Dosya yoksa `FileNotFoundError`; başlıksız veri, boş başlık, dizisiz başlık, geçersiz karakter ya da hiç kayıt olmaması durumunda dosya ve satır numarası içeren `ValueError` verilir.

```python
from pathlib import Path
from crn import DNA, RNA, Gene, read_fasta, write_fasta, read_crn_list, write_crn_list

Path("g.fa").write_text(">a insan geni\nATGC\nATGC\n; yorum\n>b\nGGCC\n")
kayitlar = read_fasta("g.fa", cls=DNA)
print([(k.name, k.seq) for k in kayitlar])  # -> [('a', 'ATGCATGC'), ('b', 'GGCC')]
print(kayitlar[0].features.first("fasta_header"))  # -> a insan geni
```

### write_fasta

`write_fasta(sequences, file_name)`

`sequences` içindeki dizileri `file_name` dosyasına yazar, varsa üzerine yazar. Başlık `>AD this sequence contains N feature(s) but not listed` biçimindedir: N, `fasta_header` ve `source_file` dışındaki anahtar sayısıdır ve okura features'ın kaybolduğunu söyler. Dizi `FASTA_LINE_WIDTH` (60) karakterlik satırlara bölünür. Adı boş ya da boşluklu dizi veya boş dizi `ValueError` verir ve bu durumda hiçbir şey yazılmaz.

```python
write_fasta([DNA("ATGC", name="a")], "cikti.fa")
print(Path("cikti.fa").read_text().splitlines())  # -> ['>a this sequence contains 0 features but not listed', 'ATGC']
```

## 2. CRN yazma ve okuma

`Sequence.to_crn`, `Sequence.from_crn`, `Gene.to_crn` ve `Gene.from_crn` aşağıdaki fonksiyonları çağırır. `crn` paketinin kökünden yalnızca liste fonksiyonları dışa verilir; diğerleri `crn.fileio` üzerinden alınır.

### write_crn

`write_crn(seq, file_name)`

`seq` dizisini `file_name` dosyasına yazar. Features'ta desteklenmeyen tipte bir değer varsa `TypeError` verir ve dosyaya dokunmaz.

```python
from crn.fileio import write_crn, read_crn, write_gene_crn, read_gene_crn

d = DNA("ATGC", name="d")
d.features.add("cds", (1, 4))
write_crn(d, "d.crn")
```

### read_crn

`read_crn(file_name, cls=Sequence)`

`write_crn` ile yazılmış `file_name` dosyasını okur. `cls` varsayılan `Sequence` ise dosyadaki tip seçilir; başka sınıf verilirse dosyadaki tip onunla uyuşmalıdır. Dosya yoksa `FileNotFoundError`; CRN değilse, bozuksa, sürümü desteklenmiyorsa ya da tip uyuşmuyorsa `ValueError` verir.

```python
geri = read_crn("d.crn", cls=DNA)
print(geri.features.first("cds"))  # -> (1, 4)
```

### write_gene_crn

`write_gene_crn(gene, file_name)`

`gene` nesnesini (ad, features, genomik DNA, transkriptler, proteinler) `file_name` dosyasına yazar.

```python
write_gene_crn(Gene("g", dna=DNA("ATG")), "g.crn")
```

### read_gene_crn

`read_gene_crn(file_name)`

`write_gene_crn` ile yazılmış `file_name` dosyasını okur. Dosya gen dosyası değilse ya da bozuksa `ValueError` verir.

```python
print(read_gene_crn("g.crn"))  # -> g: DNA <3>, 0 transkript, 0 protein
```

### write_crn_list

`write_crn_list(items, file_name)`

`items` içindeki `Gene` ve dizi nesnelerini sırayla tek `file_name` dosyasına yazar. Liste karışık olabilir, aynı addan birden çok kayıt olabilir, dosya varsa üzerine yazılır. Gen gruplarını tek dosyada tutmanın yolu budur: `DifferentialCodonUsage.from_crn` bu dosyaları okur. Eleman `Gene` ya da `Sequence` değilse `TypeError` verir ve dosyaya dokunmaz.

```python
write_crn_list([Gene("g", dna=DNA("ATG")), RNA("AUG", name="r")], "liste.crn")
```

### read_crn_list

`read_crn_list(file_name)`

`write_crn_list` ile yazılmış `file_name` dosyasını yazıldığı sırayla okur. Her eleman dosyadaki tipine uygun sınıftadır.

```python
print([type(x).__name__ for x in read_crn_list("liste.crn")])  # -> ['Gene', 'RNA']
```

## 3. CRN dosya biçimi

Bütün sayılar little-endian'dır. Dosya `CRN_MAGIC` (`b"CRN\x00"`), `uint16` `CRN_VERSION` (şu an 3), bir kayıt ve baştan itibaren bütün baytların `uint32` CRC32'sinden oluşur.

```text
dizi kaydı   type str (SEQ, DNA, RNA, PRT), name str, seq str,
             n_root uint32 + n_root × (key str, etiketli değer), features
gen kaydı    type str "GENE", name str, features,
             has_dna uint8 (1 ise DNA'nın dizi kaydı),
             n_transcripts uint32 + n × (id str, dizi kaydı),
             n_proteins uint32 + n × (id str, dizi kaydı)
liste kaydı  type str "LIST", n_items uint32 + n × (dizi ya da gen kaydı)
features     uint32 anahtar sayısı + her anahtar için
             key str, uint32 değer sayısı, etiketli değerler
str          uint32 bayt uzunluğu + UTF-8 baytları
```

Etiketli değer bir baytlık etiket ve içerikten oluşur.

| Etiket | Tip | İçerik |
|---|---|---|
| `N`, `T`, `F` | `None`, `True`, `False` | Yok |
| `i` | int | `uint8` uzunluk + işaretli little-endian baytlar |
| `f` | float | `float64` |
| `s`, `b` | str, bytes | str; `uint32` uzunluk + baytlar |
| `e` | enum | Sınıf adı ve üye adı; yalnızca `CRN_ENUMS` içindekiler (`Molecule`, `Strand`, `GenomicStrand`, `GeneticCode`) |
| `t`, `l` | tuple, list | `uint32` eleman sayısı + etiketli değerler |

Features'a yalnızca bu tipler yazılabilir. `dict`, `set` ve `numpy.int64` gibi tipler `TypeError` verir, numpy tam sayıları önce `int()` ile çevrilmelidir; `numpy.float64` `float` alt sınıfı olduğu için yazılabilir. Sürüm 3 dosyalarında kök alanlarda bulunan eski `Maps` enum'u `_LEGACY_ENUMS` sayesinde okunurken yok sayılır.

## 4. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_FASTA_KEYS` | `write_fasta`'nın başlıktaki sayıya katmadığı anahtarlar: `fasta_header`, `source_file` |
| `_pack_str` | Metni uzunluk önekli UTF-8 baytlara çevirir |
| `_pack_value` | Bir değeri etiketli baytlara çevirir; desteklenmeyen tipte `TypeError` |
| `_pack_features` | `Features` nesnesini paketler |
| `_pack_sequence` | Dizi kaydını paketler |
| `_pack_gene` | Gen kaydını paketler |
| `_write` | Magic, sürüm ve gövdeyi birleştirip CRC32 ekler ve dosyaya yazar |
| `_Reader` | Dosyayı açar, magic, CRC32 ve sürümü denetler, baytları sırayla okur |
| `_Reader.__init__` | Bütün dosyayı belleğe alır ve üç denetimi yapar |
| `_Reader.take` | `n` bayt okur; dosya kısaysa `ValueError` |
| `_Reader.uint` | `struct` biçimiyle bir sayı okur |
| `_Reader.text` | Uzunluk önekli metin okur |
| `_Reader.finish` | Dosyanın sonunda fazladan veri olmadığını denetler |
| `_Reader.value` | Etiketli değer okur; bilinmeyen etikette `ValueError` |
| `_Reader.features` | `Features` okur |
| `_Reader.sequence` | Dizi kaydı okur ve hedef sınıfı seçer; dosya gen ya da liste dosyasıysa doğru okuyucuyu bildirir |
| `_Reader.gene` | Tipi okunmuş gen kaydını okur |
