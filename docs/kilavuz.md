# crn Kullanım Kılavuzu

`crn`, DNA, RNA ve protein dizilerini nesne olarak modelleyen, aralarındaki
dönüşümleri (transkripsiyon, ters transkripsiyon, translasyon) yapan, ORF
arayan, dizileri FASTA ve kendi ikili biçimi olan CRN ile okuyup yazan ve
GENCODE ile NCBI'dan gen, transkript ve protein çeken bir Python paketidir.

Bu kılavuz veri modelinin (`crn`) ve veri kaynaklarının (`crn.fetch`) genel
(public) sınıflarını ve metotlarını tek tek anlatır. Kodon kullanımı (`crn.cub`)
ve protein özellikleri (`ProteinProperties`) için [README](../README.md)'ye bakın. Her metot için amacı, parametreleri, dönüş değeri, verdiği hatalar ve
çalıştırılmış örnekleri verilir. Örneklerdeki çıktılar gerçek çıktılardır.
Paketin sınırlılıkları en sondaki bölümde toplanmıştır.

---

## İçindekiler

1. [Kurulum ve içe aktarma](#1-kurulum-ve-içe-aktarma)
2. [Temel tasarım ilkeleri](#2-temel-tasarım-ilkeleri)
3. [Sabitler ve enum'lar](#3-sabitler-ve-enumlar) — `Molecule`, `Strand`, `GenomicStrand`, `GeneticCode`, tümleme tabloları, alfabeler
4. [`Features`](#4-features) — ek bilgilerin deposu
5. [`Sequence`](#5-sequence) — tüm dizilerin temel sınıfı
6. [`DNA`](#6-dna)
7. [`RNA`](#7-rna)
8. [`PRT`](#8-prt) — protein
9. [`Gene`](#9-gene) — DNA, RNA ve proteinler için kapsayıcı
10. [CRN dosya biçimi](#10-crn-dosya-biçimi)
11. [`NcbiSession`](#11-ncbisession) — NCBI HTTP oturumu
12. [`Gencode`](#12-gencode) — GENCODE'dan veri çekme
13. [`Ncbi`](#13-ncbi) — NCBI'dan veri çekme
14. [Uçtan uca senaryolar](#14-uçtan-uca-senaryolar)
15. [Hata mesajları](#15-hata-mesajları)
16. [Sınırlılıklar](#16-sınırlılıklar)

---

## 1. Kurulum ve içe aktarma

### Gereksinimler

- Python **3.12 veya üstü**.
- Bağımlılıklar `pip` tarafından otomatik kurulur: `tqdm` (ilerleme çubukları),
  `requests` (NCBI istekleri) ve `matplotlib` (yalnızca CUB grafikleri çizilirken
  yüklenir).

### Kurulum

Proje klasöründe:

```bash
python -m venv .venv
```

```bash
source .venv/bin/activate.fish
```

```bash
pip install -e .
```

bash/zsh kullanıyorsanız ikinci komut `source .venv/bin/activate` olur.
`-e` (düzenlenebilir kurulum) sayesinde kaynak kodda yaptığınız değişiklikler
yeniden kurulum gerektirmez.

### Paket yapısı

```text
crn/
├── constants.py      alfabeler, kodonlar, tümleme tabloları, GeneticCode
├── tables.py         kütleler, pKa, hidropati ölçekleri, tAI cezaları
├── features.py       Features
├── seq.py            Molecule, Strand, GenomicStrand, Sequence, DNA, RNA, PRT
├── gene.py           Gene
├── fileio.py         FASTA ve CRN biçimleri
├── protein.py        ProteinProperties
├── fetch/            veri kaynakları
│   ├── session.py    NcbiSession
│   ├── gencode.py    Gencode
│   └── ncbi.py       Ncbi
└── cub/              kodon kullanım yanlılığı (README'ye bakın)
sys/                  indirilen veriler (Gencode ve Ncbi yazar; CRN_SYS_DIR ile değişir)
```

### İçe aktarma

```python
# Diziler, kapsayıcılar ve sabitler
from crn import (Molecule, Strand, GenomicStrand, GeneticCode,
                 Features, Sequence, DNA, RNA, PRT, Gene)

# Veri kaynakları
from crn.fetch import Gencode, Ncbi, NcbiSession

# Alfabeler ve tümleme tabloları (gerekirse)
from crn.constants import DNA_ALPHABET, RNA_ALPHABET, PRT_ALPHABET, DNA_COMPLEMENT
```

---

## 2. Temel tasarım ilkeleri

Paketi doğru kullanmak için bu bölümdeki kuralları bilmek gerekir.

### 2.1. Ek bilgiler `features` içinde tutulur

Bir dizinin nereden geldiği, hangi gene ait olduğu, sürümü, CDS koordinatları
gibi **her türlü ek bilgi** nesnenin `features` alanına (bir `Features`
nesnesi) yazılır. Sınıflara bu bilgiler için yeni özellik (property) eklenmez.

```python
rna.features.set("gene_name", "MEFV")
rna.features.first("gene_name")     # 'MEFV'
```

### 2.2. Root alanları (istisnalar)

Bazı bilgiler features'ta değil, nesnenin kendisinde tutulur. Bunlara
**root alanı** denir. Root alanları kopyalarda korunur ve CRN dosyasına ayrı
bir bölüm olarak yazılır.

| Sınıf | Root alanları |
|---|---|
| Tüm diziler | `name`, `seq` (ve sınıf alanı `molecule`) |
| `DNA` | `strand` |
| `RNA` | `strand`, `genomic_strand` |
| `PRT` | — |
| `Gene` | `name`, `dna`, transkriptler, proteinler |

Molekül tipi sınıfın kendisindedir: `DNA("ATG").molecule is Molecule.DNA`.
Bir sınıfın root alanları `root_types` sözlüğünde, bir nesneninkiler
`root_fields()` ile okunur.

### 2.3. `strand` = yazılış yönü

`Strand`, dizinin **hangi uçtan hangi uca yazıldığını** gösterir:

- `Strand.FivePrime`: dizi soldan sağa 5'→3' yazılmıştır (normal durum).
- `Strand.ThreePrime`: dizi soldan sağa 3'→5' yazılmıştır.

Tümleme (`complement`) diziyi tersine çevirmez; tümleyici zinciri ilk zincirin
altına **hizalı** yazar. Bu yüzden tümleyen dizinin yönü 3'→5' olur:

```text
5'-ATGGCC-3'      DNA("ATGGCC")                 strand = FivePrime
3'-TACCGG-5'      DNA("ATGGCC").complementary() strand = ThreePrime
```

Yönü 5'→3' olarak okumak gereken işlemler (`translation`, `find_orfs`) 3'→5'
yazılmış diziyi önce kendileri ters çevirir.

### 2.4. `genomic_strand` yalnızca RNA'dadır

`GenomicStrand` (+, −, ., ?), bir RNA'nın referans genomun hangi zincirinden
transkript edildiğini gösterir. **Yalnızca RNA'nın** root alanıdır; DNA'da
yoktur. Bir DNA'nın kromozomdaki zinciri gerekiyorsa features'a yazılır
(`chromosome_strand`).

### 2.5. Diziler büyük harfe çevrilir ve doğrulanır

Her dizi oluşturulurken büyük harfe çevrilir. `DNA`, `RNA` ve `PRT` yalnızca
kendi alfabelerindeki karakterleri kabul eder (bkz. 3.6). `Sequence` hiçbir
karakter kontrolü yapmaz.

### 2.6. Kopyalar bağımsızdır

Dilimleme (`s[a:b]`), `complementary()`, `from_sequence()`, `transcription()`,
`reverse_transcription()` ve `translation()` yeni bir nesne döndürür ve
features'ın **bağımsız bir kopyasını** verir. Yeni nesnenin features'ını
değiştirmek eskisini etkilemez.

### 2.7. Dönüş tipleri ve editör desteği

`from_crn`, `from_fasta_file`, `from_sequence`, dilimleme ve `complementary`
gibi metotlar `Self` döndürür: metot hangi sınıftan çağrıldıysa o sınıfın
nesnesi döner. Böylece editörünüz `DNA.from_crn(...)` sonucunda DNA'nın
metotlarını önerir.

### 2.8. Hatalar sessizce yutulmaz

Bulunamayan kayıtlar için `KeyError`, geçersiz değerler için `ValueError`,
yanlış tipler için `TypeError` verilir. Hiçbir metot hata durumunda sessizce
`None` döndürmez (yalnızca `Features.get/first` gibi okuma metotları,
belgelenmiş varsayılan değeri döndürür).

---

## 3. Sabitler ve enum'lar

`Molecule` dışındaki enum üyelerinin değeri bir `tuple`'dır. Karşılaştırmada `is` kullanın
(`x is Strand.FivePrime`). Adından üye almak için `Molecule["DNA"]` yazılır.

### 3.1. `Molecule` — molekül tipi

**Ne için:** Her dizi sınıfının `molecule` alanında durur; dizinin DNA, RNA,
protein ya da genel dizi olduğunu gösterir. CRN dosyası okunurken hangi
sınıfın oluşturulacağı buna göre seçilir.

| Üye | `value` | Kullanan sınıf |
|---|---|---|
| `Molecule.DNA` | `"DNA"` | `DNA` |
| `Molecule.RNA` | `"RNA"` | `RNA` |
| `Molecule.PRT` | `"PRT"` | `PRT` |
| `Molecule.SEQ` | `"SEQ"` | `Sequence` |

```python
>>> Molecule.DNA, Molecule.DNA.value, Molecule["PRT"]
(Molecule.DNA, 'DNA', Molecule.PRT)
>>> DNA("ATG").molecule is Molecule.DNA
True
```

### 3.2. `Strand` — yazılış yönü

**Ne için:** `DNA` ve `RNA`'nın `strand` alanı. Anlamı için bkz. 2.3.

| Üye | `value` = (başlangıç ucu, bitiş ucu, açıklama) |
|---|---|
| `Strand.FivePrime` | `("5'", "3'", "Five prime")` |
| `Strand.ThreePrime` | `("3'", "5'", "Three prime")` |

`Strand.FivePrime.flipped` ters yönü (`Strand.ThreePrime`) verir.

```python
>>> start, end, label = Strand.FivePrime.value
>>> f"{start}-ATGGCC-{end}"
"5'-ATGGCC-3'"
```

### 3.3. `GenomicStrand` — genomik zincir (yalnızca RNA)

**Ne için:** `RNA`'nın `genomic_strand` alanı. GFF/GTF/BED dosyalarının
strand sütunuyla aynı notasyonu kullanır.

| Üye | `value` | Anlamı |
|---|---|---|
| `GenomicStrand.Plus` | `("+", "Plus")` | Referans genomda yazılan zincir |
| `GenomicStrand.Minus` | `("-", "Minus")` | Referans zincirin tümleyeni |
| `GenomicStrand.Unstranded` | `(".", "Unstranded")` | Zincir bilgisi anlamsız |
| `GenomicStrand.Unknown` | `("?", "Unknown")` | Zincir anlamlı ama bilinmiyor (varsayılan) |

Örnek: MEFV geni kromozom 16'nın (−) zincirindedir; mRNA'sı yine 5'→3' yazılır:
`RNA(..., strand=Strand.FivePrime, genomic_strand=GenomicStrand.Minus)`.

### 3.4. Tümleme tabloları

**Ne için:** Bazların birbirine dönüştürülmesi. `crn.constants` içindeki bu
sabitler Python'un `str.translate` ile kullanılan tablolardır. Dizi sınıfları
bunları kendisi seçer; elle kullanmak gerekmez.

| Sabit | Dönüşüm | Kullanıldığı yer |
|---|---|---|
| `DNA_COMPLEMENT` | A↔T, G↔C ve IUPAC: R↔Y, K↔M, B↔V, D↔H; S, W, N, `-` aynı kalır | `DNA` tümleme |
| `RNA_COMPLEMENT` | A↔U, G↔C ve aynı IUPAC kuralları | `RNA` tümleme |
| `DNA_TO_RNA` | Yalnızca T→U | `DNA.transcription(template=False)` |
| `RNA_TO_DNA` | Yalnızca U→T | RNA'dan kodlayan zincire |

Kalıp zincirden transkripsiyon `DNA_COMPLEMENT` + `DNA_TO_RNA`, ters
transkripsiyon `RNA_COMPLEMENT` + `RNA_TO_DNA` ile yapılır.

```python
>>> from crn.constants import DNA_COMPLEMENT, DNA_TO_RNA
>>> "ATGCRYN".translate(DNA_COMPLEMENT)
'TACGYRN'
>>> "ATGC".translate(DNA_COMPLEMENT).translate(DNA_TO_RNA)
'UACG'
```

Sürüm 3 CRN dosyalarında bulunan eski `maps` alanı okunurken yok sayılır.

### 3.5. `GeneticCode` — genetik kodlar

**Ne için:** `RNA.translation` ve `DNA.find_orfs` için kodon tablosu, başlangıç
ve stop kodonları. Kodonlar RNA alfabesiyle (U) yazılır; stop `*` ile gösterilir.

| Üye | NCBI tablosu | Başlangıç kodonları | Stop kodonları | Standarttan farkı |
|---|---|---|---|---|
| `GeneticCode.STANDARD` | 1 | AUG, UUG, CUG, GUG | UAA, UAG, UGA | — |
| `GeneticCode.VERTEBRATE_MITOCHONDRIAL` | 2 | AUG, AUU, AUC, AUA, GUG | UAA, UAG, AGA, AGG | UGA=W, AUA=M, AGA/AGG=stop |

Özellikler:

| Özellik | Döndürdüğü |
|---|---|
| `table` | `{kodon: amino asit}` sözlüğü (64 kodon) |
| `starts` | Başlangıç kodonlarının kümesi (`frozenset`) |
| `stops` | Stop kodonlarının kümesi (`frozenset`) |
| `sense_codons` | Amino asit kodlayan kodonlar (stop hariç; standartta 61) |
| `families` | `{amino asit: kodonlar}` eş anlamlı kodon aileleri |
| `synonymous_codons` | Amino asidi birden çok kodonla kodlanan kodonlar (standartta 59) |
| `fourfold_codons` | 3. konumu tümüyle eş anlamlı kodonlar (standartta 32) |

Türetilen kümeler ilk kullanımda bir kez hesaplanır.

```python
>>> GeneticCode.STANDARD.table["AUG"]
'M'
>>> sorted(GeneticCode.STANDARD.starts)
['AUG', 'CUG', 'GUG', 'UUG']
>>> sorted(GeneticCode.VERTEBRATE_MITOCHONDRIAL.stops)
['AGA', 'AGG', 'UAA', 'UAG']
>>> GeneticCode.VERTEBRATE_MITOCHONDRIAL.table["UGA"]
'W'
```

### 3.6. Alfabeler

`crn.constants` içinde `frozenset` olarak tanımlıdır:

| Sabit | Karakterler |
|---|---|
| `DNA_ALPHABET` | `A C G T` + IUPAC `R Y K M S W B D H V N` + boşluk `-` |
| `RNA_ALPHABET` | `A C G U` + IUPAC `R Y K M S W B D H V N` + `-` |
| `PRT_ALPHABET` | 20 standart amino asit + `B Z J X` (belirsiz) + `U` (selenosistein) + `O` (pirolizin) + `*` (stop) + `-` |

IUPAC kodları: R=A/G, Y=C/T, K=G/T, M=A/C, S=G/C, W=A/T, B=C/G/T, D=A/G/T,
H=A/C/T, V=A/C/G, N=herhangi.

---

## 4. `Features`

**Ne için:** Bir nesneye ait tüm ek bilgileri anahtar → değer listesi olarak
tutar. Her dizide (`seq.features`) ve her `Gene`'de (`gene.features`) bulunur.

**Veri modeli:**

- Anahtarlar boş olmayan `str` olmalıdır.
- Her anahtar bir **liste** tutar; aynı anahtara birden çok değer eklenebilir.
- Anahtarlar eklenme sırasını korur.
- Değerler herhangi bir Python nesnesi olabilir; ancak CRN dosyasına yazılacaksa
  izin verilen tiplerden olmalıdır (bkz. 10).
- Okuma metotları (`get`, `items`, `to_dict`) listelerin **kopyasını** döndürür.

### `Features()`

Boş bir depo oluşturur. Normalde kendiniz oluşturmazsınız; her dizi ve gen
kendi deposuyla gelir.

### `add(key, value) -> None`

**Amaç:** Anahtarın listesinin sonuna bir değer ekler; anahtar yoksa oluşturur.
Çok değerli bilgiler (diğer adlar, varyantlar, bölgeler) için kullanılır.

**Hatalar:** Anahtar `str` değilse `TypeError`, boşsa `ValueError`.

```python
>>> f = Features()
>>> f.add("alias", "TRIM20")
>>> f.add("alias", "FMF")
>>> f.get("alias")
['TRIM20', 'FMF']
```

### `set(key, value) -> None`

**Amaç:** Anahtarın tüm değerlerini silip yerine tek bir değer koyar. Tek
değerli bilgiler (sürüm, kaynak, gene_id) için kullanılır.

**Hatalar:** `add` ile aynı.

```python
>>> f.set("gene_id", 4210)
>>> f.set("gene_id", 4210)     # tekrar set etmek değer biriktirmez
>>> f.get("gene_id")
[4210]
```

> **Kural:** Tek değerli bilgilerde `add` yerine `set` kullanın; `add` ile
> aynı anahtara yazılan değerler birikir.

### `get(key, default=None) -> list | default`

**Amaç:** Anahtarın değer listesinin kopyasını döndürür; anahtar yoksa
`default`. Dönen listeyi değiştirmek depoyu etkilemez.

```python
>>> f.get("alias")
['TRIM20', 'FMF']
>>> f.get("yok")
>>> f.get("yok", [])
[]
```

### `first(key, default=None) -> Any`

**Amaç:** Anahtarın **ilk** değerini döndürür; anahtar yoksa `default`.
`set` ile yazılmış bilgileri okumanın en kısa yolu.

```python
>>> f.first("gene_id")
4210
>>> f.first("yok", "-")
'-'
```

### `keys() -> list[str]`

Anahtarları eklenme sırasıyla döndürür.

```python
>>> f.keys()
['alias', 'gene_id']
```

### `items() -> Iterator[(key, list)]`

`(anahtar, değer listesi)` çiftleri üzerinde dolaşır (listeler kopyadır).

```python
for key, values in f.items():
    print(key, values)
```

### `to_dict() -> dict[str, list]`

Tüm bilgileri bağımsız bir sözlük olarak verir; yazdırma, karşılaştırma ve
dışa aktarma için uygundur.

```python
>>> f.to_dict()
{'alias': ['TRIM20', 'FMF'], 'gene_id': [4210]}
```

### `copy() -> Features`

**Amaç:** Bağımsız bir derin kopya döndürür. Kopyada yapılan değişiklik
orijinali etkilemez.

```python
>>> k = f.copy()
>>> k.add("alias", "X")
>>> f.get("alias")
['TRIM20', 'FMF']
```

### `remove(key) -> bool`

Anahtarı tüm değerleriyle siler. Anahtar vardıysa `True`, yoksa `False`.

```python
>>> f.remove("alias"), f.remove("alias")
(True, False)
```

### `clear() -> None`

Tüm bilgileri siler.

### İşleçler

| Yazım | Anlamı |
|---|---|
| `"gene_id" in f` | Anahtar var mı |
| `len(f)` | Anahtar sayısı |
| `for key in f` | Anahtarlar üzerinde dolaşma |
| `repr(f)` | `Keys:alias,gene_id` |

### Paketin kendisinin yazdığı anahtarlar

| Anahtar | Yazan | Değer |
|---|---|---|
| `fasta_header`, `source_file` | `from_fasta_file` | FASTA başlığının tamamı, dosya yolu |
| `orf_start`, `orf_end`, `orf_strand`, `orf_frame` | `DNA.find_orfs` | bkz. 6 |
| `source`, `release`, `gene_id`, `gene_name`, `gene_type`, `chromosome`, `start`, `end`, `chromosome_strand`, `transcript_id`, `transcript_name`, `transcript_type`, `protein_id`, `utr5`, `cds`, `utr3` | `Gencode` | bkz. 12 |
| `source`, `gene_id`, `gene_name`, `description`, `organism`, `chromosome`, `map_location`, `aliases`, `genomic_region`, `transcript_id`, `title`, `predicted`, `protein_id`, `utr5`, `cds`, `utr3` | `Ncbi` | bkz. 13 |

Kendi bilgileriniz için küçük harfli, `snake_case` anahtarlar kullanmanız önerilir.

---

## 5. `Sequence`

**Ne için:** Tüm dizilerin temel sınıfı. Molekül tipi belli olmayan ya da
alfabesi denetlenmeyecek diziler için doğrudan kullanılabilir (`molecule` =
`Molecule.SEQ`). `DNA`, `RNA` ve `PRT` bu sınıftan türer ve bu bölümdeki her
şey onlar için de geçerlidir.

### `Sequence(seq, name="")`

| Parametre | Açıklama |
|---|---|
| `seq` | Dizi metni; büyük harfe çevrilir. |
| `name` | Dizinin adı (FASTA başlığı, CRN kaydı ve `repr` için). |

```python
>>> s = Sequence("atgcgtacg", name="ornek")
>>> s
ornek: ATGCGTACG <9> [SEQ]
```

`repr` biçimi: `ad: dizi <uzunluk> [tip]`. Dizi 80 karakter ya da daha
uzunsa ilk 40 ve son 40 karakter gösterilir:

```python
>>> Sequence("A" * 50 + "C" * 50, "uzun")
uzun: AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC <100> [SEQ]
```

### Özellikler

| Özellik | Okuma | Yazma |
|---|---|---|
| `seq` | Dizi metni | Yeni değer büyük harfe çevrilir ve (alt sınıflarda) alfabe denetlenir |
| `name` | Ad | Serbest |
| `features` | `Features` deposu | Yalnızca `Features` nesnesi atanabilir; aksi halde `TypeError` |
| `molecule` | `Molecule` üyesi (sınıf alanı) | Değiştirmeyin |

```python
>>> s.seq, s.name, s.molecule, len(s)
('ATGCGTACG', 'ornek', Molecule.SEQ, 9)
```

### İndeksleme ve dilimleme: `s[i]`, `s[a:b]`

- `s[i]` tek harfi **`str`** olarak döndürür.
- `s[a:b]` (ve adımlı `s[a:b:c]`) **aynı sınıftan yeni bir dizi** döndürür.
  Yeni dizi aynı adı, root alanlarını (ör. DNA'nın strand'i) ve features'ın
  bağımsız bir kopyasını taşır.
- Adım 0 ise `ValueError`.

```python
>>> s[0]
'A'
>>> s[0:3]
ornek: ATG <3> [SEQ]
```

Biyolojik konumdan dilim almak için (1 tabanlı, iki ucu dahil `start..end`):
`s[start - 1:end]`.

### Karşılaştırma: `==`

- İki dizi, **molekül tipleri ve dizileri aynıysa** eşittir. Ad, features ve
  strand karşılaştırılmaz.
- Bir `str` ile büyük/küçük harf farkı gözetmeden karşılaştırılabilir.

```python
>>> s == "ATGCGTACG"
True
>>> Sequence("ATG") == DNA("ATG")     # tipleri farklı
False
```

### Alt dizi arama: `in`

`str` ya da başka bir dizi aranabilir; büyük/küçük harf fark etmez. Başka bir
tip verilirse `TypeError`.

```python
>>> "CGT" in s, "cgt" in s
(True, True)
```

### Dolaşma ve uzunluk

```python
>>> list(Sequence("ATG"))
['A', 'T', 'G']
>>> len(s)
9
```

### Hash

Diziler değiştirilebilir (`seq` atanabilir) olduğu için **hash'lenemez**;
küme elemanı ya da sözlük anahtarı olamaz (`TypeError: unhashable type`).
Sözlükte saklamak için anahtar olarak adı ya da erişim numarasını kullanın.

### `from_fasta_file(file_name) -> list[Self]` (classmethod)

**Amaç:** FASTA dosyasındaki tüm kayıtları, metodun çağrıldığı sınıfın
nesneleri olarak okur.

| Kural | Ayrıntı |
|---|---|
| Dönüş | **Her zaman liste** (tek kayıt olsa bile). |
| Sıkıştırma | `.gz` uzantılı dosyalar doğrudan okunur. |
| Ad | Başlığın ilk kelimesi. |
| features | `fasta_header` (başlığın tamamı) ve `source_file` (dosya yolu). |
| Atlananlar | Boş satırlar, `;` ile başlayan yorum satırları. |
| Dizi satırları | Birleştirilir; içlerindeki boşluklar silinir. |
| Kodlama | UTF-8; BOM ve Windows satır sonları sorun değildir. |
| Sınıf | `DNA.from_fasta_file` her diziyi DNA alfabesine göre denetler. |

**Hatalar:**

| Durum | Hata |
|---|---|
| Dosya yok | `FileNotFoundError` |
| `>` başlığı olmadan dizi verisi | `ValueError: dosya:satır: '>' header'ı olmadan sequence verisi bulundu` |
| Yalnızca `>` olan başlık | `ValueError: dosya:satır: boş FASTA başlığı` |
| Dizisi olmayan başlık | `ValueError: dosya:satır: '>başlık' bir seq bloğu olmadan bitti` |
| Sınıfın alfabesi dışında karakter | `ValueError: dosya:satır: DNA dizisinde geçersiz karakter(ler): ...` |
| Hiç kayıt yok | `ValueError: dosya: geçerli bir FASTA kaydı bulunamadı` |

```python
>>> prts = PRT.from_fasta_file("ornek.fasta")
>>> prts
[sp|P01308|INS_HUMAN: MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALY...GPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN <110> [PRT] N->C, seq2: MAKT <4> [PRT] N->C]
>>> prts[0].features.to_dict()
{'fasta_header': ['sp|P01308|INS_HUMAN Insulin OS=Homo sapiens'], 'source_file': ['ornek.fasta']}
```

Dosya tipini bilmiyorsanız `Sequence.from_fasta_file` kullanın, sonra
`DNA.from_sequence` / `RNA.from_sequence` / `PRT.from_sequence` ile çevirin.

### `to_fasta(file_name) -> None`

**Amaç:** Diziyi FASTA dosyasına yazar. Dosya varsa üzerine yazılır.

- Başlık: `>AD this sequence contains N feature(s) but not listed`. `N`,
  `fasta_header` ve `source_file` dışındaki features anahtarlarının sayısıdır;
  features'ın içeriği dosyaya **yazılmaz**.
- Dizi 60 karakterlik satırlara bölünür.

**Hatalar:** Ad boşsa ya da boşluk içeriyorsa, veya dizi boşsa `ValueError`
(dosya yazılmaz; çünkü bu durumda dosya aynı haliyle geri okunamaz).

```python
>>> prts[0].name = "INS_HUMAN"
>>> prts[0].to_fasta("ins.fasta")
```

`ins.fasta`:

```text
>INS_HUMAN this sequence contains 0 features but not listed
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```

```python
>>> Sequence("ATG", "BRCA1 exon 2").to_fasta("x.fa")
ValueError: FASTA adı boş olamaz ve boşluk içeremez: 'BRCA1 exon 2'
```

> Tip, strand, genomic_strand ve features FASTA'da saklanamaz. Bunları
> kaybetmeden saklamak için `to_crn` kullanın.

### `to_crn(file_name) -> None`

**Amaç:** Dizinin **tüm içeriğini** (tip, ad, dizi, root alanları, features)
CRN ikili dosyasına yazar. pickle kullanılmaz. Dosya varsa üzerine yazılır.

**Hatalar:** Features'ta CRN'nin desteklemediği tipte bir değer varsa
(ör. `dict`) `TypeError`; bu durumda dosyaya hiçbir şey yazılmaz.

```python
>>> d = DNA("ATGGTGCATCTGACTCCTGAG", "HBB_parca")
>>> d.features.set("gene", "HBB")
>>> d.to_crn("hbb.crn")
```

### `from_crn(file_name) -> Self` (classmethod)

**Amaç:** `to_crn` ile yazılmış dosyayı okur.

- `Sequence.from_crn` dosyadaki tipe uygun sınıfı döndürür: DNA dosyası için
  `DNA`, RNA için `RNA`, protein için `PRT`.
- `DNA.from_crn` gibi bir alt sınıftan çağrılırsa dosyadaki tip o sınıfla
  aynı olmalıdır.

**Hatalar:**

| Durum | Hata |
|---|---|
| Dosya yok | `FileNotFoundError` |
| CRN dosyası değil | `ValueError: ...: CRN dosyası değil` |
| Dosya bozuk | `ValueError: ...: CRN dosyası bozuk (CRC32 uyuşmuyor)` |
| Desteklenmeyen sürüm | `ValueError: ...: desteklenmeyen CRN sürümü: N` |
| Tip uyuşmazlığı | `ValueError: ...: dosyadaki dizi tipi DNA, RNA sınıfı RNA tipi bekliyor` |
| Gen dosyası | `ValueError: ...: bu bir gen dosyası; Gene.from_crn ile açın` |

```python
>>> DNA.from_crn("hbb.crn")
HBB_parca: ATGGTGCATCTGACTCCTGAG <21> [DNA] 5'->3'
>>> type(Sequence.from_crn("hbb.crn")).__name__
'DNA'
>>> DNA.from_crn("hbb.crn").features.to_dict()
{'gene': ['HBB']}
>>> RNA.from_crn("hbb.crn")
ValueError: hbb.crn: dosyadaki dizi tipi DNA, RNA sınıfı RNA tipi bekliyor
```

---

## 6. `DNA`

**Ne için:** DNA dizileri. Tümleme, transkripsiyon ve ORF arama yapar.
`Sequence`'in tüm özelliklerine sahiptir.

### `DNA(seq, name="", strand=Strand.FivePrime)`

| Parametre | Açıklama |
|---|---|
| `seq` | Dizi; büyük harfe çevrilir, `DNA_ALPHABET` ile denetlenir. |
| `name` | Ad. |
| `strand` | Yazılış yönü (bkz. 2.3). |

**Hatalar:** Alfabe dışı karakterde `ValueError`; `strand` bir `Strand` değilse `TypeError`.

```python
>>> d = DNA("ATGGTGCATCTGACTCCTGAG", "HBB_parca")
>>> d
HBB_parca: ATGGTGCATCTGACTCCTGAG <21> [DNA] 5'->3'
>>> DNA("AUG")
ValueError: DNA dizisinde geçersiz karakter(ler): U
>>> DNA("ATG1")
ValueError: DNA dizisinde geçersiz karakter(ler): 1
```

`repr`'ın sonundaki `5'->3'` strand'i gösterir.

### Özellikler

| Özellik | Açıklama |
|---|---|
| `strand` | Okunur ve yazılır; yalnızca `Strand` kabul eder. |

### `from_sequence(sequence) -> Self` (classmethod)

**Amaç:** Bir `Sequence`'i (ya da başka bir dizi nesnesini) aynı ad, dizi ve
features kopyasıyla DNA'ya çevirir. Kaynağın features'ında `strand`
anahtarında bir `Strand` varsa DNA'nın strand'i o olur; yoksa `FivePrime`.

**Hatalar:** DNA alfabesi dışı karakterde `ValueError`.

```python
>>> s = Sequence("ATGGCC", "s")
>>> s.features.set("strand", Strand.ThreePrime)
>>> DNA.from_sequence(s)
s: ATGGCC <6> [DNA] 3'->5'
```

### `complement() -> Self`

**Amaç:** DNA'yı **yerinde** tümleyene çevirir ve strand'i ters çevirir.
`self` döndürür, böylece zincirleme kullanılabilir.

```python
>>> x = DNA("AACC", "x")
>>> x.complement() is x
True
>>> x
x: TTGG <4> [DNA] 3'->5'
```

### `complementary() -> Self`

**Amaç:** Tümleyici zinciri **yeni** bir DNA olarak döndürür; orijinal
değişmez. Yeni DNA aynı adı ve features kopyasını taşır.

```python
>>> d.complementary()
HBB_parca: TACCACGTAGACTGAGGACTC <21> [DNA] 3'->5'
>>> d
HBB_parca: ATGGTGCATCTGACTCCTGAG <21> [DNA] 5'->3'
>>> DNA("ATGNRYKMBVDH").complementary().seq     # IUPAC kodları da tümlenir
'TACNYRMKVBHD'
```

### `reverse_complement() -> Self`

**Amaç:** Karşı zinciri bu diziyle **aynı yönde** okunacak biçimde, yeni bir
nesne olarak döndürür (5'→3' `ATGGCC` → 5'→3' `GGCCAT`). Strand değişmez.
`RNA`'da da vardır.

```python
>>> DNA("ATGGCC", "d").reverse_complement()
d: GGCCAT <6> [DNA] 5'->3'
```

Yalnızca dizinin 5'→3' okunuşu gerekiyorsa (boşluklar atılmış, 3'→5' yazılmışsa
çevrilmiş): `d.seq_5to3()`.

### `transcription(template=True, genomic_strand=GenomicStrand.Unknown) -> RNA`

**Amaç:** DNA'dan RNA üretir (transkripsiyon).

| Parametre | Açıklama |
|---|---|
| `template=True` | Bu DNA **kalıp** (antisense) zincirdir. RNA polimeraz gibi tümleyici RNA üretilir (A→U, T→A, G→C, C→G) ve strand ters çevrilir. |
| `template=False` | Bu DNA **kodlayan** (sense) zincirdir. RNA aynı diziye sahiptir; yalnızca T→U yapılır ve strand korunur. NCBI/GENCODE'dan gelen gen ve CDS dizileri kodlayan zincirdir. |
| `genomic_strand` | Oluşan RNA'nın `genomic_strand` alanı. |

Yeni RNA aynı adı ve features kopyasını taşır.

```python
>>> cds = DNA("ATGGTGCATCTGACTCCTGAGGAGAAGTAA", "HBB_cds")
>>> cds.transcription(template=False)
HBB_cds: AUGGUGCAUCUGACUCCUGAGGAGAAGUAA <30> [RNA] 5'->3' (?)
>>> cds.transcription()                   # bu DNA'yı kalıp kabul eder
HBB_cds: UACCACGUAGACUGAGGACUCCUCUUCAUU <30> [RNA] 3'->5' (?)
>>> cds.complementary().transcription() == cds.transcription(template=False)
True
>>> cds.transcription(template=False, genomic_strand=GenomicStrand.Minus)
HBB_cds: AUGGUGCAUCUGACUCCUGAGGAGAAGUAA <30> [RNA] 5'->3' (-)
```

### `find_orfs(minimum_length_of_orf=-1, genetic_code=GeneticCode.STANDARD, alternative_starts=False) -> list[Self]`

**Amaç:** DNA'nın **iki zincirindeki** tüm geçerli açık okuma çerçevelerini
(ORF) bulur.

**Geçerli ORF:** Bir başlangıç kodonuyla başlayıp aynı okuma çerçevesindeki
ilk stop kodonuyla biten dizi.

| Kural | Ayrıntı |
|---|---|
| Taranan çerçeveler | Bu zincir ve tümleyici zincir, 5'→3' yönünde, 3'er çerçeve (toplam 6). DNA 3'→5' yazılmışsa önce ters çevrilir. |
| Kısmi ORF | Stop kodonuna ulaşmadan dizinin sonuna gelen ORF'ler **alınmaz**. |
| İç içe ORF'ler | Aynı stop'u paylaşan, sonraki başlangıç kodonlarından başlayan ORF'lerin **hepsi** ayrı ORF olarak döner. |
| Boşluklar | `-` karakterleri taramadan önce atlanır. |
| Dairesel diziler | Desteklenmez. |

| Parametre | Açıklama |
|---|---|
| `minimum_length_of_orf` | Stop dahil en küçük ORF uzunluğu (nükleotit). `-1` (varsayılan) ise sınır yok. `-1`'den küçükse `ValueError`. |
| `genetic_code` | Stop (ve alternatif başlangıç) kodonlarının alınacağı kod. |
| `alternative_starts` | `False`: yalnızca ATG. `True`: genetik kodun tüm başlangıç kodonları (standartta TTG, CTG, GTG da). |

**Dönüş:** Konuma göre sıralı liste. Her ORF, **kodlayan yönde 5'→3'
yazılmış** bir DNA'dır; adı `<ad>_ORF<n>`, features kopyasını
taşır. Features'a ek olarak:

| Anahtar | Değer |
|---|---|
| `orf_start`, `orf_end` | ORF'un bu DNA'nın `seq` dizisinde kapladığı bölge (1 tabanlı, ilk nükleotitten uzaklık). Tümleyici zincirdeki ORF'lerde de `orf_start <= orf_end`. |
| `orf_strand` | `"+"` bu zincirde, `"-"` tümleyici zincirde. |
| `orf_frame` | O zincirin 5' ucuna göre okuma çerçevesi (0, 1, 2). |

```python
>>> x = DNA("CCATGGTGCATCTGACTCCTGAGGAGAAGTAAGCTTCTAGGGTTTCATGG", "test")
>>> for orf in x.find_orfs():
...     print(orf, orf.features.to_dict())
test_ORF1: ATGGTGCATCTGACTCCTGAGGAGAAGTAA <30> [DNA] 5'->3' {'orf_start': [3], 'orf_end': [32], 'orf_strand': ['+'], 'orf_frame': [2]}
test_ORF2: ATGAAACCCTAG <12> [DNA] 5'->3' {'orf_start': [37], 'orf_end': [48], 'orf_strand': ['-'], 'orf_frame': [2]}

>>> [o.name for o in x.find_orfs(minimum_length_of_orf=20)]
['test_ORF1']
>>> [o.seq for o in DNA("AACTGAAATAACC", "a").find_orfs(alternative_starts=True)]
['CTGAAATAA']
>>> [o.seq for o in DNA("ATGAAATGAGGGTAA").find_orfs()]
['ATGAAATGA']
>>> [o.seq for o in DNA("ATGAAATGAGGGTAA").find_orfs(genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL)]
['ATGAAATGAGGGTAA', 'ATGAGG']
```

ORF'u proteine çevirmek: `orf.transcription(template=False).translation()`.

Bir `+` ORF'u orijinal diziden geri almak:
`x.seq[orf_start - 1:orf_end]`; `-` ORF için bu bölgenin ters tümleyeni ORF'u verir.

---

## 7. `RNA`

**Ne için:** RNA dizileri. Tümleme, ters transkripsiyon ve translasyon yapar.

### `RNA(seq, name="", strand=Strand.FivePrime, genomic_strand=GenomicStrand.Unknown)`

| Parametre | Açıklama |
|---|---|
| `seq` | Dizi; büyük harfe çevrilir, `RNA_ALPHABET` ile denetlenir (**T kabul edilmez**). |
| `strand` | Yazılış yönü. |
| `genomic_strand` | Transkript edildiği genomik zincir. |

**Hatalar:** Alfabe dışı karakterde `ValueError`; `strand`/`genomic_strand`
yanlış tipteyse `TypeError`.

```python
>>> m = RNA("AUGGUGCAUCUGACUCCUGAGGAGAAGUAA", "HBB_mRNA", genomic_strand=GenomicStrand.Minus)
>>> m
HBB_mRNA: AUGGUGCAUCUGACUCCUGAGGAGAAGUAA <30> [RNA] 5'->3' (-)
```

`repr`'ın sonundaki `(-)` genomic_strand'i gösterir (`?` = bilinmiyor).

### Özellikler

| Özellik | Açıklama |
|---|---|
| `strand` | Okunur/yazılır; yalnızca `Strand`. |
| `genomic_strand` | Okunur/yazılır; yalnızca `GenomicStrand`. |

### `from_sequence(sequence) -> Self` (classmethod)

DNA'daki gibi; ayrıca kaynağın features'ındaki `genomic_strand` anahtarında
bir `GenomicStrand` varsa onu kullanır (yoksa `Unknown`).

### `complement() -> Self` ve `complementary() -> Self`

DNA'dakiyle aynı; yalnızca `strand` ters çevrilir, **`genomic_strand`
değişmez** (RNA'nın geldiği genomik zincir tümlemeyle değişmez).

```python
>>> m.complementary()
HBB_mRNA: UACCACGUAGACUGAGGACUCCUCUUCAUU <30> [RNA] 3'->5' (-)
```

### `reverse_transcription() -> DNA`

**Amaç:** Ters transkripsiyon: RNA'ya tümleyici, antiparalel ilk zincir cDNA'yı
üretir (A→T, U→A, G→C, C→G). Sonuç RNA'nın altına hizalı yazılır; strand ters
çevrilir. Yeni DNA aynı adı ve features kopyasını taşır; DNA'da anlamı olmayan
`genomic_strand` anahtarı kopyadan silinir.

```python
>>> m.reverse_transcription()
HBB_mRNA: TACCACGTAGACTGAGGACTCCTCTTCATT <30> [DNA] 3'->5'
>>> m.reverse_transcription().complementary().seq     # mRNA ile aynı anlamlı ikinci zincir
'ATGGTGCATCTGACTCCTGAGGAGAAGTAA'
>>> m.reverse_transcription().transcription(genomic_strand=GenomicStrand.Minus) == m
True
```

### `translation(genetic_code=GeneticCode.STANDARD) -> PRT`

**Amaç:** RNA'yı seçilen genetik koda göre proteine çevirir.

| Kural | Ayrıntı |
|---|---|
| Başlangıç | Dizi **ilk nükleotitten** başlanarak okunur. AUG, ORF ya da stop **aranmaz**. |
| Stop kodonları | `*` olarak yazılır ve okuma devam eder. |
| Yön | RNA 3'→5' yazılmışsa önce ters çevrilir (ribozom 5'→3' okur). |
| Boşluk | `-` karakterleri atılır. |
| Belirsiz bazlar | N, R, Y gibi kod içeren kodonlar `X` olur. |
| Son | 3'ten kısa kalan son parça okunmaz. |
| Başlangıç kodonu | Özel işlem yapılmaz; her kodon tablodaki amino asidine çevrilir (ör. GUG her yerde V). |
| features | Kopyalanır; proteinde anlamı olmayan `strand` ve `genomic_strand` anahtarları silinir. |

```python
>>> m.translation()
HBB_mRNA: MVHLTPEEK* <10> [PRT] N->C
>>> RNA("AUGUGAAUAAGA").translation().seq
'M*IR'
>>> RNA("AUGUGAAUAAGA").translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL).seq
'MWM*'
>>> RNA("AUGNNNUGGCC").translation().seq       # NNN -> X, sondaki "CC" okunmaz
'MXW'
>>> RNA("CCAUGGUGCAU").translation().seq       # AUG'yi aramaz, ilk nükleotitten okur
'PWC'
>>> RNA(m.seq[::-1], "ters", strand=Strand.ThreePrime).translation().seq
'MVHLTPEEK*'
```

Bir mRNA'nın yalnızca kodlayan bölgesini çevirmek için CDS'i keserek çevirin:

```python
start, end = rna.features.first("cds")     # Gencode/Ncbi bu bilgiyi yazar
protein = rna[start - 1:end].translation()
```

### `match_orf(protein, genetic_code=GeneticCode.STANDARD) -> RNA | None`

**Amaç:** Çevirisi `protein` ile birebir aynı olan ORF'yi (AUG … stop) bulur;
bir transkriptin kodlayan bölgesini (CDS) proteinden doğrulayarak çıkarmanın
yoludur. Yalnızca RNA'nın kendi zincirine bakılır; sondaki stop simgesi iki
tarafta da yok sayılır. Birden çok ORF uyarsa en baştaki seçilir.

| Parametre | Açıklama |
|---|---|
| `protein` | `PRT` ya da tek harfli amino asit dizisi. |
| `genetic_code` | Translasyon için genetik kod. |

**Dönüş:** ORF'yi (stop dahil) 5'→3' yeni bir RNA olarak, uyan yoksa `None`.
Yeni RNA bu RNA'nın adını, `genomic_strand`'ini ve features kopyasını taşır;
`orf_start`, `orf_end` (bu RNA üzerinde 1 tabanlı) ve `cds = (1, uzunluk)`
yazılır, `utr5` ve `utr3` silinir.

```python
>>> m = RNA("GGAUGGUGCAUCUGACUCCUGAGGAGAAGUAAGC", "HBB_mRNA", genomic_strand=GenomicStrand.Minus)
>>> orf = m.match_orf("MVHLTPEEK")
>>> orf
HBB_mRNA: AUGGUGCAUCUGACUCCUGAGGAGAAGUAA <30> [RNA] 5'->3' (-)
>>> orf.features.first("orf_start"), orf.features.first("orf_end")
(3, 32)
>>> m.match_orf("MVHL") is None      # stop'tan önce bitmiyor
True
```

---

## 8. `PRT`

**Ne için:** Protein dizileri (tek harfli amino asit kodları, N→C yönünde).
Strand ve tümleme kavramı yoktur; root alanı yoktur.

### `PRT(seq, name="")`

Dizi büyük harfe çevrilir ve `PRT_ALPHABET` ile denetlenir.

```python
>>> PRT("malwmrllpll", "INS_parca")
INS_parca: MALWMRLLPLL <11> [PRT] N->C
>>> PRT("MALW1")
ValueError: PRT dizisinde geçersiz karakter(ler): 1
```

### `from_sequence(sequence) -> Self` (classmethod)

Bir `Sequence`'i aynı ad, dizi ve features kopyasıyla proteine çevirir.

```python
>>> PRT.from_sequence(Sequence("MAKT", "x"))
x: MAKT <4> [PRT] N->C
```

`Sequence`'ten gelen tüm özellikler (dilimleme, `in`, `to_fasta`, `to_crn`,
`from_crn`, `from_fasta_file` ...) proteinde de çalışır.

---

## 9. `Gene`

**Ne için:** Bir gene ait dizileri bir arada tutan **kapsayıcı**: tek genomik
DNA, istenen sayıda RNA (transkript varyantları) ve istenen sayıda protein.
Bir gen birden çok transkript varyantı, dolayısıyla birden çok protein
izoformu verebilir; `Gene` bunları tek nesnede toplar.

```text
Gene (ör. MEFV)
 ├── dna          : DNA            genomik dizi (tek)
 ├── transcripts  : {id: RNA}
 ├── proteins     : {id: PRT}
 └── features     : Features       gene_id, locus ...
```

> **Önemli:** Transkriptler ve proteinler **birbirinden bağımsız** eklenir.
> Hangi RNA'nın hangi proteini kodladığı gibi ilişkiler kullanıcının
> sorumluluğundadır. `Gencode` ve `Ncbi` bu bağı proteinin
> `transcript_id` özelliğine yazar.

### `Gene(name, dna=None)`

| Parametre | Açıklama |
|---|---|
| `name` | Genin adı (ör. sembol). |
| `dna` | Genomik DNA; verilmezse sonradan `gene.dna = ...` ile atanabilir. |

```python
>>> g = Gene("MEFV", DNA("ATGGCTAAGACC", "MEFV_genomic"))
>>> g.features.set("gene_id", 4210)
```

### Özellikler

| Özellik | Açıklama |
|---|---|
| `name` | Okunur/yazılır. |
| `dna` | Okunur/yazılır; yalnızca `DNA` (aksi halde `TypeError`). Yoksa `None`. |
| `features` | Genin `Features` deposu; yalnızca `Features` atanabilir. |
| `transcripts` | `{kimlik: RNA}`, ekleme sırasıyla; **kopya sözlük** (sözlüğü değiştirmek geni etkilemez, ama içindeki RNA nesneleri aynıdır). |
| `proteins` | `{kimlik: PRT}`, ekleme sırasıyla; kopya sözlük. |

### `add_transcript(rna, transcript_id=None) -> RNA`

**Amaç:** Bir RNA'yı transkript olarak ekler ve onu döndürür. Kimlik
verilmezse RNA'nın adı kullanılır.

**Hatalar:** `rna` bir RNA değilse `TypeError`; kimlik boşsa ya da bu kimlikte
bir transkript zaten varsa `ValueError`.

```python
>>> g.add_transcript(RNA("AUGGCUAAGACCUGA", "NM_000243.3"))
NM_000243.3: AUGGCUAAGACCUGA <15> [RNA] 5'->3' (?)
>>> g.add_transcript(RNA("AUGCCC"), "lnc1")          # adsız RNA: kimliği ayrıca verin
>>> g.add_transcript(RNA("AUG", "NM_000243.3"))
ValueError: MEFV: 'NM_000243.3' transkripti zaten var
>>> g.add_transcript(RNA("AUG"))
ValueError: Transkript kimliği boş olamaz
```

### `transcript(transcript_id) -> RNA`

Kimliği verilen transkripti döndürür; yoksa `KeyError`.

```python
>>> g.transcript("lnc1")
: AUGCCC <6> [RNA] 5'->3' (?)
>>> g.transcript("yok")
KeyError: "MEFV: 'yok' transkripti yok"
```

### `remove_transcript(transcript_id) -> bool`

Transkripti siler; vardıysa `True`, yoksa `False`. Proteinlere dokunmaz.

### `add_protein(prt, protein_id=None) -> PRT`

**Amaç:** Bir proteini ekler ve onu döndürür. Kimlik verilmezse proteinin adı
kullanılır. Transkriptle bağ kurmaz.

**Hatalar:** `prt` bir PRT değilse `TypeError`; kimlik boşsa ya da zaten
varsa `ValueError`.

```python
>>> p = g.add_protein(PRT("MAKT", "NP_000234.1"))
>>> p.features.set("transcript_id", "NM_000243.3")     # bağı siz kurarsınız
>>> g
MEFV: DNA <12>, 2 transkript, 1 protein
>>> g.add_protein(RNA("AUG"), "x")
TypeError: protein bir PRT olmalı, RNA verildi
```

### `protein(protein_id) -> PRT`

Kimliği verilen proteini döndürür; yoksa `KeyError`.

### `remove_protein(protein_id) -> bool`

Proteini siler; vardıysa `True`, yoksa `False`.

### Bir transkriptin proteinini bulmak

```python
tid = "NM_000243.3"
proteinler = [p for p in g.proteins.values() if p.features.first("transcript_id") == tid]
```

### `isoforms(genetic_code=GeneticCode.STANDARD, unique_cds=True, priority=("NM_", "XM_")) -> list[Gene]`

**Amaç:** Geni izoformlarına ayırır: her biri **bir genomik DNA, bir
transkript, bir protein** içeren yeni genler. Kodon kullanımı hesabı
(`crn.cub.CUB`) bu yapıyı ister. Bu gen değişmez.

Her transkript için çevirisi bir proteinle birebir aynı ORF aranır
(`RNA.match_orf`); önce kaynağın o transkripte bağladığı protein
(`transcript_id` özelliği), sonra diğerleri denenir.

| Parametre | Açıklama |
|---|---|
| `genetic_code` | ORF arama ve translasyon için genetik kod. |
| `unique_cds` | `True` ise aynı CDS'i veren izoformlardan yalnızca biri tutulur. |
| `priority` | Tutulacak izoformu seçerken transkript kimliği önekleri, öncelik sırasıyla; eşitlikte ilk eklenen. |

**Dönüş:** Transkriptlerin ekleme sırasıyla izoform genleri. Adları
`<gen>_<transkript>`; features genin kopyasıdır, ayrıca `isoform_of`,
`transcript_id` ve `protein_id` yazılır. Transkript kırpılmaz; ORF konumu
transkriptin `orf_start`/`orf_end` özelliklerindedir.

**Uyarılar ve hatalar:** Proteini olmayan transkriptler, transkripti olmayan
proteinler ve çıkarılan kopya izoformlar `UserWarning` ile bildirilir.
Genomik DNA yoksa `ValueError`.

```python
>>> g = Gene("ORNEK", DNA("ATGAAAGCCTAGGG"))
>>> g.add_transcript(RNA("CCAUGAAAGCCUAGGG"), "T1")
>>> g.add_transcript(RNA("UUAUGAAAGCCUAGAA"), "T2")
>>> g.add_transcript(RNA("GGGCCCUUU"), "T3")
>>> g.add_protein(PRT("MKA"), "P1")
>>> [(i.name, list(i.transcripts), list(i.proteins)) for i in g.isoforms()]
UserWarning: ORNEK: proteini olmayan transkriptler: T3
UserWarning: ORNEK: aynı CDS'i veren izoformlar çıkarıldı (T2 -> T1)
[('ORNEK_T1', ['T1'], ['P1'])]
```

### `to_crn(file_name) -> None`

**Amaç:** Geni tüm içeriğiyle (ad, features, genomik DNA, transkriptler ve
proteinler kimlikleriyle) tek bir CRN dosyasına yazar. Her dizi kendi tipi,
root alanları ve features'ıyla yazılır.

**Hatalar:** Desteklenmeyen tipte bir feature değeri varsa `TypeError`
(dosya yazılmaz).

### `from_crn(file_name) -> Self` (classmethod)

**Amaç:** `Gene.to_crn` ile yazılmış dosyayı okur; tüm diziler kendi
sınıflarıyla ve ekleme sırasıyla geri gelir.

**Hatalar:** `Sequence.from_crn` ile aynı; ayrıca tek bir dizinin dosyası
verilirse `ValueError: ...: gen dosyası değil (dosya tipi: DNA)`.

```python
>>> g.to_crn("mefv.crn")
>>> h = Gene.from_crn("mefv.crn")
>>> h, h.features.to_dict()
(MEFV: DNA <12>, 2 transkript, 1 protein, {'gene_id': [4210]})
>>> Gene.from_crn("hbb.crn")
ValueError: hbb.crn: gen dosyası değil (dosya tipi: DNA)
>>> Sequence.from_crn("mefv.crn")
ValueError: mefv.crn: bu bir gen dosyası; Gene.from_crn ile açın
```

---

## 10. CRN dosya biçimi

CRN, paketin kendi ikili dosya biçimidir. FASTA'dan farkı, nesnenin **her
şeyini** (tip, root alanları, features) kayıpsız saklamasıdır. pickle
kullanılmaz; dosyadan hiçbir kod çalıştırılmaz.

| Özellik | Ayrıntı |
|---|---|
| Tanıma | Dosya `CRN\0` ile başlar. |
| Sürüm | Şu anki sürüm **3**; yalnızca bu sürüm okunur. |
| Bütünlük | Dosyanın sonunda CRC32 vardır; bozulmuş dosya okunmaz. |
| Bir dizi dosyası | Tek bir dizi (`to_crn`/`from_crn`). |
| Bir gen dosyası | Bir `Gene` ve içindeki tüm diziler (`Gene.to_crn`/`Gene.from_crn`). |
| Ayrıntılı tanım | `crn/fileio.py` dosyasının başındaki açıklama. |

**Features ve root alanlarında izin verilen değer tipleri:**

- `None`, `bool`, `int` (her büyüklükte), `float`, `str`, `bytes`
- `Molecule`, `Strand`, `GenomicStrand`, `GeneticCode` enum üyeleri (eski
  dosyalardaki `Maps` değerleri `None` okunur)
- Bunlardan oluşan `tuple` ve `list` (iç içe olabilir)

`dict`, `set`, kendi sınıflarınız gibi diğer tipler `TypeError` verir.

---

## 11. `NcbiSession`

**Ne için:** NCBI E-utilities isteklerini **hız sınırı**, **otomatik yeniden
deneme** ve **bağlantı havuzu** ile yapan HTTP oturumu. `Ncbi` sınıfı
kullanır; `Ncbi(session=...)` ile kendi oturumunuzu verebilirsiniz.

### `NcbiSession(api_key=None, e_mail=None, timeout=30, max_retries=5, pool_connections=10, pool_maxsize=20)`

| Parametre | Açıklama |
|---|---|
| `api_key` | NCBI API anahtarı. Yoksa saniyede **3** istek, varsa **10** istek. |
| `e_mail` | NCBI'ın istediği iletişim adresi; User-Agent'a yazılır. |
| `timeout` | Her istek için zaman aşımı (saniye). |
| `max_retries` | 429, 500, 502, 503, 504 yanıtlarında yeniden deneme sayısı (artan bekleme ile). |
| `pool_connections`, `pool_maxsize` | Bağlantı havuzu ayarları. |

API anahtarını NCBI hesabınızın ayarlarından ücretsiz alabilirsiniz.

```python
>>> NcbiSession()
NcbiSession(Anonim, hız=3 istek/sn, durum=Açık)
>>> NcbiSession(api_key="abcd1234wxyz", e_mail="ben@ornek.com")
NcbiSession(api_key=***wxyz, hız=10 istek/sn, durum=Açık)
```

`repr` API anahtarının yalnızca son 4 karakterini gösterir.

### `get(url, **kwargs) -> requests.Response`

**Amaç:** Hız sınırına uyarak GET isteği yapar. `kwargs` doğrudan
`requests.Session.get`'e iletilir (ör. `params=`, `headers=`).

**Hatalar:** Oturum kapalıysa `RuntimeError`; HTTP hatası ya da ağ hatası
olursa `requests` hatası.

```python
with NcbiSession() as session:
    response = session.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={"db": "gene", "term": "MEFV[Gene Name] AND 9606[Taxonomy ID]", "retmode": "json"},
    )
    print(response.json()["esearchresult"]["idlist"])     # ['4210']
```

### `is_open() -> bool`, `close() -> None`

Oturumun açık olup olmadığını döndürür; oturumu kapatır.

```python
>>> s = NcbiSession(); s.is_open()
True
>>> s.close(); s.is_open()
False
>>> s.get("https://eutils.ncbi.nlm.nih.gov/")
RuntimeError: Kapalı bir oturum üzerinden istek yapılamaz.
```

### Context manager

`with NcbiSession() as session:` bloğunun sonunda oturum kendiliğinden kapanır.

---

## 12. `Gencode`

**Ne için:** GENCODE'un insan ve fare verilerinden genomik DNA, transkript
(RNA), protein ve `Gene` çıkarmak. Veriler bir kez indirilir, `sys/`
klasöründe hazırlanır ve sonraki tüm isteklerde ağa çıkılmadan oradan okunur.

### Çalışma mantığı

1. Bir dizi istendiğinde önce `sys/` içindeki hazır sürüme bakılır.
2. Hazır sürüm yoksa `update()` kendiliğinden çağrılır: dosyalar indirilir,
   MD5 ile doğrulanır ve hazırlanır.
3. İlk istekte tablolar belleğe yüklenir (insan v50 için yaklaşık 5 saniye,
   en çok yaklaşık 1,3 GB bellek); aynı nesneyle sonraki istekler anlıktır
   (ör. `gene()` ~0,03 sn).

### İndirilen dosyalar (insan, v50 örneği)

| Dosya | Boyut | İçerik |
|---|---|---|
| `GRCh38.primary_assembly.genome.fa.gz` | ~806 MB | Genom |
| `gencode.v50.transcripts.fa.gz` | ~175 MB | Tüm transkriptlerin dizileri |
| `gencode.v50.primary_assembly.annotation.gtf.gz` | ~119 MB | Gen ve transkript koordinatları |
| `gencode.v50.pc_transcripts.fa.gz` | ~124 MB | CDS/UTR koordinatları |
| `gencode.v50.pc_translations.fa.gz` | ~20 MB | Protein dizileri |

Hazırlıkla birlikte `sys/gencode/human/v50` klasörü yaklaşık **6 GB**
tutar (açılmış genom ~3 GB).

### `sys/` yapısı

```text
sys/gencode/<tür>/<sürüm>/        ör. sys/gencode/human/v50/
    raw/                          indirilen .gz dosyaları
    genome.fa, genome.fai         açılmış genom ve konum dizini
    transcripts.fa, .idx          transkript dizileri ve dizini
    translations.fa, .idx         protein dizileri ve dizini
    genes.tsv, transcripts.tsv    gen ve transkript koordinatları
    regions.tsv                   UTR5/CDS/UTR3 koordinatları
    manifest.json                 hazırlık tamamlanınca yazılır
```

### `Gencode(species="human", release=None, sys_dir=None, base_url=GENCODE_URL, verbose=True)`

| Parametre | Açıklama |
|---|---|
| `species` | `"human"` ya da `"mouse"`. Başka bir değer `ValueError`. |
| `release` | `"v50"`, `"50"`, fare için `"vM39"`/`"M39"`. `None` ise `sys/` içindeki en yeni hazır sürüm; o da yoksa sunucudaki en son sürüm. |
| `sys_dir` | Veri klasörü; varsayılan proje kökündeki `sys/`. |
| `base_url` | GENCODE sunucusu (yansı sunucu kullanmak için). |
| `verbose` | `True` ise indirme, hazırlık ve yükleme adımları tqdm ilerleme çubuklarıyla gösterilir; `False` tamamen sessizdir. |

```python
>>> gc = Gencode(verbose=False)
>>> Gencode(release="49").release, Gencode(species="mouse", release="M39").release
('v49', 'vM39')
```

### Özellikler

| Özellik | Açıklama | Örnek |
|---|---|---|
| `release` | Kullanılan sürüm | `'v50'` |
| `species_dir` | Türün klasörü | `.../sys/gencode/human` |
| `release_dir` | Sürümün klasörü; sürüm belli değilse `ValueError` | `.../sys/gencode/human/v50` |

### `installed_releases() -> list[str]`

`sys/` içinde hazırlığı tamamlanmış sürümler, eskiden yeniye.

```python
>>> gc.installed_releases()
['v50']
```

### `is_ready() -> bool`

Kullanılan sürümün tüm hazırlık dosyaları `sys/` içinde var mı.

### `latest_release() -> str`

GENCODE sunucusundaki en son sürümü döndürür. **Ağ bağlantısı gerekir.**

```python
>>> gc.latest_release()
'v50'
```

### `update(force=False) -> str`

**Amaç:** Sürümü indirir ve hazırlar; hazırlanan sürümü döndürür.

| Kural | Ayrıntı |
|---|---|
| Hangi sürüm | Nesne bir sürümle oluşturulduysa o sürüm; oluşturulmadıysa sunucudaki **en son** sürüm (sys/ içinde eskisi olsa bile). |
| Zaten hazırsa | Hiçbir şey yapmaz (`force=False`). |
| Doğrulama | Her dosya sunucunun MD5SUMS listesiyle doğrulanır. |
| Tekrar indirme | İndirilmiş ve doğrulanmış dosyalar yeniden indirilmez; bozuk yerel dosya yeniden indirilir. |
| Artımlı hazırlık | Önceden hazırlanmış dosyalar yeniden hazırlanmaz; yalnızca eksikler üretilir. |
| Güvenli yazma | Her çıktı önce `.tmp` olarak yazılır, bitince yerine konur; yarıda kalan iş bir sonraki çağrıda tamamlanır. |
| `force=True` | Tüm dosyalar yeniden indirilir ve yeniden hazırlanır. |

**Hatalar:** MD5 tutmazsa ya da beklenen dosya sunucuda yoksa `RuntimeError`;
ağ hatalarında `urllib` hatası.

```python
gc = Gencode()
gc.update()                      # sunucudaki en son sürümü indirir/hazırlar
Gencode(release="49").update()   # belirli bir sürüm
```

Komut satırından:

```bash
.venv/bin/python -c "from crn.fetch import Gencode; Gencode().update()"
```

### `dna(gene) -> DNA`

**Amaç:** Genin genomik DNA'sını **genin kendi (kodlayan) zincirinde, 5'→3'**
döndürür. (−) zincirdeki genler için referans bölgenin ters tümleyeni alınır.

**`gene`:** Gen sembolü (`"MEFV"`) ya da gene_id (`"ENSG00000103313"`,
sürümlü ya da sürümsüz).

**features:** `source`, `release`, `chromosome`, `start`, `end`,
`chromosome_strand`, `gene_id`, `gene_name`, `gene_type`.

**Hatalar:** Gen yoksa `KeyError`; aynı adda birden çok gen varsa `ValueError`
(gene_id kullanın).

```python
>>> dna = gc.dna("MEFV")
>>> dna
MEFV: AGATCCTGGCAGGAAGGAAGATTGGAGGGGGTGTCTGGAA...CTAATACCAGAACATTTTATCATATGTATTAGCAGTCACT <16563> [DNA] 5'->3'
>>> dna.features.to_dict()
{'source': ['GENCODE'], 'release': ['v50'], 'chromosome': ['chr16'], 'start': [3240392], 'end': [3256954], 'chromosome_strand': ['-'], 'gene_id': ['ENSG00000103313.16'], 'gene_name': ['MEFV'], 'gene_type': ['protein_coding']}
```

### `rna(transcript) -> RNA`

**Amaç:** Transkriptin olgun (splice edilmiş) dizisini RNA olarak döndürür.

**`transcript`:** transcript_id (`"ENST00000219596.6"` ya da sürümsüz
`"ENST00000219596"`) veya transkript adı (`"MEFV-201"`).

**features:** `source`, `release`, `transcript_id`, `transcript_name`,
`transcript_type`, `gene_id`, `gene_name`, `chromosome`, `start`, `end`;
protein kodlayan transkriptlerde ayrıca `utr5`, `cds`, `utr3` (transkript
üzerinde 1 tabanlı, iki ucu dahil `(başlangıç, bitiş)`). `genomic_strand`
genin zincirinden doldurulur.

**Hatalar:** Transkript yoksa ya da dizisi yoksa `KeyError`.

```python
>>> r = gc.rna("MEFV-201")
>>> r
ENST00000219596.6: CUACCAGAAGCCAGACAGCUGGCUCGAGCCUCUCCUGCUC...AAUUGACUGGGCCUCAUGAAUAAAAAAAUUAGAAAAUCUG <3506> [RNA] 5'->3' (-)
>>> r.features.to_dict()
{'source': ['GENCODE'], 'release': ['v50'], 'transcript_id': ['ENST00000219596.6'], 'transcript_name': ['MEFV-201'], 'transcript_type': ['protein_coding'], 'gene_id': ['ENSG00000103313.16'], 'gene_name': ['MEFV'], 'chromosome': ['chr16'], 'start': [3242027], 'end': [3256633], 'utr5': [(1, 46)], 'cds': [(47, 2392)], 'utr3': [(2393, 3506)]}
>>> s, e = r.features.first("cds")
>>> r[s - 1:e].translation().seq[:10]
'MAKTPSDHLL'
```

### `protein(identifier) -> PRT`

**Amaç:** Proteini döndürür.

**`identifier`:** protein_id (`"ENSP00000339639"`), transcript_id
(`"ENST00000339854.9"`) ya da transkript adı; sürüm numarası olmadan da
verilebilir.

**features:** `source`, `release`, `protein_id`, `transcript_id`, `gene_id`,
`transcript_name`, `gene_name`.

**Hatalar:** Protein yoksa (ör. kodlamayan transkript) `KeyError`.

```python
>>> gc.protein("ENST00000339854.9")
ENSP00000339639.4: MAKTPSDHLLSTLEELVPYDFEKFKFKLQNTSVQKEHSRI...TFASCSFSGPLQPIFSPGTRDGGKNTAPLTICPVGGQGPD <601> [PRT] N->C
>>> gc.protein("ENST00000697124.1")
KeyError: "GENCODE v50: 'ENST00000697124.1' için protein bulunamadı"
```

### `region(chromosome, start, end, strand="+", name="") -> DNA`

**Amaç:** Genomdan istenen bölgeyi DNA olarak döndürür.

| Parametre | Açıklama |
|---|---|
| `chromosome` | `"chr1"` ... `"chrX"`, `"chrY"`, `"chrM"` ya da scaffold adı. |
| `start`, `end` | 1 tabanlı, iki ucu dahil (GTF ile aynı). |
| `strand` | `"+"`: referans zincir; `"-"`: ters tümleyen. Sonuç her durumda 5'→3'. |
| `name` | Boşsa `chr16:3242028-3242087(-)` biçiminde. |

**features:** `source`, `release`, `chromosome`, `start`, `end`, `chromosome_strand`.

**Hatalar:** Kromozom yoksa `KeyError`; konumlar dışarıdaysa ya da `strand`
geçersizse `ValueError`.

```python
>>> gc.region("chr16", 3242028, 3242087, "+").seq
'AGATTTTCTAATTTTTTTATTCATGAGGCCCAGTCAATTCTCTTAAGAATCATGGCCCGG'
>>> gc.region("chr16", 3242028, 3242087, "-")
chr16:3242028-3242087(-): CCGGGCCATGATTCTTAAGAGAATTGACTGGGCCTCATGAATAAAAAAATTAGAAAATCT <60> [DNA] 5'->3'
>>> gc.region("chrZ", 1, 2)
KeyError: "GENCODE v50: 'chrZ' kromozomu yok"
>>> gc.region("chr16", 0, 5)
ValueError: chr16:0-5 kromozomun dışında (uzunluk 90338345)
```

### `gene(gene) -> Gene`

**Amaç:** Geni tüm içeriğiyle döndürür: genomik DNA (`dna()`), genin dizisi
olan **tüm** transkriptleri (kimlik = transcript_id) ve protein kodlayan
transkriptlerin proteinleri (kimlik = protein_id; bağ proteinin
`transcript_id` özelliğinde). Genin features'ı: `source`, `release`,
`gene_id`, `gene_name`, `gene_type`, `chromosome`, `start`, `end`,
`chromosome_strand`.

```python
>>> g = gc.gene("MEFV")
>>> g
MEFV: DNA <16563>, 17 transkript, 15 protein
>>> list(g.transcripts)[:4]
['ENST00000339854.9', 'ENST00000956137.1', 'ENST00001141017.1', 'ENST00000219596.6']
>>> gc.gene("ENSG00000103313")          # gene_id ile de
MEFV: DNA <16563>, 17 transkript, 15 protein
```

### Gen tipine göre liste (ör. rRNA genleri)

Hazır bir metot yoktur; `genes.tsv` doğrudan okunabilir:

```python
import csv
gc = Gencode()
gc.gene("MEFV")                      # veri yoksa indirir/hazırlar
with open(gc.release_dir / "genes.tsv") as f:
    rrna = [r for r in csv.DictReader(f, delimiter="\t")
            if r["gene_type"] in ("rRNA", "Mt_rRNA")]
for r in rrna:
    print(r["gene_name"], r["gene_id"], r["chromosome"], r["start"], r["end"], r["strand"])
```

`genes.tsv` sütunları: `gene_id, gene_name, gene_type, chromosome, start, end, strand`.
`transcripts.tsv` sütunları: `transcript_id, gene_id, transcript_name, transcript_type, chromosome, start, end, strand`.

---

## 13. `Ncbi`

**Ne için:** NCBI'dan (RefSeq) genomik DNA, transkript, protein ve `Gene`
çıkarmak. `Gencode` ile aynı mantıkla çalışır; ancak tüm veritabanını değil,
**istenen genin** kayıtlarını indirir ve `sys/ncbi/` içinde saklar.

### Çalışma mantığı

1. Gen sembolü NCBI Gene ID'sine çevrilir (eşleşme saklanır).
2. Genin kaydı (`sys/ncbi/genes/<gene_id>.json`) ve bağlı tüm dosyalar varsa
   **ağa çıkılmaz**.
3. Yoksa `update()` kendiliğinden çağrılır.

İndirilenler: gen özeti, genin genomik bölgesi, gene bağlı tüm RefSeq RNA'ları
(NM_, NR_, varsa XM_, XR_) ile özellik tabloları ve bu RNA'ların kodladığı
proteinler.

### `sys/` yapısı

```text
sys/ncbi/
    symbols/<organizma>.json        sembol -> Gene ID
    genes/<gene_id>.json            özet + RNA/protein listesi (update sonunda yazılır)
    nuccore/<erişim>.fa, .ft        transkript dizisi ve özellik tablosu
    protein/<erişim>.fa             protein dizisi
    region/<erişim>_<baş>_<son>_<plus|minus>.fa   genomik bölge
```

### `Ncbi(organism="human", session=None, sys_dir=None, verbose=True)`

| Parametre | Açıklama |
|---|---|
| `organism` | `human`, `mouse`, `rat`, `zebrafish`, `fruit_fly`, `worm`. Başka bir değer `ValueError`. Yalnızca sembol aramasını etkiler. |
| `session` | `NcbiSession`; verilmezse API anahtarsız bir oturum açılır. |
| `sys_dir` | Veri klasörü; varsayılan proje kökündeki `sys/`. |
| `verbose` | İlerleme çubukları. |

```python
ncbi = Ncbi()
ncbi = Ncbi(session=NcbiSession(api_key="...", e_mail="ben@ornek.com"))   # daha hızlı
ncbi = Ncbi(organism="mouse")
```

### `ncbi_dir` (özellik)

`sys/ncbi` klasörünün yolu.

### `search_gene(gene) -> str`

**Amaç:** Gen sembolünü NCBI Gene ID'sine çevirir. Rakamlardan oluşan değer
zaten Gene ID kabul edilir. Bulunan eşleşme saklanır; aynı sembol tekrar
sorulduğunda ağa çıkılmaz. Yalnızca yürürlükteki (silinmemiş) kayıtlar aranır.

**Hatalar:** Sembol bu organizmada yoksa `KeyError`.

```python
>>> ncbi.search_gene("MEFV"), ncbi.search_gene("4210")
('4210', '4210')
>>> ncbi.gene("BOYLEBIRGENYOK123")
KeyError: "NCBI: 'BOYLEBIRGENYOK123' geni human organizmasında bulunamadı"
```

### `update(gene, force=False) -> str`

**Amaç:** Genin tüm kayıtlarını indirip `sys/` içine yazar; Gene ID'yi döndürür.
`sys/` içinde zaten bulunan RNA, protein ve bölge kayıtları yeniden
indirilmez; `force=True` ise hepsi yenilenir. En sonda gen kaydı yazılır; bu
dosya yoksa gen hazır sayılmaz.

```python
ncbi.update("MEFV")               # gerekmez; ilk istekte kendiliğinden yapılır
ncbi.update("MEFV", force=True)   # NCBI'daki güncel kayıtlarla yenile
```

### `dna(gene) -> DNA`

**Amaç:** Genin genomik DNA'sını genin kendi zincirinde, 5'→3' döndürür.

**features:** `source`, `chromosome` (erişim, ör. `NC_000016.10`), `start`,
`end`, `chromosome_strand`, `gene_id`, `gene_name`, `description`.

**Hatalar:** Gen ya da genomik koordinatı yoksa `KeyError`.

```python
>>> ncbi.dna("MEFV").features.to_dict()
{'source': ['NCBI'], 'chromosome': ['NC_000016.10'], 'start': [3242027], 'end': [3256633], 'chromosome_strand': ['-'], 'gene_id': ['4210'], 'gene_name': ['MEFV'], 'description': ['MEFV innate immunity regulator, pyrin']}
```

### `rna(accession, genomic_strand=GenomicStrand.Unknown) -> RNA`

**Amaç:** RefSeq transkriptini RNA olarak döndürür (`sys/` içinde yoksa
indirir). Tek başına çağrıldığında genin zinciri bilinmediği için
`genomic_strand` verilmezse `Unknown` olur; `gene()` bunu doldurur.

**`accession`:** `"NM_000243.3"` ya da sürümsüz `"NM_000243"` (sürümsüzde
`sys/` içindeki en yeni sürüm, yoksa NCBI'daki güncel sürüm).

**features:** `source`, `transcript_id`, `title` (NCBI başlığı), `predicted`
(XM_/XR_ ise `True`), `gene_id`, `gene_name`, `protein_id`; protein
kodlayanlarda `utr5`, `cds`, `utr3`.

```python
>>> r = ncbi.rna("NM_000243")
>>> r.features.to_dict()
{'source': ['NCBI'], 'transcript_id': ['NM_000243.3'], 'title': ['Homo sapiens MEFV innate immunity regulator, pyrin (MEFV), transcript variant 1, mRNA'], 'predicted': [False], 'gene_id': ['4210'], 'gene_name': ['MEFV'], 'protein_id': ['NP_000234.1'], 'utr5': [(1, 46)], 'cds': [(47, 2392)], 'utr3': [(2393, 3506)]}
```

### `protein(accession) -> PRT`

**Amaç:** Proteini döndürür.

**`accession`:** Protein erişimi (NP_, XP_) ya da proteini kodlayan
transkriptin erişimi (NM_, XM_); sürümsüz de olur. Transkript erişimiyle
istendiğinde proteinin features'ına `transcript_id` de yazılır.

**features:** `source`, `protein_id`, `title`, (varsa) `transcript_id`.

**Hatalar:** Protein yoksa ya da transkript protein kodlamıyorsa `KeyError`.

```python
>>> ncbi.protein("NM_000243")
NP_000234.1: MAKTPSDHLLSTLEELVPYDFEKFKFKLQNTSVQKEHSRI...TFASCSFSGPLQPIFSPGTRDGGKNTAPLTICPVGGQGPD <781> [PRT] N->C
>>> ncbi.protein("NR_002819.5")
KeyError: 'NCBI: NR_002819.5 protein kodlamıyor'
```

### `region(accession, start, end, strand="+", name="") -> DNA`

**Amaç:** Bir genomik bölgeyi DNA olarak döndürür (`sys/` içinde yoksa
indirir). Parametreler `Gencode.region` ile aynıdır; kromozom adı yerine
NCBI erişimi (`NC_000016.10`) verilir.

**Hatalar:** `strand` ya da konumlar geçersizse `ValueError`.

```python
>>> ncbi.region("NC_000016.10", 3242027, 3242086, "+").seq
'CAGATTTTCTAATTTTTTTATTCATGAGGCCCAGTCAATTCTCTTAAGAATCATGGCCCG'
```

### `gene(gene) -> Gene`

**Amaç:** Geni tüm içeriğiyle döndürür: genomik DNA, gene bağlı tüm RefSeq
RNA'ları (kimlik = erişim; önce doğrulanmış NM_/NR_, sonra tahmini XM_/XR_)
ve kodladıkları proteinler (kimlik = erişim; proteinin `transcript_id`
özelliği doldurulur). Genin features'ı: `source`, `gene_id`, `gene_name`,
`description`, `organism`, `chromosome`, `map_location`, `aliases` (her biri
ayrı değer), `genomic_region` (`(erişim, başlangıç, bitiş, zincir)`).

```python
>>> g = ncbi.gene("MEFV")
>>> g
MEFV: DNA <14607>, 2 transkript, 2 protein
>>> g.features.to_dict()
{'source': ['NCBI'], 'gene_id': ['4210'], 'gene_name': ['MEFV'], 'description': ['MEFV innate immunity regulator, pyrin'], 'organism': ['Homo sapiens'], 'chromosome': ['16'], 'map_location': ['16p13.3'], 'aliases': ['FMF', 'MEF', 'PAAND', 'TRIM20'], 'genomic_region': [('NC_000016.10', 3242027, 3256633, '-')]}
>>> ncbi.gene("MALAT1")
MALAT1: DNA <8779>, 3 transkript, 0 protein
```

---

## 14. Uçtan uca senaryolar

### 14.1. Santral dogma

```python
from crn import DNA, GeneticCode, GenomicStrand

kodlayan = DNA("ATGGTGCATCTGACTCCTGAGGAGAAGTAA", "HBB_cds")

mrna = kodlayan.transcription(template=False, genomic_strand=GenomicStrand.Minus)   # DNA -> RNA
protein = mrna.translation()                                                        # RNA -> protein
cdna = mrna.reverse_transcription()                                                 # RNA -> DNA

print(mrna)       # HBB_cds: AUGGUGCAUCUGACUCCUGAGGAGAAGUAA <30> [RNA] 5'->3' (-)
print(protein)    # HBB_cds: MVHLTPEEK* <10> [PRT] N->C
print(cdna)       # HBB_cds: TACCACGTAGACTGAGGACTCCTCTTCATT <30> [DNA] 3'->5'
```

### 14.2. Bir genin tüm izoformlarının CDS'ini çevirip doğrulama

```python
from crn.fetch import Gencode

gc = Gencode()
g = gc.gene("MEFV")
for tid, rna in g.transcripts.items():
    cds = rna.features.first("cds")
    if cds is None:
        print(tid, rna.features.first("transcript_type"), "kodlamıyor")
        continue
    start, end = cds
    protein = rna[start - 1:end].translation()
    resmi = [p for p in g.proteins.values() if p.features.first("transcript_id") == tid]
    print(tid, len(protein) - 1, "aa", protein.seq.rstrip("*") == resmi[0].seq)
```

### 14.3. CDS koordinatı bilinmeyen bir mRNA'da proteini bulmak

```python
rna = gc.rna("ENST00000339854.9")
hedef = gc.protein("ENST00000339854.9").seq
dna = rna.reverse_transcription().complementary()        # mRNA ile aynı yönde DNA
for orf in dna.find_orfs():
    if orf.features.first("orf_strand") == "+" and \
       orf.transcription(template=False).translation().seq.rstrip("*") == hedef:
        print(orf.features.first("orf_start"), orf.features.first("orf_end"))   # 368 2173
```

### 14.4. Genomik DNA'da ORF taramak

```python
dna = gc.dna("MEFV")
uzunlar = dna.find_orfs(minimum_length_of_orf=300)
for orf in uzunlar[:5]:
    f = orf.features
    print(orf.name, f.first("orf_start"), f.first("orf_end"), f.first("orf_strand"), len(orf))
```

### 14.5. GENCODE ve NCBI'ı karşılaştırmak

```python
from crn.fetch import Gencode, Ncbi

e = Gencode().gene("MEFV")
n = Ncbi().gene("MEFV")
print(e)    # MEFV: DNA <16563>, 17 transkript, 15 protein
print(n)    # MEFV: DNA <14607>, 2 transkript, 2 protein
print(e.dna.features.first("start"), n.dna.features.first("start"))
```

İki kaynağın gen sınırları ve transkript sayıları farklıdır (bkz. 16).
Her nesnenin `source` özelliği kaynağını gösterir.

### 14.6. Kaydetmek ve geri okumak

```python
from crn import Gene, RNA

g = Gencode().gene("MEFV")
g.to_crn("MEFV_gencode.crn")
g2 = Gene.from_crn("MEFV_gencode.crn")

r = g2.transcripts["ENST00000219596.6"]
r.to_crn("MEFV-201.crn")
RNA.from_crn("MEFV-201.crn").features.first("cds")    # (47, 2392)

r.to_fasta("MEFV-201.fasta")                          # yalnızca ad ve dizi
```

### 14.7. Bir gen listesini toplu indirip kaydetmek

```python
from pathlib import Path
from crn.fetch import Ncbi

ncbi = Ncbi(verbose=False)
cikti = Path("genler"); cikti.mkdir(exist_ok=True)
for sembol in ["MEFV", "TP53", "BRCA1", "CFTR"]:
    hedef = cikti / f"{sembol}.crn"
    if hedef.exists():
        continue
    try:
        ncbi.gene(sembol).to_crn(hedef)
    except KeyError as e:
        print(sembol, "atlandı:", e)
```

---

## 15. Hata mesajları

| Mesaj (özet) | Tip | Nedeni | Çözüm |
|---|---|---|---|
| `Anahtar str olmalı` / `Anahtar boş olamaz` | TypeError / ValueError | Features anahtarı geçersiz | Boş olmayan `str` kullanın |
| `X dizisinde geçersiz karakter(ler): ...` | ValueError | Sınıfın alfabesi dışında karakter | Diziyi temizleyin ya da `Sequence` kullanın |
| `strand bir Strand olmalı` | TypeError | Yanlış tip | `Strand.FivePrime` gibi enum verin |
| `features bir Features olmalı` | TypeError | `features`'a başka tip atandı | `Features()` ya da `.copy()` atayın |
| `FASTA adı boş olamaz ve boşluk içeremez` | ValueError | `to_fasta` için ad uygun değil | Adı tek kelime yapın |
| `boş dizi FASTA olarak yazılamaz` | ValueError | Dizi boş | — |
| `'>' header'ı olmadan sequence verisi` | ValueError | Bozuk FASTA | Dosyayı düzeltin |
| `CRN dosyası değil` / `bozuk` / `desteklenmeyen CRN sürümü` | ValueError | Dosya CRN değil, bozuk ya da eski | Dosyayı yeniden üretin |
| `dosyadaki dizi tipi X, Y sınıfı Y tipi bekliyor` | ValueError | Yanlış sınıfla `from_crn` | `Sequence.from_crn` ya da doğru sınıfı kullanın |
| `bu bir gen dosyası; Gene.from_crn ile açın` | ValueError | Gen dosyası dizi olarak açıldı | `Gene.from_crn` |
| `gen dosyası değil` | ValueError | Dizi dosyası gen olarak açıldı | `Sequence.from_crn` |
| `CRN dosyasına yazılamayan değer tipi` | TypeError | Features'ta desteklenmeyen tip | Değeri `str`/`tuple` gibi bir tipe çevirin |
| `transkripti zaten var` / `proteini zaten var` | ValueError | Aynı kimlik ikinci kez | Farklı kimlik verin ya da önce silin |
| `minimum_length_of_orf -1 ya da 0'dan büyük olmalı` | ValueError | Geçersiz sınır | `-1` ya da ≥0 |
| `... geni bulunamadı` / `transkripti bulunamadı` / `protein bulunamadı` | KeyError | Kayıt yok | Kimliği kontrol edin |
| `adında birden çok gen var` | ValueError | Sembol birden çok gene karşılık geliyor | gene_id kullanın |
| `kromozomun dışında` / `kromozomu yok` | ValueError / KeyError | Bölge hatalı | Koordinatları kontrol edin |
| `MD5 doğrulaması başarısız` | RuntimeError | İndirme bozuk | `update()` tekrar çalıştırın |
| `Desteklenmeyen tür` / `Desteklenmeyen organizma` | ValueError | Listede olmayan tür | Desteklenen değerlerden birini seçin |
| `Kapalı bir oturum üzerinden istek yapılamaz` | RuntimeError | Kapatılmış `NcbiSession` | Yeni oturum açın |

---

## 16. Sınırlılıklar

### Diziler ve dönüşümler

1. **FASTA kayıplıdır.** `to_fasta` yalnızca adı ve diziyi yazar; tip,
   `strand`, `genomic_strand` ve features dosyaya geçmez. FASTA'dan
   okunan her dizi `FivePrime` kabul edilir; 3'→5' yazılmış bir dizi FASTA'ya
   yazılıp okunursa yönü yanlış yorumlanır. Kayıpsız saklama için CRN kullanın.
2. **FASTA adı tek kelime olmalıdır.** Boşluklu adla `to_fasta` hata verir;
   `from_fasta_file` adı başlığın ilk kelimesi olarak alır.
3. **Eşitlik yalnızca tip ve diziye bakar.** Ad, features ve strand
   karşılaştırılmaz; aynı harf dizisine sahip 5'→3' ve 3'→5' iki DNA eşit
   sayılır.
4. **Diziler hash'lenemez;** küme ya da sözlük anahtarı olamaz.
5. **Ters dilim yönü değiştirmez.** `s[::-1]` harfleri ters çevirir ama
   `strand` aynı kalır; yönü de değiştirmek gerekiyorsa `strand`'i elle atayın.
6. **Alfabe kısıtları:** RNA'da T kabul edilmez (T içeren RNA dosyaları
   `RNA.from_fasta_file` ile okunamaz; `Sequence` ile okuyup T→U çevirin).
   Değiştirilmiş nükleotitler (ör. inozin, psödouridin) desteklenmez.
7. **Küçük harf bilgisi kaybolur.** Tüm diziler büyük harfe çevrildiği için
   tekrar bölgelerini gösteren küçük harf maskesi (soft-mask) korunmaz.
8. **Belirsiz kodonlar her zaman `X` olur.** `UAR` (her durumda stop) ya da
   `GCN` (her durumda alanin) gibi aslında belirli olan kodonlar da `X` olur;
   `find_orfs` belirsiz baz içeren stop kodonlarını stop olarak tanımaz.
9. **Yalnızca iki genetik kod vardır:** standart (NCBI 1) ve omurgalı
   mitokondriyal (NCBI 2).
10. **`translation` ilk nükleotitten okur.** Başlangıç kodonu, ORF ya da stop
    aramaz; alternatif başlangıç kodonunu (ör. CUG) M'ye çevirmez. Selenosistein
    (U) ve pirolizin (O) hiçbir zaman üretilmez.
11. **`find_orfs`:** dairesel diziler (plazmid, mitokondriyal DNA) desteklenmez;
    başlangıç noktasının üzerinden geçen ORF'ler bulunmaz. İç içe ORF'lerin
    hepsi döner (her stop için yalnızca en uzununu seçme seçeneği yoktur).
    Kısmi ORF'ler alınmaz. Uzun genomik dizilerde çok sayıda kısa ORF dönebilir;
    `minimum_length_of_orf` kullanın.
12. **`genomic_strand` yalnızca RNA'dadır.** DNA'nın kromozomdaki zinciri
    yalnızca features'ta (`chromosome_strand`) tutulur.

### Gene ve CRN

13. **`Gene` ilişki denetlemez.** Transkript–protein bağı, genomik DNA ile
    transkriptlerin uyumu gibi hiçbir tutarlılık kontrol edilmez; silme
    işlemleri birbirini etkilemez. Ekson/intron yapısı ve splicing modellenmez.
14. **CRN değer tipleri sınırlıdır** (bkz. 10); `dict`, `set` ve kendi
    sınıflarınız saklanamaz.
15. **CRN sürümü:** Yalnızca sürüm 3 okunur. Geliştirme sırasında daha eski
    biçimlerle yazılmış dosyalar (ör. gen dosyasının transkript–protein bağlı ilk
    biçimi) açılmaz; yeniden üretilmeleri gerekir.
16. **Bir CRN dosyası tek bir dizi ya da tek bir gen içerir;** birden çok
    bağımsız diziyi tek dosyada saklama biçimi yoktur.

### Gencode

17. **Yalnızca insan ve fare** desteklenir.
18. **Disk ve bellek:** Bir insan sürümü yaklaşık 6 GB disk kaplar; ilk
    istekte tablolar belleğe yüklenir (~5 sn, ~1,3 GB).
19. **İlk indirme ~1,2 GB'tır ve kaldığı yerden devam etmez;** yarıda kesilen
    dosya bir sonraki `update()`'te baştan indirilir (tamamlanmış dosyalar
    yeniden indirilmez).
20. **`latest_release()` ve sürüm verilmemiş `update()` ağ gerektirir.**
21. **Dizisi olmayan transkriptler:** GENCODE yalnızca referans
    kromozomlardaki transkriptlerin dizisini verir; scaffold'lardaki
    transkriptler `gene()` sonucuna eklenmez, `rna()` ile istenirse `KeyError`
    verir.
22. **Aynı adlı genler:** PAR bölgesindeki genlerde X kopyası seçilir; başka
    çoklu eşleşmelerde gene_id gerekir. Bazı rRNA genlerinin adı resmi sembol
    değildir (ör. birden çok gende tekrarlanan `5_8S_rRNA`).
23. **18S, 28S ve 45S rRNA genleri** GRCh38 birincil montajında yer almadığı
    için GENCODE'da yoktur.
24. **Başlangıcı eksik transkriptler:** `utr5` olmayan ve `cds`'i 1'den
    başlayan transkriptlerde CDS AUG ile başlamayabilir ve ilk kodon çerçeveden
    kayık olabilir.
25. **Gen tipine göre listeleme** için hazır metot yoktur (`genes.tsv` okunur).

### Ncbi

26. **Yalnızca RefSeq** RNA'ları ve proteinleri getirilir; GenBank kayıtları
    alınmaz. Altı organizma desteklenir (sembol araması için).
27. **İlk istek ağ gerektirir** ve NCBI'ın hız sınırına tabidir (anahtarsız
    3, anahtarla 10 istek/sn). Uzun gen listelerinde API anahtarı kullanın.
28. **Önbellek kendiliğinden yenilenmez.** Bir gen bir kez indirildikten sonra
    NCBI'daki değişiklikler `update(gene, force=True)` çağrılana kadar
    görünmez. Sembol → Gene ID eşlemeleri hiç yenilenmez; değiştiyse
    `sys/ncbi/symbols/<organizma>.json` dosyasını silin.
29. **Bölge önbelleği koordinata bağlıdır;** aynı bölgenin farklı
    koordinatlarla istenmesi yeni indirme yapar.
30. **Gen tipi bilgisi yoktur** (NCBI özeti bunu vermez); `gene_type`
    features'ı yazılmaz.
31. **Tek başına çağrılarda bağlam eksik olabilir:** `rna()` doğrudan
    çağrıldığında `genomic_strand` `Unknown` olur; `protein()` NP_ erişimiyle
    çağrıldığında `transcript_id` yazılmaz. Bu bilgiler `gene()` ile tam gelir.
32. **Yalnızca ilk CDS okunur;** birden çok CDS içeren nadir transkriptlerde
    (ör. çift kodlayan) ikincisi alınmaz.

### Kaynaklar arası farklar

33. **GENCODE ve NCBI aynı gen için farklı sonuç verir.** Gen sınırları
    farklıdır (MEFV: GENCODE 16 563 nt, NCBI 14 607 nt); transkript sayıları
    farklıdır (MEFV: 17'ye karşı 2); kimlikler farklıdır (ENST/ENSP'ye karşı
    NM/NP). Koordinat sistemleri de farklıdır: GENCODE kromozom adı
    (`chr16`), NCBI erişim (`NC_000016.10`) kullanır.

### Genel

34. **Eşzamanlılık:** `NcbiSession`'ın hız sınırlayıcısı iş parçacığı
    güvenlidir; diğer sınıflar aynı nesnenin birden çok iş parçacığından
    kullanılması için tasarlanmamıştır. Aynı `sys/` klasörüne aynı anda iki
    süreçten `update()` çalıştırmayın.
35. **Python 3.12 ve üstü** gerekir.
