# `crn.seq`: dizi sınıfları

Kaynak: [`crn/seq.py`](../crn/seq.py)

```text
Sequence                 ad, dizi, features; alfabe denetimi; FASTA/CRN
 ├── NucleicAcid         strand (yazılış yönü), tümleme
 │    ├── DNA            transcription, find_orfs
 │    └── RNA            genomic_strand, translation, reverse_transcription, match_orf
 └── PRT                 protein (N->C)
```

Sınıfın kendisinde yalnızca dizinin biyolojik anlamını değiştiren alanlar bulunur: `strand` ve RNA'da `genomic_strand`. Bunlara kök alan (root field) denir. Geri kalan her bilgi [`features`](features.md) içindedir. Yeni nesne döndüren işlemler (dilimleme, `complementary`, `transcription`) features'ın bağımsız bir kopyasını taşır.

<!-- toc -->
- [1. Numaralandırmalar](#1-numaralandırmalar): [`Molecule`](#molecule) · [`Strand`](#strand) · [`flipped`](#strandflipped) · [`GenomicStrand`](#genomicstrand)
- [2. `Sequence`](#2-sequence): [`Sequence`](#sequence) · [`seq`](#sequenceseq) · [`features`](#sequencefeatures) · [`root_fields`](#sequenceroot_fields) · [`from_sequence`](#sequencefrom_sequence) · [`sliding`](#sequencesliding) · [`from_fasta_file`](#sequencefrom_fasta_file) · [`to_fasta`](#sequenceto_fasta) · [`to_crn`](#sequenceto_crn) · [`from_crn`](#sequencefrom_crn) · [`Python protokolleri`](#python-protokolleri)
- [3. `NucleicAcid`](#3-nucleicacid): [`NucleicAcid`](#nucleicacid) · [`strand`](#nucleicacidstrand) · [`seq_5to3`](#nucleicacidseq_5to3) · [`complement`](#nucleicacidcomplement) · [`complementary`](#nucleicacidcomplementary) · [`reverse_complement`](#nucleicacidreverse_complement)
- [4. `DNA`](#4-dna): [`DNA`](#dna) · [`transcription`](#dnatranscription) · [`find_orfs`](#dnafind_orfs)
- [5. `RNA`](#5-rna): [`RNA`](#rna) · [`genomic_strand`](#rnagenomic_strand) · [`reverse_transcription`](#rnareverse_transcription) · [`translation`](#rnatranslation) · [`match_orf`](#rnamatch_orf)
- [6. `PRT`](#6-prt): [`PRT`](#prt)
<!-- /toc -->

## 1. Numaralandırmalar

### Molecule

Dizinin tipini gösterir: `DNA`, `RNA`, `PRT`, `SEQ`. CRN dosyasındaki tip etiketi budur.

```python
from crn import DNA, RNA, PRT, Sequence, Molecule, Strand, GenomicStrand, GeneticCode

print(DNA("ATG").molecule, Sequence("xyz").molecule)  # -> Molecule.DNA Molecule.SEQ
```

### Strand

Dizinin yazılış yönüdür: `Strand.FivePrime` 5'→3', `Strand.ThreePrime` 3'→5'. Harf dizisiyle yönü ayrı tutar: `DNA.complement` hizalı tümleyiciyi yazar ve yönü çevirir.

### Strand.flipped

`Strand.flipped  (özellik)`

Karşı yönü verir.

```python
print(Strand.FivePrime.flipped)  # -> Strand.ThreePrime
```

### GenomicStrand

RNA'nın transkript edildiği genomik zincirdir (GFF/GTF/BED notasyonu): `Plus` (`+`), `Minus` (`-`), `Unstranded` (`.`), `Unknown` (`?`). Yalnızca RNA özelliğidir ve yazılış yönünden bağımsızdır. MEFV 16p13.3'te `-` zincirdedir; mRNA'sı `-` zincirin 5'→3' dizisidir ama `strand` değeri yine `FivePrime` olur.

## 2. `Sequence`

### Sequence

`Sequence(seq, name='')`

Adı, dizisi ve features'ı olan genel dizidir. `seq` büyük harfe çevrilir, `name` isteğe bağlıdır. Alfabe denetimi yapmaz; DNA, RNA ya da protein olduğu bilinmeyen bir FASTA okunurken işe yarar. Alt sınıflar `alphabet`, `molecule` ve `root_types` sınıf değişkenlerini tanımlar. `SEQUENCE_CLASSES` sözlüğü her `Molecule` değerini sınıfına bağlar ve CRN okuyucusu doğru sınıfı buradan bulur; alt sınıf tanımlanırken `Sequence.__init_subclass__` kaydı kendiliğinden yapar.

```python
s = Sequence("acgt-x", name="ornek")
print(s)  # -> ornek: ACGT-X <6> [SEQ]
print(s.molecule, s.alphabet, s.root_types)  # -> Molecule.SEQ None {}
```

### Sequence.seq

`Sequence.seq  (özellik)`

Dizi metnidir. Atama yapılırsa alfabe yeniden denetlenir.

```python
d = DNA("ATG", name="d")
d.seq = "attg"
print(d.seq)  # -> ATTG
```

### Sequence.features

`Sequence.features  (özellik)`

Dizinin `Features` deposudur. Yalnızca `Features` nesnesi atanabilir, başka tip `TypeError` verir.

```python
d.features.add("gene_name", "ornek")
print(d.features.first("gene_name"))  # -> ornek
```

### Sequence.root_fields

`Sequence.root_fields()`

Kök alanları ve değerlerini sözlük olarak verir. `Sequence` ve `PRT` için boştur.

```python
print(DNA("ATG").root_fields())  # -> {'strand': <Strand.FivePrime: ("5'", "3'", 'Five prime')>}
```

### Sequence.from_sequence

`sınıf metodu: Sequence.from_sequence(sequence)`

Başka bir diziden (`sequence`) aynı ad, dizi ve features kopyasıyla bu sınıftan yeni nesne üretir. En sık kullanım FASTA'dan genel `Sequence` olarak okunan kaydı tipine çevirmektir. Kök alanlar features'ta doğru tipte değer varsa oradan, yoksa varsayılandan alınır. Dizide hedef sınıfın alfabesinde olmayan harf varsa `ValueError` verir.

```python
genel = Sequence("ATGC", name="g")
print(DNA.from_sequence(genel))  # -> g: ATGC <4> [DNA] 5'->3'
```

### Sequence.sliding

`Sequence.sliding(window_size, step_size=1)`

Diziyi kayan pencerelere böler. `window_size` pencere uzunluğu, `step_size` iki pencerenin başlangıçları arasındaki uzaklıktır. Her pencere `s[a:b]` ile aynı sınıftandır, aynı adı, kök alanları ve features kopyasını taşır. Yalnızca tam pencereler döner; dizi `window_size`'dan kısaysa liste boştur. İki değerden biri 1'den küçükse `ValueError` verilir.

```python
d = DNA("ATGCGTA", name="d")
print([p.seq for p in d.sliding(3, step_size=2)])  # -> ['ATG', 'GCG', 'GTA']
print(type(d.sliding(3)[0]).__name__)  # -> DNA
gc = [sum(b in "GC" for b in p) / 4 for p in DNA("GGCCAATT").sliding(4)]
print(gc)  # -> [1.0, 0.75, 0.5, 0.25, 0.0]
```

### Sequence.from_fasta_file

`sınıf metodu: Sequence.from_fasta_file(file_name)`

FASTA dosyasındaki bütün kayıtları bu sınıftan liste olarak okur. `file_name` `.gz` olabilir. Ayrıntı: [`read_fasta`](fileio.md#read_fasta).

```python
DNA("ATGC", name="a").to_fasta("a.fa")
print(DNA.from_fasta_file("a.fa")[0].seq)  # -> ATGC
```

### Sequence.to_fasta

`Sequence.to_fasta(file_name)`

Diziyi `file_name` dosyasına FASTA olarak yazar. FASTA kayıplıdır: yalnızca ad ve dizi geçer. Dizinin adı boş olmamalıdır. Ayrıntı: [`write_fasta`](fileio.md#write_fasta).

### Sequence.to_crn

`Sequence.to_crn(file_name)`

Dizinin tamamını (tip, ad, dizi, kök alanlar, features) `file_name` dosyasına CRN biçiminde yazar.

```python
d = DNA("ATGC", name="d")
d.features.add("kaynak", "elle")
d.to_crn("d.crn")
```

### Sequence.from_crn

`sınıf metodu: Sequence.from_crn(file_name)`

`to_crn` ile yazılmış `file_name` dosyasını okur. `Sequence.from_crn` dosyadaki tipe uygun sınıfı döndürür; `DNA.from_crn` gibi alt sınıftan çağrılırsa dosyadaki tip o sınıf olmalıdır.

```python
geri = Sequence.from_crn("d.crn")
print(type(geri).__name__, geri.features.first("kaynak"))  # -> DNA elle
```

### Python protokolleri

| Ad | Davranış |
|---|---|
| `Sequence.__len__` | Harf sayısı |
| `Sequence.__eq__` | Aynı tip ve aynı dizi ise eşittir; `str` ile büyük/küçük harf gözetmeden karşılaştırır. Ad, features ve strand karşılaştırılmaz |
| `Sequence.__hash__` | `None`: diziler küme ya da sözlük anahtarı olamaz |
| `Sequence.__contains__` | `"ATG" in s` ya da `diger in s` alt dizi arar |
| `Sequence.__iter__` | Harfler üzerinde gezinir |
| `Sequence.__getitem__` | `s[i]` tek harf, `s[a:b]` aynı sınıftan yeni dizi |
| `Sequence.__repr__` | `ad: dizi <uzunluk> [tip]`; 80 harften uzun dizilerde ortası kısaltılır |
| `Sequence._check` | Büyük harfe çevirir ve alfabeyi denetler |
| `Sequence._new` | Aynı sınıf, ad, kök alan ve features kopyasıyla nesne üretir |

```python
d = DNA("ATGCATGC", name="d")
print(len(d), d[0], d[2:5].seq, "GCA" in d, d == "atgcatgc")  # -> 8 A GCA True True
```

## 3. `NucleicAcid`

### NucleicAcid

`NucleicAcid(seq, name='', strand=Strand.FivePrime)`

DNA ve RNA'nın ortak tabanıdır, doğrudan kullanılmaz. `seq`, `name` ve `strand` alır. `complement_table` alt sınıfın tümleme tablosudur.

### NucleicAcid.strand

`NucleicAcid.strand  (özellik)`

Yazılış yönüdür; yalnızca `Strand` atanabilir.

```python
print(DNA("ATG").strand)  # -> Strand.FivePrime
```

### NucleicAcid.seq_5to3

`NucleicAcid.seq_5to3()`

Dizinin 5'→3' okunuşunu boşluklar atılmış metin olarak verir. Yön 3'→5' ise harfleri ters çevirir.

```python
t = DNA("ATGC").complementary()
print(t.seq, t.strand.value[0], t.seq_5to3())  # -> TACG 3' GCAT
```

### NucleicAcid.complement

`NucleicAcid.complement()`

Diziyi yerinde tümleyicisine çevirir, yönü de döndürür ve `self` verir. Tümleyici özgün dizinin altına hizalı yazılır: 5'→3' `ATGC`, 3'→5' `TACG` olur.

```python
d = DNA("ATGC")
d.complement()
print(d.seq, d.strand.value[0])  # -> TACG 3'
```

### NucleicAcid.complementary

`NucleicAcid.complementary()`

`complement` işleminin yeni nesne döndüren biçimidir; özgün dizi değişmez.

```python
print(DNA("ATGC").complementary().seq)  # -> TACG
```

### NucleicAcid.reverse_complement

`NucleicAcid.reverse_complement()`

Karşı zinciri bu diziyle aynı yönde okunacak biçimde verir. 5'→3' `ATGC`, 5'→3' `GCAT` olur. Yeni nesnedir, yön değişmez.

```python
r = DNA("ATGC").reverse_complement()
print(r.seq, r.strand.value[0])  # -> GCAT 5'
```

`NucleicAcid.__repr__` diziye yön bilgisini ekler.

## 4. `DNA`

### DNA

`DNA(seq, name='', strand=Strand.FivePrime)`

`seq`, `name` ve `strand` alır, `DNA_ALPHABET` ile denetler. `molecule`, `alphabet` ve `complement_table` sınıf değişkenleridir. Geçersiz harf `ValueError` verir.

```python
print(DNA("ATGNN-", name="d"))  # -> d: ATGNN- <6> [DNA] 5'->3'
```

### DNA.transcription

`DNA.transcription(template=True, genomic_strand=GenomicStrand.Unknown)`

DNA'dan RNA üretir. `template=True` (varsayılan) bu zincirin kalıp olduğunu söyler: tümleyici RNA üretilir ve yön döner, böylece `rna.reverse_transcription().transcription()` ilk RNA'yı geri verir. `template=False` kodlayan zincir demektir, yalnızca T→U yapılır; NCBI ve GENCODE gen dizileri kodlayan zincirdedir, bu yüzden onlarda `False` kullanılır. `genomic_strand` oluşan RNA'nın genomik zinciridir.

```python
print(DNA("ATGC").transcription(template=False).seq)  # -> AUGC
print(DNA("ATGC").transcription().seq)  # -> UACG
print(DNA("ATG").transcription(False, genomic_strand=GenomicStrand.Minus).genomic_strand)  # -> GenomicStrand.Minus
```

### DNA.find_orfs

`DNA.find_orfs(minimum_length_of_orf=-1, genetic_code=GeneticCode.STANDARD, alternative_starts=False)`

İki zincirdeki bütün açık okuma çerçevelerini bulur. ORF, başlangıç kodonundan aynı çerçevedeki ilk stopa (stop dahil) kadardır. Bu zincir ve tümleyicisi toplam altı çerçevede taranır, stop'suz kısmi ORF alınmaz, aynı stop'u paylaşan iç içe ORF'lerin hepsi döner. `minimum_length_of_orf` stop dahil en küçük uzunluktur (nt), `-1` sınırsızdır. `genetic_code` başlangıç ve stop kodonlarını belirler. `alternative_starts=True` kodun bütün başlangıç kodonlarını, `False` yalnızca ATG'yi kullanır.

Sonuç konuma göre sıralı, kodlayan yönde 5'→3' DNA listesidir. Adı `<ad>_ORF<n>` olur; features'a `orf_start`, `orf_end` (bu dizide 1 tabanlı), `orf_strand` ve `orf_frame` yazılır.

```python
d = DNA("CCATGAAATGATAACC", name="d")
for orf in d.find_orfs():
    print(orf.name, orf.seq, orf.features.first("orf_start"), orf.features.first("orf_frame"))
print(len(d.find_orfs(minimum_length_of_orf=12)))  # -> 0
print(len(d.find_orfs(genetic_code=GeneticCode.STANDARD, alternative_starts=True)))  # -> 1
```

## 5. `RNA`

### RNA

`RNA(seq, name='', strand=Strand.FivePrime, genomic_strand=GenomicStrand.Unknown)`

`seq`, `name`, `strand` ve `genomic_strand` alır. `root_types` CRN okuyucusuna hangi kök alanların bulunduğunu bildirir. Alfabede T yoktur.

```python
r = RNA("AUGGCCUAA", name="r", genomic_strand=GenomicStrand.Plus)
print(r)  # -> r: AUGGCCUAA <9> [RNA] 5'->3' (+)
```

### RNA.genomic_strand

`RNA.genomic_strand  (özellik)`

RNA'nın transkript edildiği genomik zincirdir. `complement` bu değeri değiştirmez, yalnızca `GenomicStrand` atanabilir.

```python
print(r.genomic_strand.value[0])  # -> +
```

### RNA.reverse_transcription

`RNA.reverse_transcription()`

Tümleyici, antiparalel ilk zincir cDNA üretir: A→T, U→A, G→C, C→G; yön döner. mRNA ile aynı anlamlı ikinci zincir için sonucun `complementary()`'si kullanılır. `genomic_strand` anahtarı features'tan çıkarılır.

```python
c = RNA("AUGC").reverse_transcription()
print(c.seq, c.strand.value[0], c.complementary().seq)  # -> TACG 3' ATGC
```

### RNA.translation

`RNA.translation(genetic_code=GeneticCode.STANDARD)`

İlk nükleotitten son tam kodona kadar üçerli okur. Başlangıç ya da stop aranmaz, stoplar `*` yazılır ve okuma sürer. 3'→5' yazılmış RNA önce çevrilir, boşluklar atılır, belirsiz baz içeren kodon `X` olur. `genetic_code` kullanılacak kodu seçer.

```python
print(RNA("AUGGCCUAAGG").translation().seq)  # -> MA*
print(RNA("AUAUGA").translation(genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL).seq)  # -> MW
```

### RNA.match_orf

`RNA.match_orf(protein, genetic_code=GeneticCode.STANDARD)`

Çevirisi `protein` ile birebir aynı olan ORF'yi (AUG ... stop) bulur. Yalnızca RNA'nın kendi zincirine bakılır, sondaki stop simgesi iki tarafta da yok sayılır, birden çok ORF uyarsa en baştaki seçilir. `genetic_code` çeviri kodudur. Sonuç stop dahil 5'→3' RNA'dır, uyan ORF yoksa `None` döner. Yeni RNA'ya `orf_start`, `orf_end` ve `cds = (1, uzunluk)` yazılır, `utr5` ve `utr3` silinir. `Gene.isoforms` ve `CUB` bu metoda dayanır.

```python
r = RNA("GGAUGGCCUAAGG", name="r")
orf = r.match_orf("MA")
print(orf.seq, orf.features.first("orf_start"), orf.features.first("orf_end"))  # -> AUGGCCUAA 3 11
print(r.match_orf(PRT("MK")))  # -> None
```

`RNA.__repr__` diziye genomik zinciri ekler.

## 6. `PRT`

### PRT

`PRT(seq, name='')`

Tek harfli kodlarla N-ucundan C-ucuna yazılmış proteindir. `seq` ve `name` alır, `PRT_ALPHABET` ile denetler. `PRT.__repr__` çıktıya `N->C` ekler.

```python
p = RNA("AUGGCCUAA").translation()
print(p)  # -> : MA* <3> [PRT] N->C
print(PRT("MKV", name="p").seq)  # -> MKV
```
