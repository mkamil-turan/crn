# `crn.fetch.gencode`: GENCODE

Kaynak: [`crn/fetch/gencode.py`](../../crn/fetch/gencode.py)

Bir GENCODE sürümünü bütünüyle indirir, açar ve dizinler; sonra herhangi bir geni dosyadan bayt konumuyla okur. İlk kurulum uzun sürer (v50 için 1,3 GB indirme, açılmış halde 5 GB), sonraki her sorgu ağsız ve hızlıdır. Bu yüzden genom ölçekli çalışmalarda birincil kaynak budur.

```text
sys/gencode/<tür>/<sürüm>/          ör. sys/gencode/human/v50/
    raw/                            GENCODE'dan indirilen .gz dosyaları
    genome.fa, genome.fai           açılmış genom ve konum dizini
    transcripts.fa, transcripts.idx transkript dizileri ve dizini
    translations.fa, translations.idx
                                    protein dizileri ve dizini
    genes.tsv, transcripts.tsv      GTF'den çıkarılan koordinatlar
    regions.tsv                     protein kodlayan transkriptlerin
                                    UTR5/CDS/UTR3 koordinatları
    manifest.json                   hazırlık tamamlandığında yazılır
    details.tsv                     etiketler, destek düzeyi ve ekzonlar
                                    (ilk gerektiğinde GTF'den üretilir)
```

<!-- toc -->
- [1. Sabitler ve hata](#1-sabitler-ve-hata): [`AmbiguousGeneError`](#ambiguousgeneerror)
- [2. Sürüm ve kurulum](#2-sürüm-ve-kurulum): [`Gencode`](#gencode) · [`species_dir`](#gencodespecies_dir) · [`release`](#gencoderelease) · [`release_dir`](#gencoderelease_dir) · [`installed_releases`](#gencodeinstalled_releases) · [`is_ready`](#gencodeis_ready) · [`latest_release`](#gencodelatest_release) · [`update`](#gencodeupdate)
- [3. Sorgu: gen ve transkript listeleri](#3-sorgu-gen-ve-transkript-listeleri): [`gene_types`](#gencodegene_types) · [`transcript_types`](#gencodetranscript_types) · [`genes`](#gencodegenes) · [`transcripts`](#gencodetranscripts) · [`details`](#gencodedetails)
- [4. Diziler](#4-diziler): [`gene`](#gencodegene) · [`canonical`](#gencodecanonical) · [`dna`](#gencodedna) · [`rna`](#gencoderna) · [`protein`](#gencodeprotein) · [`region`](#gencoderegion)
- [5. İç yardımcılar](#5-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler ve hata

`GENCODE_URL` sunucu adresidir. `SPECIES` desteklenen türleri (`human`, `mouse`) sunucu klasörü ve sürüm önekiyle (`v`, `vM`) eşler. `Gencode.SOURCE` features'a yazılan kaynak adıdır: `GENCODE`.

### AmbiguousGeneError

Aynı adı taşıyan birden çok gen olduğunda verilir, `ValueError` alt sınıfıdır. Örneğin `PDE8B` ve `SHOX` adları birden çok kayda uyar; bu durumda `gene_id` kullanılmalıdır. Mesajda aday kimlikler yazılır.

```python
from crn.fetch import Gencode
from crn.fetch.gencode import AmbiguousGeneError

gc = Gencode(verbose=False)
try:
    gc.gene("PDE8B")
except AmbiguousGeneError as hata:
    mesaj = str(hata)
print(mesaj[:34])  # -> 'PDE8B' adında birden çok gen var:
```

## 2. Sürüm ve kurulum

### Gencode

`Gencode(species='human', release=None, sys_dir=None, base_url=GENCODE_URL, verbose=True)`

`species` `"human"` ya da `"mouse"`, başka değer `ValueError` verir. `release` `"v50"`, `"50"` ya da fare için `"vM39"` gibi bir sürümdür; `None` ise `sys/` içindeki en yeni hazır sürüm kullanılır, hiç yoksa `update()` sunucudaki en son sürümü indirir. `sys_dir` veri klasörü, `base_url` sunucu adresi, `verbose` ilerleme çubuklarıdır.

```python
print(Gencode(species="human", release="50", verbose=False).release)  # -> v50
```

### Gencode.species_dir

`Gencode.species_dir  (özellik)`

Türün klasörü: `sys/gencode/human`.

### Gencode.release

`Gencode.release  (özellik)`

Kullanılan sürüm. Belirtilmediyse `sys/` içindeki en yeni hazır sürümdür.

```python
print(gc.release)  # -> v50
```

### Gencode.release_dir

`Gencode.release_dir  (özellik)`

Sürümün klasörü. Sürüm belirli değilse `ValueError` verir.

```python
print(gc.release_dir.name)  # -> v50
```

### Gencode.installed_releases

`Gencode.installed_releases()`

`sys/` içinde hazırlığı tamamlanmış sürümleri eskiden yeniye verir.

```python
print(gc.installed_releases())  # -> ['v50']
```

### Gencode.is_ready

`Gencode.is_ready()`

Kullanılan sürümün bütün hazırlık dosyaları (`_OUTPUTS`) yerinde mi söyler.

```python
print(gc.is_ready())  # -> True
```

### Gencode.latest_release

`Gencode.latest_release()`

Sunucudaki en son sürümü verir, ağ gerekir.

<!-- calistirma -->
```python
print(gc.latest_release())
```

### Gencode.update

`Gencode.update(force=False)`

Seçili sürümü indirir ve hazırlar. Nesne bir sürümle oluşturulduysa o sürüm, oluşturulmadıysa sunucudaki en son sürüm alınır (`sys/` içinde daha eski sürüm olsa bile). Her dosya MD5SUMS ile doğrulanır. Doğrulanmış dosyalar yeniden indirilmez, hazırlanmış dosyalar yeniden üretilmez, yalnızca eksikler tamamlanır. Bitince `manifest.json` yazılır. `force=True` hepsini yeniler. Hazırlanan sürümü döndürür. Veri yoksa ilk sorgu bunu kendisi çağırır.

<!-- calistirma -->
```python
Gencode(release="v50").update()
Gencode(species="mouse", release="vM39").update(force=True)
```

## 3. Sorgu: gen ve transkript listeleri

Bu metotlar dizi değil tablo satırı döndürür, yani gen seti kurmanın başlangıç noktasıdır.

### Gencode.gene_types

`Gencode.gene_types()`

`{gene_type: gen sayısı}` sözlüğünü çoktan aza verir.

```python
print(list(gc.gene_types().items())[:2])  # -> [('lncRNA', 35020), ('protein_coding', 20131)]
```

### Gencode.transcript_types

`Gencode.transcript_types()`

`{transcript_type: transkript sayısı}` sözlüğünü çoktan aza verir.

```python
print(list(gc.transcript_types().items())[0])  # -> ('protein_coding', 278484)
```

### Gencode.genes

`Gencode.genes(gene_type=None, chromosome=None, start=None, end=None, name=None)`

Koşullara uyan gen satırlarını GTF sırasıyla verir. Satır `gene_id, gene_name, gene_type, chromosome, start, end, strand` içerir, `start` ve `end` sayıdır. Verilmeyen koşul uygulanmaz. `gene_type` tip ya da tip listesi, `chromosome` kromozom adıdır. `start` ve `end` bu aralıkla örtüşen genleri seçer ve `chromosome` gerektirir, yoksa `ValueError` verilir. `name` gen adında aranacak düzenli ifadedir.

```python
hox = gc.genes(gene_type="protein_coding", name=r"^HOXA\d+$")
print(len(hox), hox[0]["gene_name"], hox[0]["chromosome"])  # -> 11 HOXA1 chr7
print(len(gc.genes(gene_type=["miRNA"], chromosome="chr21")))  # -> 29
print([g["gene_name"] for g in gc.genes(chromosome="chr17", start=7668402, end=7687550, gene_type="protein_coding")])  # -> ['TP53', 'WRAP53']
```

### Gencode.transcripts

`Gencode.transcripts(gene=None, transcript_type=None, chromosome=None, start=None, end=None, tag=None)`

Koşullara uyan transkript satırlarını verir: `transcript_id, gene_id, transcript_name, transcript_type, chromosome, start, end, strand`. `gene` gen adı ya da `gene_id`'dir, `transcript_type` tip ya da tiplerdir. `chromosome`, `start` ve `end` `genes` ile aynıdır. `tag` bir GTF etiketidir (`MANE_Select`, `basic`) ve `details.tsv` gerektirir, yoksa ilk kullanımda üretilir. Gen bulunamazsa `KeyError`, ad belirsizse `ValueError` verir.

```python
print(gc.transcripts(gene="TP53", tag="MANE_Select")[0]["transcript_name"])  # -> TP53-201
print(len(gc.transcripts(gene="TP53", transcript_type="protein_coding")) > 5)  # -> True
```

### Gencode.details

`Gencode.details(transcript)`

Transkriptin GTF ayrıntılarını verir: `level` GENCODE kanıt düzeyidir (1 doğrulanmış, 2 elle, 3 otomatik), `support_level` TSL değeridir (`"1"` en güçlü, `""` yok), `tags` etiketlerdir, `exons` genomik 1 tabanlı koordinatlardır ve transkriptin 5'→3' sırasıyla gelir: (-) zincirde koordinatlar azalır. `transcript` kimlik ya da addır, yoksa `KeyError`.

```python
d = gc.details("ENST00000269305")
print(d["level"], d["support_level"], "MANE_Select" in d["tags"], d["exons"][0])  # -> 2 1 True (7687377, 7687490)
```

## 4. Diziler

### Gencode.gene

`Gencode.gene(gene)`

Geni bütün içeriğiyle verir: genomik DNA, dizisi olan bütün transkriptler (kimlik = `transcript_id`) ve proteinler (kimlik = `protein_id`). `gene` gen adı ya da `gene_id`'dir.

```python
tp53 = gc.gene("TP53")
print(tp53.name, tp53.features.first("gene_type"), len(tp53.proteins) > 5)  # -> TP53 protein_coding True
```

### Gencode.canonical

`Gencode.canonical(gene)`

Geni yalnızca ana transkripti ve proteiniyle verir. Ana transkript `MANE_Select`, yoksa `Ensembl_canonical` etiketli protein kodlayan transkripttir; SHANK3 gibi genlerde ve mitokondride MANE kaydı yoktur. Hangisinin seçildiği `canonical_tag` özelliğindedir. Her gen için bir temsilci seçen analizlerde bu kullanılır, çünkü izoformlar ayrı sayılırsa çok izoformlu genler sonucu domine eder. `gene` bulunamazsa `KeyError`, bu etiketlerde protein kodlayan transkript yoksa `ValueError` verir.

```python
ana = gc.canonical("TP53")
print(ana, ana.features.first("canonical_tag"))  # -> TP53: DNA <25768>, 1 transkript, 1 protein MANE_Select
print(gc.canonical("SHANK3").features.first("canonical_tag"))  # -> Ensembl_canonical
```

### Gencode.dna

`Gencode.dna(gene)`

Genin (`gene`) genomik DNA'sını kodlayan zincirde, 5'→3' verir. (-) zincirdeki genlerde referans bölgenin ters tümleyicisi alınır.

```python
print(len(gc.dna("TP53")))  # -> 25768
```

### Gencode.rna

`Gencode.rna(transcript)`

Transkriptin olgun, birleştirilmiş dizisini `RNA` olarak verir. `transcript` sürümlü ya da sürümsüz `transcript_id` ya da transkript adıdır (`TP53-201`). Protein kodlayan transkriptlerde `utr5`, `cds` ve `utr3` features'tadır.

```python
rna = gc.rna("ENST00000269305")
bas, son = rna.features.first("cds")
print(rna[bas - 1:son].translation().seq[:10], rna.genomic_strand.value[0])  # -> MEEPQSDPSV -
```

### Gencode.protein

`Gencode.protein(identifier)`

Proteini verir. `identifier` ENSP, ENST ya da transkript adıdır, sürümsüz olabilir. Kodlamayan transkriptte `KeyError`.

```python
print(gc.protein("ENST00000269305").seq[:10])  # -> MEEPQSDPSV
```

### Gencode.region

`Gencode.region(chromosome, start, end, strand='+', name='')`

Genomdan bölge çıkarır. `chromosome` kromozom adı, `start` ve `end` 1 tabanlı ve iki ucu dahildir (GTF ile aynı). `strand="-"` ters tümleyiciyi alır, iki durumda da DNA 5'→3' yazılır. `name` boşsa `chr16:3242028-3256627(-)` biçiminde konur. Kromozom yoksa `KeyError`, konum dışarıdaysa ya da zincir geçersizse `ValueError` verir.

```python
print(gc.region("chr17", 7668402, 7668411, "-", name="p53_ucu").seq, gc.region("chr17", 7668402, 7668411).name)  # -> CTGAGGGGTG chr17:7668402-7668411(+)
```

## 5. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_ATTRIBUTE` | GTF öznitelik sütununun düzenli ifadesi |
| `_VERSION` | Kimliklerdeki `.13` gibi sürüm eki |
| `_GENES_COLUMNS` | `genes.tsv` sütunları |
| `_TRANSCRIPTS_COLUMNS` | `transcripts.tsv` sütunları |
| `_REGION` | pc_transcripts başlığındaki `UTR5:1-367` kalıbı |
| `_LEVEL` | GTF `level` değeri |
| `_TAG` | GTF `tag` değeri |
| `_OUTPUTS` | Hazırlığın ürettiği dosyalar; `is_ready` bunlara bakar |
| `_unversioned` | `ENSG00000103313.13` → `ENSG00000103313`; `_PAR_Y` gibi ekleri korur |
| `_reverse_complement` | Metin üzerinde ters tümleyici; nesne kurmadan genom okumak için |
| `Gencode._normalize_release` | `"50"` → `"v50"`, farede `"vM39"` |
| `Gencode._remote_dir` | Sürümün sunucu klasörü |
| `Gencode._prepare` | Sürümü indirir, eksik hazırlık dosyalarını üretir, `manifest.json` yazar |
| `Gencode._ensure` | Veri hazır değilse `update` çağırır; sürüm belirtilmedi ve yerel sürüm varsa ağa çıkmaz |
| `Gencode._select_files` | MD5SUMS listesinden gerekli beş dosyayı seçer |
| `Gencode._read_checksums` | MD5SUMS dosyasını okur |
| `Gencode._read_url` | URL'yi belleğe okur |
| `Gencode._download` | `.part` dosyasına indirir, MD5 doğrular, yerine koyar; doğru dosyayı atlar |
| `Gencode._md5` | Dosyanın MD5'ini ilerleme çubuğuyla hesaplar |
| `Gencode._gzip_lines` | gzip satırlarını verir; çubuk sıkıştırılmış bayta göre ilerler |
| `Gencode._file_lines` | Metin satırlarını okunan bayta göre çubukla verir |
| `Gencode._build_genome` | Genomu açar ve samtools faidx biçiminde `genome.fai` yazar |
| `Gencode._build_fasta` | FASTA'yı açar ve kimlik, konum, uzunluk, başlık dizini yazar |
| `Gencode._build_tables` | GTF'den `genes.tsv` ve `transcripts.tsv` üretir |
| `Gencode._build_regions` | pc_transcripts başlıklarından `regions.tsv` üretir |
| `Gencode._tmp` | Yolun `.tmp` karşılığı |
| `Gencode._commit` | `.tmp` dosyalarını gerçek adlarına taşır |
| `Gencode._load` | Tabloları ve dizinleri belleğe bir kez yükler |
| `Gencode._find_gene` | Gen adı ya da kimliğinden satır bulur; `_PAR_Y` kopyalarını eler, belirsizse `AmbiguousGeneError` |
| `Gencode._find_transcript` | Transkript satırını bulur |
| `Gencode._read_record` | Dizinden bir kaydı bayt konumuyla okur |
| `Gencode._features` | Features'a `release` değerini de ekler |
| `Gencode._select` | `genes` ve `transcripts`'in ortak filtre döngüsü |
| `Gencode._load_details` | `details.tsv`'yi bir kez yükler, yoksa üretir |
| `Gencode._build_details` | Ham GTF'den `details.tsv` yazar |
| `Gencode._gene_object` | Gen satırından, verilen transkriptleri içeren `Gene` kurar |
