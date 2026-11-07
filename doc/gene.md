# `crn.gene`: `Gene` kapsayıcısı

Kaynak: [`crn/gene.py`](../crn/gene.py)

```text
Gene (ör. MEFV)
 ├── dna          : DNA            genomik dizi (tek)
 ├── transcripts  : {id: RNA}      transkript varyantları
 ├── proteins     : {id: PRT}      proteinler
 └── features     : Features       gene_id, locus ...
```

Transkriptler ve proteinler bağımsız eklenir, aralarındaki eşleşme `isoforms()` ile kurulur. `Ncbi` ve `Gencode` her proteinin `transcript_id` özelliğine kaynağın bildirdiği transkripti yazar.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`Gene`](#gene) · [`dna`](#genedna) · [`features`](#genefeatures) · [`add_transcript`](#geneadd_transcript) · [`add_protein`](#geneadd_protein) · [`transcript`](#genetranscript) · [`protein`](#geneprotein) · [`remove_transcript`](#generemove_transcript) · [`remove_protein`](#generemove_protein) · [`transcripts`](#genetranscripts) · [`proteins`](#geneproteins) · [`isoforms`](#geneisoforms) · [`to_crn`](#geneto_crn) · [`from_crn`](#genefrom_crn)
- [2. İç yardımcılar](#2-iç-yardımcılar)
<!-- /toc -->

## 1. Sınıf

### Gene

`Gene(name, dna=None)`

`name` genin adı, `dna` isteğe bağlı genomik dizidir.

```python
from crn import DNA, RNA, PRT, Gene

gen = Gene("ornek", dna=DNA("CCATGGCCTAACC", name="ornek"))
print(gen)  # -> ornek: DNA <13>, 0 transkript, 0 protein
```

### Gene.dna

`Gene.dna  (özellik)`

Genomik DNA'dır, yalnızca `DNA` atanabilir.

```python
print(gen.dna.seq)  # -> CCATGGCCTAACC
```

### Gene.features

`Gene.features  (özellik)`

Genin `Features` deposudur. `gene_name`, `gene_id`, `group` gibi bilgiler burada durur.

```python
gen.features.set("group", "kontrol")
print(gen.features.first("group"))  # -> kontrol
```

### Gene.add_transcript

`Gene.add_transcript(rna, transcript_id=None)`

`rna` nesnesini transkript olarak ekler ve döndürür. `transcript_id` verilmezse RNA'nın adı kimlik olur. `rna` bir `RNA` değilse `TypeError`, kimlik boşsa ya da zaten varsa `ValueError` verilir.

```python
gen.add_transcript(RNA("GGAUGGCCUAAGG", name="t1"))
gen.add_transcript(RNA("AUGGCCUAA", name="ikinci"), transcript_id="t2")
print(list(gen.transcripts))  # -> ['t1', 't2']
```

### Gene.add_protein

`Gene.add_protein(prt, protein_id=None)`

`prt` nesnesini ekler ve döndürür. `protein_id` verilmezse proteinin adı kimlik olur; hatalar `add_transcript` ile aynıdır.

```python
p = gen.add_protein(PRT("MA", name="p1"))
p.features.set("transcript_id", "t1")
print(list(gen.proteins))  # -> ['p1']
```

### Gene.transcript

`Gene.transcript(transcript_id)`

`transcript_id` kimlikli RNA'yı verir, yoksa `KeyError`.

```python
print(gen.transcript("t2").seq)  # -> AUGGCCUAA
```

### Gene.protein

`Gene.protein(protein_id)`

`protein_id` kimlikli proteini verir, yoksa `KeyError`.

```python
print(gen.protein("p1").seq)  # -> MA
```

### Gene.remove_transcript

`Gene.remove_transcript(transcript_id)`

`transcript_id` kimlikli transkripti siler; bulunduysa `True` döner.

```python
print(gen.remove_transcript("t2"), gen.remove_transcript("t2"))  # -> True False
```

### Gene.remove_protein

`Gene.remove_protein(protein_id)`

`protein_id` kimlikli proteini siler; bulunduysa `True` döner.

```python
print(gen.remove_protein("yok"))  # -> False
```

### Gene.transcripts

`Gene.transcripts  (özellik)`

`{kimlik: RNA}` sözlüğünü ekleme sırasıyla verir. Sözlük kopyadır, nesneler paylaşılır.

### Gene.proteins

`Gene.proteins  (özellik)`

`{kimlik: PRT}` sözlüğünü ekleme sırasıyla verir. Sözlük kopyadır.

### Gene.isoforms

`Gene.isoforms(genetic_code=GeneticCode.STANDARD, unique_cds=True, priority=TRANSCRIPT_PRIORITY)`

Geni izoformlarına ayırır: her biri bir DNA, bir transkript ve bir protein içeren genlerin listesi. `CUB` tam olarak bu biçimi ister. Her transkript için çevirisi bir proteinle birebir aynı ORF aranır; önce kaynağın o transkripte bağladığı protein, sonra diğerleri denenir. Eşleşen her transkript proteiniyle yeni bir gene konur, özgün gen değişmez.

Yeni genin adı `<gen>_<transkript>` olur. features genin kopyasıdır, ayrıca `isoform_of`, `transcript_id` ve `protein_id` yazılır. Transkript kırpılmaz; ORF konumu `orf_start` ve `orf_end` özelliklerindedir. `genetic_code` çeviri kodudur.

`unique_cds=True` aynı CDS'i veren izoformlardan yalnızca birini tutar. Aynı CDS'in iki kez sayılması kodon sayımlarını ve referans setlerini çarpıtacağından bu ayar açık kalmalıdır. Kimliği `priority` içinde önce gelen önekle başlayan kayıt seçilir (önce `NM_`, sonra `XM_`), eşitlikte ilk eklenen kalır. Çıkarılanlar ve eşleşmeyen kayıtlar `UserWarning` ile bildirilir. Genomik DNA yoksa `ValueError` verilir.

```python
import warnings
from crn import GeneticCode

gen = Gene("ornek", dna=DNA("CCATGGCCTAACC"))
gen.add_transcript(RNA("GGAUGGCCUAAGG", name="NM_1"))
gen.add_transcript(RNA("AUGGCCUAA", name="XM_1"))
gen.add_protein(PRT("MA", name="p1"))
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    bir = gen.isoforms()
    iki = gen.isoforms(genetic_code=GeneticCode.STANDARD, unique_cds=False)
    ters = gen.isoforms(priority=("XM_", "NM_"))
print([g.name for g in bir], len(iki), ters[0].name)  # -> ['ornek_NM_1'] 2 ornek_XM_1
print(bir[0].features.first("isoform_of"), bir[0].features.first("protein_id"))  # -> ornek p1
```

### Gene.to_crn

`Gene.to_crn(file_name)`

Geni bütün içeriğiyle `file_name` dosyasına yazar.

```python
gen.to_crn("ornek.crn")
```

### Gene.from_crn

`sınıf metodu: Gene.from_crn(file_name)`

`to_crn` ile yazılmış `file_name` dosyasını okur.

```python
print(Gene.from_crn("ornek.crn"))  # -> ornek: DNA <13>, 2 transkript, 1 protein
```

## 2. İç yardımcılar

| Ad | Görevi |
|---|---|
| `Gene._add` | `add_transcript` ve `add_protein`'in ortak denetimi: tip, boş kimlik, yinelenen kimlik |
| `Gene._unique_cds` | Aynı CDS'i veren izoformlardan önceliği en yüksek olanı tutar, sırayı korur ve çıkarılanları uyarır |
| `Gene.__repr__` | `TP53: DNA <25768>, 1 transkript, 1 protein` biçiminde özet |
