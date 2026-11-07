# `crn.fetch.ncbi`: NCBI RefSeq

Kaynak: [`crn/fetch/ncbi.py`](../../crn/fetch/ncbi.py)

Bir geni NCBI'dan tüm içeriğiyle çeker: genomik DNA, gene bağlı RefSeq RNA'ları (NM_, NR_, XM_, XR_) ve bunların proteinleri. Her kayıt `sys/ncbi/` içinde saklanır, sonraki isteklerde ağa çıkılmaz. Bu sayfadaki `MEFV`, `TP53` gibi örnekler `sys/` içinde kayıt varsa çevrimdışı çalışır.

```text
sys/ncbi/
    symbols/<organizma>.json      gen sembolü -> NCBI Gene ID
    genes/<gene_id>.json          genin özeti ve RNA/protein listesi
                                  (update tamamlandığında yazılır)
    nuccore/<erişim>.fa, .ft      transkript dizisi ve özellik tablosu
    protein/<erişim>.fa           protein dizisi
    region/<erişim>_<baş>_<son>_<zincir>.fa
                                  genomik bölge
```

<!-- toc -->
- [1. Sabitler](#1-sabitler)
- [2. Sınıf](#2-sınıf): [`Ncbi`](#ncbi) · [`ncbi_dir`](#ncbincbi_dir) · [`search_gene`](#ncbisearch_gene) · [`update`](#ncbiupdate) · [`gene`](#ncbigene) · [`dna`](#ncbidna) · [`rna`](#ncbirna) · [`protein`](#ncbiprotein) · [`region`](#ncbiregion)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler

`NCBI_URL` E-utilities kök adresidir. `ORGANISMS` desteklenen organizmaları NCBI taxon kimliklerine bağlar: `human`, `mouse`, `rat`, `zebrafish`, `fruit_fly`, `worm`.

## 2. Sınıf

### Ncbi

`Ncbi(organism='human', session=None, sys_dir=None, verbose=True)`

`organism` gen sembolü aramasında kullanılacak organizmadır, `ORGANISMS` dışındaki değer `ValueError` verir. `session` verilmezse API anahtarsız bir `NcbiSession` açılır. `sys_dir` veri klasörü, `verbose` ilerleme çubuklarıdır.

```python
from crn.fetch import Ncbi
from crn import GenomicStrand

ncbi = Ncbi(verbose=False)
print(ncbi.organism, ncbi.SOURCE)  # -> human NCBI
```

### Ncbi.ncbi_dir

`Ncbi.ncbi_dir  (özellik)`

Verilerin bulunduğu klasör: `sys/ncbi`.

### Ncbi.search_gene

`Ncbi.search_gene(gene)`

Gen sembolünü (`gene`) NCBI Gene ID'sine çevirir. Rakamlardan oluşan değer zaten Gene ID kabul edilir. NCBI araması takma adları da eşleştirdiği için birden çok sonuç varsa resmi sembolü tam tutan gen seçilir (`MIP` için TNPO1 değil MIP), yoksa ilk sonuç alınır. Bulunan eşleşme saklanır. Sembol bulunamazsa `KeyError` verir.

```python
print(ncbi.search_gene("MEFV"), ncbi.search_gene("4210"))  # -> 4210 4210
```

### Ncbi.update

`Ncbi.update(gene, force=False)`

Genin bütün kayıtlarını indirip `sys/` içine yazar: özet, genomik bölge, RefSeq RNA'ları, özellik tabloları ve proteinler. `gene` sembol ya da Gene ID'dir. Var olan kayıtlar yeniden indirilmez, `force=True` hepsini yeniler. En sonda `genes/<gene_id>.json` yazılır; bu dosya yoksa gen hazır sayılmaz, yani kesilen indirme sonraki çağrıda tamamlanır. Gene ID döndürür. `dna`, `gene` gibi metotlar gerektiğinde bunu kendileri çağırır, elle çağırmak isteğe bağlıdır.

```python
print(ncbi.update("MEFV"))  # -> 4210
```

### Ncbi.gene

`Ncbi.gene(gene)`

Geni bütün içeriğiyle `Gene` olarak döndürür: genomik DNA (koordinat yoksa boş), bütün RefSeq RNA'lar (kimlik = erişim) ve proteinler. Proteinin hangi transkriptten geldiği `transcript_id` özelliğindedir. `gene` sembol ya da Gene ID'dir. `CUB`'a vermeden önce `isoforms()` ile ayrılmalıdır.

```python
mefv = ncbi.gene("MEFV")
print(mefv.name, len(mefv.transcripts) > 0, mefv.features.first("source"))  # -> MEFV True NCBI
```

### Ncbi.dna

`Ncbi.dna(gene)`

Genin (`gene`) genomik DNA'sını kodlayan zincirde, 5'→3' verir. Gen ya da koordinat bulunamazsa `KeyError`.

```python
print(ncbi.dna("MEFV")[:12].seq)  # -> CTACCAGAAGCC
```

### Ncbi.rna

`Ncbi.rna(accession, genomic_strand=GenomicStrand.Unknown)`

RefSeq transkriptini `RNA` olarak verir. `accession` sürümlü ya da sürümsüz olabilir. Protein kodlayan transkriptlerde features'a `utr5`, `cds` ve `utr3` yazılır (transkript üzerinde 1 tabanlı, iki ucu dahil). `genomic_strand` RNA'nın genomik zinciridir; `gene()` bunu kendisi doldurur. Kayıt yoksa `KeyError`.

```python
rna = ncbi.rna("NM_000243", genomic_strand=GenomicStrand.Minus)
bas, son = rna.features.first("cds")
print(rna.name, rna[bas - 1:son].translation().seq[:10])  # -> NM_000243.3 MAKTPSDHLL
```

### Ncbi.protein

`Ncbi.protein(accession)`

Proteini `PRT` olarak verir. `accession` protein erişimi (NP_, XP_) ya da onu kodlayan transkriptin erişimidir (NM_, XM_). Kodlamayan RNA'da `KeyError`.

```python
print(ncbi.protein("NM_000243").seq[:10])  # -> MAKTPSDHLL
```

### Ncbi.region

`Ncbi.region(accession, start, end, strand='+', name='')`

Bir genomik bölgeyi `DNA` olarak verir. `accession` kromozom ya da contig erişimidir (`NC_000016.10`), `start` ve `end` 1 tabanlı ve iki ucu dahildir. `strand="-"` ters tümleyiciyi alır, iki durumda da DNA 5'→3' yazılır. `name` boşsa `NC_000016.10:3242028-3256627(-)` biçiminde ad konur. Geçersiz konum ya da zincir `ValueError` verir.

<!-- calistirma -->
```python
bolge = ncbi.region("NC_000016.10", 3242028, 3242100, "-", name="MEFV_3")
print(len(bolge))
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_BATCH` | Tek `efetch` isteğindeki en çok kayıt: 100 |
| `_CDS` | Özellik tablosundaki CDS satırının düzenli ifadesi |
| `_QUALIFIER` | Özellik tablosundaki niteleyici satırının düzenli ifadesi |
| `_split_feature_table` | Çok kayıtlı özellik tablosunu (erişim, metin) çiftlerine böler |
| `_parse_feature_table` | Gen adını, GeneID'yi, ilk CDS'i ve proteinini çıkarır; eksik uç işaretlerini (`<`, `>`) atar |
| `Ncbi._gene_file` | `genes/<gene_id>.json` yolu |
| `Ncbi._symbols_file` | Organizmanın sembol önbelleği yolu |
| `Ncbi._record_file` | `db/erişim.uzantı` yolu |
| `Ncbi._cached_accession` | `sys/` içindeki tam erişimi bulur; sürümsüz verilirse en yeni sürümü seçer |
| `Ncbi._get` | E-utilities uç noktasına oturumla istek atar |
| `Ncbi._efetch` | Kayıtları `_BATCH`'lik gruplarla indirip metinleri birleştirir |
| `Ncbi._summary` | `esummary` ile gen özetini alır; bulunamazsa `KeyError` |
| `Ncbi._elink` | `elink` ile gene bağlı kayıt kimliklerini alır |
| `Ncbi._download_rnas` | RNA FASTA'larını ve eksik özellik tablolarını indirir. NCBI UID ile çalıştığı için hangi kaydın yerelde olduğu indirmeden bilinemez, bu yüzden FASTA her seferinde iner |
| `Ncbi._download_proteins` | Eksik protein dizilerini indirir |
| `Ncbi._region_file` | Bölge dosyasının yolu |
| `Ncbi._download_region` | Bölgeyi indirir ve yolunu verir |
| `Ncbi._gene_region` | Özetteki koordinatları 1 tabanlıya çevirir; `chrstart > chrstop` ise gen (-) zincirdedir |
| `Ncbi._gene_record` | Genin yerel kaydını verir, yoksa `update` çağırır |
| `Ncbi._features` | `LocalSource._features`'ı boş liste de boş sayılacak biçimde çağırır |
