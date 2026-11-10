# `crn.fetch.hgnc`: HGNC gen tablosu ve ortologlar

Kaynak: [`crn/fetch/hgnc.py`](../../crn/fetch/hgnc.py)

HGNC insan genlerinin onaylı sembollerini ve gen ailelerini, HCOP ise bu genlerin 19 türdeki ortologlarını verir. Bir aileyi ya da yolağı gen listesine çevirmek ve hasta genlerinin fare, zebra balığı ya da maya karşılıklarını bulmak için kullanılır.

HGNC yalnızca insan genlerini adlandırır; başka bir organizmanın "HGNC kaydı" yoktur. Başka bir türle bağ HCOP ortologları üzerinden kurulur ve tür HCOP adıyla (`species="chicken"`) ya da NCBI Taxonomy ID'siyle (`taxid=9031`) seçilir. HCOP'un kapsamadığı bir türün genleri `Ncbi(taxid=...)` ile alınır; ortolog bulunan genin `ortholog_entrez_gene` değeri doğrudan `Ncbi.gene` ile kullanılabilir.

```text
sys/hgnc/
    hgnc_complete_set.txt           HGNC'nin tüm genleri (sekmeyle ayrılmış)
    raw/human_<tür>_hcop_fifteen_column.txt.gz
                                    tür başına HCOP dosyaları
    orthologs.tsv                   hazırlanan türlerin ortologları tek tabloda
                                    (varsayılan 10 tür; istenen türler eklenir)
    manifest.json                   hazırlık tamamlandığında yazılır
```

Sorgularda her anahtar kelime bir sütundur. Değer dört biçimde olabilir: metin ya da sayı (hücrenin `|` ile ayrılmış parçalarından biri, büyük/küçük harf duyarsız, eşitse), liste ya da küme (öğelerinden biri tutarsa), `re.compile(...)` (parçalardan birinde desen bulunursa), fonksiyon (hücrenin tam metniyle çağrılır, `True` dönerse). Birden çok filtre VE ile bağlanır. `None` verilen filtre yok sayılır.

<!-- toc -->
- [1. Sabitler](#1-sabitler)
- [2. Sınıf](#2-sınıf): [`Hgnc`](#hgnc) · [`hgnc_dir`](#hgnchgnc_dir) · [`is_ready`](#hgncis_ready) · [`update`](#hgncupdate) · [`prepared_species`](#hgncprepared_species) · [`species_name`](#hgncspecies_name) · [`columns`](#hgnccolumns) · [`query`](#hgncquery) · [`symbols`](#hgncsymbols) · [`values`](#hgncvalues) · [`orthologs`](#hgncorthologs) · [`write_tsv`](#hgncwrite_tsv)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler

`HGNC_URL` dosya sunucusudur. `HCOP_SPECIES` HCOP'un sunduğu 19 türü NCBI Taxonomy ID'leriyle eşler; taxid'ler HCOP'un birleşik dosyasındaki değerlerdir. HCOP iki mayada suş taxid'i kullanır (fisyon mayası 284812, tomurcuklanan maya 4932). `HCOP_TAXID_ALIASES` bu türlerin öteki yaygın taxid'lerini (4896, 559292) aynı HCOP türüne bağlar. `ORTHOLOG_SPECIES` `update`'in varsayılan olarak hazırladığı, insana yakın on memelidir. `ORTHOLOG_COLUMNS` ortolog tablosunun sütunlarıdır: `species`, `taxon_id`, insan sütunları (`_HCOP_HUMAN`), ortolog sütunları (`_HCOP_ORTHOLOG`) ve `support`. `Filter` tip takma adıdır.

```python
import re
from crn.fetch import Hgnc
from crn.fetch.hgnc import HCOP_SPECIES, ORTHOLOG_SPECIES

print(list(ORTHOLOG_SPECIES))  # -> ['chimpanzee', 'macaque', 'mouse', 'rat', 'dog', 'cat', 'cattle', 'pig', 'horse', 'opossum']
print(len(HCOP_SPECIES), [n for n in HCOP_SPECIES if n not in ORTHOLOG_SPECIES])  # -> 19 ['platypus', 'chicken', 'anole_lizard', 'xenopus', 'zebrafish', 'fruitfly', 'c.elegans', 's.cerevisiae', 's.pombe']
```

## 2. Sınıf

### Hgnc

`Hgnc(sys_dir=None, base_url=HGNC_URL, verbose=True)`

`sys_dir` veri klasörü, `base_url` sunucu, `verbose` ilerleme çubuklarıdır.

```python
hgnc = Hgnc(verbose=False)
print(hgnc.SOURCE, hgnc.is_ready())  # -> HGNC True
```

### Hgnc.hgnc_dir

`Hgnc.hgnc_dir  (özellik)`

Veri klasörü: `sys/hgnc`.

### Hgnc.is_ready

`Hgnc.is_ready()`

Gen ve ortolog tabloları hazırsa `True` verir.

### Hgnc.update

`Hgnc.update(force=False, species=None)`

HGNC tablosunu ve HCOP ortologlarını indirip hazırlar. Varsayılan olarak `ORTHOLOG_SPECIES` içindeki on tür hazırlanır, daha önce eklenmiş türler korunur. `species` ayrıca eklenecek HCOP türleridir; ad ya da taxid olabilir (`["chicken", 7955]`). HCOP'ta olmayan tür `ValueError` verir. Her dosya sunucunun MD5'i ile doğrulanır, var olanlar yeniden indirilmez, `orthologs.tsv` yalnızca tür kümesi değiştiğinde yeniden üretilir. `force=True` her şeyi yeniler. Veri klasörünün yolunu döndürür.

<!-- calistirma -->
```python
Hgnc().update()
Hgnc().update(species=["chicken", 7955])      # tavuk ve zebra balığı ortologlarını ekler
```

### Hgnc.prepared_species

`Hgnc.prepared_species()`

Ortolog tablosunda hazır olan türleri (ad → taxid) verir; veri hazır değilse boş sözlük döner.

```python
print(len(hgnc.prepared_species()) >= 10, hgnc.prepared_species()["mouse"])  # -> True 10090
```

### Hgnc.species_name

`statik: Hgnc.species_name(value)`

HCOP tür adını ya da taxid'i (`value`) HCOP tür adına çevirir. HCOP'ta olmayan tür, desteklenen türleri taxid'leriyle sayan `ValueError` verir.

```python
print(Hgnc.species_name(9031), Hgnc.species_name("taxid:7955"), Hgnc.species_name(4896))  # -> chicken zebrafish s.pombe
```

### Hgnc.columns

`Hgnc.columns  (özellik)`

Gen tablosunun sütunlarını verir.

```python
print(hgnc.columns[:5])  # -> ['hgnc_id', 'symbol', 'name', 'locus_group', 'locus_type']
```

### Hgnc.query

`Hgnc.query(fields=None, **filters)`

Gen tablosunda filtrelere uyan satırları verir. `fields` döndürülecek sütunlardır, `None` ise hepsi. `**filters` sütun=değer çiftleridir. Bilinmeyen sütun `ValueError` verir.

```python
print(hgnc.query(gene_group="Aquaporins", fields=["symbol", "entrez_id"])[0])  # -> {'symbol': 'AQP1', 'entrez_id': '358'}
print(len(hgnc.query(symbol=re.compile("^AQP"), status="Approved")) > 10)  # -> True
```

### Hgnc.symbols

`Hgnc.symbols(**filters)`

Filtrelere (`**filters`) uyan genlerin onaylı sembollerini liste olarak verir.

```python
print(hgnc.symbols(symbol=["TP53", "BRCA1", "YOKGEN"]))  # -> ['BRCA1', 'TP53']
print(len(hgnc.symbols(locus_type="gene with protein product", status="Approved")))  # -> 19299
```

### Hgnc.values

`Hgnc.values(column, **filters)`

`column` sütunundaki (`|` ile ayrılmış) değerleri ve gen sayılarını çoktan aza verir. Hangi değerlerle sorgu yapılabileceğini görmek için kullanılır. `**filters` satırları önceden süzer.

```python
print(list(hgnc.values("locus_group").items())[0])  # -> ('protein-coding gene', 19299)
```

### Hgnc.orthologs

`Hgnc.orthologs(fields=None, **filters)`

HCOP tablosunda filtrelere uyan satırları verir. `fields` ve `**filters` `query` ile aynıdır, sütunlar `ORTHOLOG_COLUMNS`'tur. Tür `species` (HCOP adı) ya da `taxid` (NCBI Taxonomy ID; `taxon_id` sütununa bakar) ile seçilir; ikisi de tek değer ya da liste olabilir. İstenen tür henüz hazırlanmamışsa yalnızca o türün HCOP dosyası indirilir ve tabloya eklenir; HCOP'ta olmayan tür `ValueError` verir. Bir insan genine birden çok satır dönebilir: aynı ortolog farklı kaynaklarca ayrı satırlarda desteklenebilir (`support`), bir türde gen kopyası varsa (zebra balığında `mstna`, `mstnb`) iki ortolog çıkar.

```python
satir = hgnc.orthologs(human_symbol="AQP1", species=["mouse"], fields=["species", "ortholog_symbol"])[0]
print(satir)  # -> {'species': 'mouse', 'ortholog_symbol': 'Aqp1'}
```

Hazırlanmamış bir tür ilk istekte indirilir (tavuk dosyası birkaç MB'tır):

<!-- calistirma -->
```python
for r in hgnc.orthologs(human_symbol="MSTN", taxid=[9031, 7955],
                        fields=["species", "ortholog_symbol", "ortholog_entrez_gene"]):
    print(r)
# {'species': 'chicken', 'ortholog_symbol': 'MSTN', 'ortholog_entrez_gene': '373964'}  (birden çok destek satırı)
# {'species': 'zebrafish', 'ortholog_symbol': 'mstnb', 'ortholog_entrez_gene': '798441'}
# {'species': 'zebrafish', 'ortholog_symbol': 'mstna', 'ortholog_entrez_gene': '445816'}
```

### Hgnc.write_tsv

`statik: Hgnc.write_tsv(rows, file_name)`

Sorgu sonucunu (`rows`) `file_name` dosyasına sekmeyle ayrılmış yazar.

```python
Hgnc.write_tsv(hgnc.query(gene_group="Aquaporins", fields=["symbol"]), "aqp.tsv")
print(len(open("aqp.tsv").readlines()) > 10)  # -> True
```

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_GENES` | Gen tablosunun dosya adı |
| `_ORTHOLOGS` | Ortolog tablosunun dosya adı |
| `_OUTPUTS` | `is_ready`'nin baktığı dosyalar |
| `_values` | Tek değeri listeye sarar; liste, küme ve demetleri listeye çevirir |
| `_matches` | Bir hücrenin sorgu değerine uyup uymadığına karar verir (kurallar sayfa başında) |
| `Hgnc._download` | `.part` dosyasına indirir ve `x-goog-hash` başlığındaki MD5 ile doğrular |
| `Hgnc._merge_orthologs` | Tür başına HCOP dosyalarını ortak sütun adlarıyla birleştirir |
| `Hgnc._rows` | Tabloyu satır satır okur; veri yoksa önce `update` çağırır |
| `Hgnc._select` | Sütun denetimi ve filtre döngüsü |
| `Hgnc.__repr__` | `Hgnc('sys/hgnc')` |
