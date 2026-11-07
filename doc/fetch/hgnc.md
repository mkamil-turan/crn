# `crn.fetch.hgnc`: HGNC gen tablosu ve ortologlar

Kaynak: [`crn/fetch/hgnc.py`](../../crn/fetch/hgnc.py)

HGNC insan genlerinin onaylı sembollerini ve gen ailelerini, HCOP ise on türdeki ortologlarını verir. Bir aileyi ya da yolağı gen listesine çevirmek ve hasta genlerinin fare karşılıklarını bulmak için kullanılır.

```text
sys/hgnc/
    hgnc_complete_set.txt           HGNC'nin tüm genleri (sekmeyle ayrılmış)
    raw/human_<tür>_hcop_fifteen_column.txt.gz
                                    tür başına HCOP dosyaları
    orthologs.tsv                   10 türün ortologları tek tabloda
    manifest.json                   hazırlık tamamlandığında yazılır
```

Sorgularda her anahtar kelime bir sütundur. Değer dört biçimde olabilir: metin ya da sayı (hücrenin `|` ile ayrılmış parçalarından biri, büyük/küçük harf duyarsız, eşitse), liste ya da küme (öğelerinden biri tutarsa), `re.compile(...)` (parçalardan birinde desen bulunursa), fonksiyon (hücrenin tam metniyle çağrılır, `True` dönerse). Birden çok filtre VE ile bağlanır. `None` verilen filtre yok sayılır.

<!-- toc -->
- [1. Sabitler](#1-sabitler)
- [2. Sınıf](#2-sınıf): [`Hgnc`](#hgnc) · [`hgnc_dir`](#hgnchgnc_dir) · [`is_ready`](#hgncis_ready) · [`update`](#hgncupdate) · [`columns`](#hgnccolumns) · [`query`](#hgncquery) · [`symbols`](#hgncsymbols) · [`values`](#hgncvalues) · [`orthologs`](#hgncorthologs) · [`write_tsv`](#hgncwrite_tsv)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Sabitler

`HGNC_URL` dosya sunucusudur. `ORTHOLOG_SPECIES` on türü taxon kimlikleriyle eşler. `ORTHOLOG_COLUMNS` ortolog tablosunun sütunlarıdır: `species`, `taxon_id`, insan sütunları (`_HCOP_HUMAN`), ortolog sütunları (`_HCOP_ORTHOLOG`) ve `support`. `Filter` tip takma adıdır.

```python
import re
from crn.fetch import Hgnc
from crn.fetch.hgnc import ORTHOLOG_SPECIES

print(list(ORTHOLOG_SPECIES))  # -> ['chimpanzee', 'macaque', 'mouse', 'rat', 'dog', 'cat', 'cattle', 'pig', 'horse', 'opossum']
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

`Hgnc.update(force=False)`

HGNC tablosunu ve on türün HCOP ortologlarını indirip hazırlar. Her dosya sunucunun MD5'i ile doğrulanır, var olanlar yeniden indirilmez, `orthologs.tsv` ham dosyalardan bir kez üretilir. `force=True` her şeyi yeniler. Veri klasörünün yolunu döndürür.

<!-- calistirma -->
```python
Hgnc().update()
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

On türün HCOP tablosunda filtrelere uyan satırları verir. `fields` ve `**filters` `query` ile aynıdır, sütunlar `ORTHOLOG_COLUMNS`'tur, tür `species` sütunuyla seçilir.

```python
satir = hgnc.orthologs(human_symbol="AQP1", species=["mouse"], fields=["species", "ortholog_symbol"])[0]
print(satir)  # -> {'species': 'mouse', 'ortholog_symbol': 'Aqp1'}
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
| `_matches` | Bir hücrenin sorgu değerine uyup uymadığına karar verir (kurallar sayfa başında) |
| `Hgnc._download` | `.part` dosyasına indirir ve `x-goog-hash` başlığındaki MD5 ile doğrular |
| `Hgnc._merge_orthologs` | Tür başına HCOP dosyalarını ortak sütun adlarıyla birleştirir |
| `Hgnc._rows` | Tabloyu satır satır okur; veri yoksa önce `update` çağırır |
| `Hgnc._select` | Sütun denetimi ve filtre döngüsü |
| `Hgnc.__repr__` | `Hgnc('sys/hgnc')` |
