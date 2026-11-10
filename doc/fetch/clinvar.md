# `crn.fetch.clinvar`: klinik varyantlardan gene ve hastalığa

Kaynak: [`crn/fetch/clinvar.py`](../../crn/fetch/clinvar.py)

ClinVar, varyantların klinik yorumlarını laboratuvarlardan ve uzman panellerinden toplar. `Clinvar` bu veriyi varyant düzeyinde sorgular: klinik anlama, inceleme düzeyine, hastalığa, varyant tipine, protein etkisine ve konuma göre süzer. Süzülen varyantlardan gen ve hastalık listeleri çıkarır, genleri `Hpo` ile aynı yoldan `Gene`, `DNA`, `RNA` ve `PRT` nesnelerine çevirir. İki kaynak birbirini tamamlar: HPO bir genin hangi fenotiplerle ilişkili olduğunu, ClinVar o gende hangi varyantların hastalığa yol açtığının bildirildiğini söyler.

Sürüm, NCBI'ın aylık `variant_summary` arşivinin ayıdır (`2026-11`). Aynı sürüm her zaman aynı veriyi verir; haftalık güncellenen dosya kullanılmaz. Arşiv yaklaşık 450 MB'tır. 4,7 milyon varyant içeren veri bir kez SQLite veritabanına hazırlanır (GRCh38 için 1 GB), sorgular bu veritabanından yapılır ve belleğe yalnızca sonuç alınır.

```text
sys/clinvar/<sürüm>/                       ör. sys/clinvar/2026-11/
    raw/variant_summary_<sürüm>.txt.gz     NCBI arşiv dosyası
    clinvar_<derleme>.sqlite               hazırlanmış veritabanı (GRCh38 ya da GRCh37)
    manifest.json                          hazırlık tamamlandığında yazılır
```

<!-- toc -->
- [1. Yorumlarken bilinmesi gerekenler](#1-yorumlarken-bilinmesi-gerekenler)
- [2. Filtreler](#2-filtreler)
- [3. Sabitler ve dönüştürücüler](#3-sabitler-ve-dönüştürücüler): [`significance_classes`](#significance_classes) · [`review_stars`](#review_stars) · [`protein_effect`](#protein_effect)
- [4. Sürüm ve kurulum](#4-sürüm-ve-kurulum): [`Clinvar`](#clinvar) · [`clinvar_dir`](#clinvarclinvar_dir) · [`release`](#clinvarrelease) · [`release_dir`](#clinvarrelease_dir) · [`database`](#clinvardatabase) · [`installed_releases`](#clinvarinstalled_releases) · [`is_ready`](#clinvaris_ready) · [`latest_release`](#clinvarlatest_release) · [`update`](#clinvarupdate) · [`summary`](#clinvarsummary)
- [5. Sorgular](#5-sorgular): [`variants`](#clinvarvariants) · [`genes`](#clinvargenes) · [`symbols`](#clinvarsymbols) · [`conditions`](#clinvarconditions)
- [6. crn nesnelerine dönüşüm](#6-crn-nesnelerine-dönüşüm): [`annotate`](#clinvarannotate) · [`to_genes`](#clinvarto_genes) · [`to_sequences`](#clinvarto_sequences) · [`pack`](#clinvarpack)
- [7. İç yardımcılar](#7-iç-yardımcılar)
<!-- /toc -->

## 1. Yorumlarken bilinmesi gerekenler

Klinik anlam, ClinVar'ın germline toplam sınıflandırmasıdır ve on sınıfa ayrılır (`SIGNIFICANCES`). Birleşik değerler birden çok sınıfa girer: `Pathogenic/Likely pathogenic` hem `pathogenic` hem `likely_pathogenic` sayılır. Gönderenlerin çeliştiği varyantlar (`Conflicting classifications of pathogenicity`) yalnızca `conflicting` sınıfındadır ve patojen sayılmaz. Bu nedenle `significance="pathogenic"` sonucu, laboratuvarların uzlaştığı patojen varyantlarla sınırlıdır. Somatik ve onkojenite sınıflandırmaları ayrı sütunlarda tutulur ve bu süzgece girmez.

İnceleme yıldızı (0-4) sınıflandırmanın ne kadar desteklendiğini gösterir: 0 ölçüt bildirilmemiş, 1 tek gönderen ya da çelişki, 2 çok gönderen ve çelişkisiz, 3 uzman paneli, 4 uygulama kılavuzu. Çalışmaya alınacak varyantlar için en az 1, tercihen 2 yıldız istenmelidir.

Varyantın gene atanması ClinVar'ın gen listesine dayanır. Büyük kopya sayısı değişimleri onlarca, yüzlerce gen kapsar ve gen düzeyindeki sayımları şişirir. ClinVar bu varyantlarda gen alanına serbest metin yazar: `subset of 107 genes: CRKL:SMARCB1:TBX1` yalnızca doz duyarlı genleri listeler, `covers 17 genes, none of which ...` hiçbirini listelemez. Hazırlıkta listelenen semboller varyanta bağlanır, `n_genes` ise metindeki toplamdır (107, 17). Bu yüzden `max_genes` varsayılan olarak 1'dir, yani yalnızca tek gen içeren varyantlar sayılır. `LOC` ile başlayan tanımlanmamış lokuslar sayıma girmez. Bu varsayılan, iki komşu geni birlikte listeleyen (ATXN8 ve ATXN8OS gibi) varyantları da dışarıda bırakır; böyle genlerde `max_genes=None` ya da daha büyük bir değer verilmelidir.

`effect`, varyantın adındaki HGVS protein yazımından (`p.Met694Val`) kurallarla çıkarılır (bkz. [`protein_effect`](#protein_effect)). Protein değişimi yazılmamış varyantlarda (kodlamayan, intronik, büyük yapısal varyantlar) boştur. Bir varyant etki tahmin aracının (VEP gibi) yerine geçmez: ekleme bölgesini bozan bir eşanlamlı varyant yine `synonymous` görünür.

Hastalık adları ClinVar'a gönderildiği biçimdedir ve aynı hastalık farklı adlarla geçebilir. Adın bir parçasıyla arama geniş sonuç verir (2026-11'de `"fever"` 14 hastalığa, `"Marfan"` araması patojen varyantlarıyla 31 gene uyar; `OMIM:154700` 6 gen verir). Kesin seçim için kimlik kullanılmalıdır: `OMIM:249100`, `MONDO:0018088`, `MedGen:C0031069`.

Bir varyant iki derlemede de konumlanmışsa yalnızca seçilen derlemenin konumu tutulur. Yalnızca öteki derlemede konumu olan varyantlar o derlemenin konumuyla kalır (`assembly` sütunu). Bölge filtresi yalnızca seçilen derlemedeki konumlara uygulanır.

## 2. Filtreler

`variants` ve onu kullanan bütün metotlar (`genes`, `symbols`, `conditions`, `annotate`, `to_genes`, `to_sequences`, `pack`) aynı filtreleri alır. Hepsi birlikte sağlanmalıdır (VE); liste verilen filtrelerde öğelerden biri tutması yeter (VEYA).

| Filtre | Anlamı |
|---|---|
| `gene` | Sembol (büyük/küçük harf duyarsız) ya da NCBI Gene ID; bulunamayan gen `KeyError` verir |
| `condition` | Hastalık adının bir parçası, `re.compile` deseni ya da kimlik (`OMIM:249100`) |
| `significance` | `SIGNIFICANCES` sınıfları |
| `min_stars` | En az bu inceleme yıldızı (0-4) |
| `type` | ClinVar varyant tipi, ör. `"single nucleotide variant"`, `"Deletion"`, `"copy number loss"` |
| `origin` | Köken, parça eşleşmesiyle: `germline`, `somatic`, `de novo` ... |
| `effect` | `EFFECTS` protein etkileri |
| `chromosome`, `start`, `end` | Seçili derlemede bölgeyle örtüşen varyantlar |
| `max_genes` | En çok bu kadar gen içeren varyantlar (varsayılan 1); `None` sınırsız |

## 3. Sabitler ve dönüştürücüler

`CLINVAR_URL` NCBI dosya sunucusudur. `ASSEMBLIES` desteklenen derlemeler, `SIGNIFICANCES` klinik anlam sınıfları, `EFFECTS` protein etkileridir. `REVIEW_STARS` inceleme durumlarını yıldıza çevirir. `UNINFORMATIVE_CONDITIONS` hastalık sayılmayan koşul adlarıdır (`not provided`, `not specified`, `see cases`). `SCHEMA_VERSION` veritabanı biçiminin sürümüdür; değiştiğinde eski veritabanı hazır sayılmaz ve yeniden kurulur. `ConditionQuery` ve `Filter` tip takma adlarıdır.

```python
import re
from crn.fetch import Clinvar, Gencode
from crn.fetch.clinvar import SIGNIFICANCES, EFFECTS, significance_classes, review_stars, protein_effect

print(SIGNIFICANCES[:3], EFFECTS[:3])  # -> ('pathogenic', 'likely_pathogenic', 'uncertain') ('missense', 'nonsense', 'synonymous')
```

### significance_classes

`significance_classes(text)`

ClinVar'ın toplam sınıflandırma metnini (`text`) `SIGNIFICANCES` sınıflarına çevirir. Değer `/` ve `;` ile birleşik gelir. Tanınmayan terim `other` olur; hazırlık sırasında tanınmayan terimler sayılır ve manifest'e yazılır.

```python
print(significance_classes("Pathogenic/Likely pathogenic; risk factor"))  # -> ['pathogenic', 'likely_pathogenic', 'risk']
print(significance_classes("Conflicting classifications of pathogenicity"))  # -> ['conflicting']
```

### review_stars

`review_stars(status)`

İnceleme durumunu (`status`) ClinVar yıldızına (0-4) çevirir; tanınmayan durum 0 sayılır.

```python
print(review_stars("criteria provided, multiple submitters, no conflicts"), review_stars("reviewed by expert panel"))  # -> 2 3
```

### protein_effect

`protein_effect(hgvs_p)`

HGVS protein yazımından (`hgvs_p`, `p.` sonrası) etkiyi çıkarır. Kurallar sırayla uygulanır: `=` eşanlamlı, `fs` çerçeve kayması, `ext` ya da `Ter`'in başka amino aside dönmesi stop kaybı, `Ter`'le biten tek değişim ya da `Ter` içeren ekleme anlamsız (nonsense), `Met1` değişimi başlangıç kaybı, tek amino asit değişimi yanlış anlamlı (missense), `del`/`ins`/`dup` çerçeve içi indel, geri kalanı `other`. Yazım yoksa ya da `?` ise boş döner.

```python
print([protein_effect(p) for p in ("Met694Val", "Glu314=", "Gln232Ter", "Leu473fs", "Met1?", "Glu315del")])  # -> ['missense', 'synonymous', 'nonsense', 'frameshift', 'start_lost', 'inframe_indel']
```

## 4. Sürüm ve kurulum

### Clinvar

`Clinvar(release=None, assembly='GRCh38', sys_dir=None, verbose=True)`

`release` `"2026-11"` biçiminde arşiv ayıdır. Verilmezse `sys/` içindeki en yeni hazır sürüm kullanılır, hiç yoksa `update()` en son arşivi indirir. Biçim yanlışsa `ValueError` verir. `assembly` konumların alınacağı derlemedir (`"GRCh38"` ya da `"GRCh37"`). İki derleme aynı arşivden ayrı veritabanlarına hazırlanır. `sys_dir` veri klasörü, `verbose` ilerleme çubuklarıdır.

```python
cv = Clinvar(verbose=False)
print(cv.release, cv.assembly)  # -> 2026-11 GRCh38
```

### Clinvar.clinvar_dir

`Clinvar.clinvar_dir  (özellik)`

ClinVar sürümlerinin tutulduğu klasör (`sys/clinvar`).

### Clinvar.release

`Clinvar.release  (özellik)`

Kullanılan sürüm. Belirtilmediyse `sys/` içindeki, seçili derleme için hazırlanmış en yeni sürümdür; hiç yoksa `None`.

### Clinvar.release_dir

`Clinvar.release_dir  (özellik)`

Sürümün klasörü; sürüm belirli değilse `ValueError` verir.

### Clinvar.database

`Clinvar.database  (özellik)`

Seçili sürüm ve derlemenin SQLite dosyası.

```python
print(cv.database.name)  # -> clinvar_GRCh38.sqlite
```

### Clinvar.installed_releases

`Clinvar.installed_releases()`

`sys/` içinde seçili derleme için hazırlanmış sürümler, eskiden yeniye.

### Clinvar.is_ready

`Clinvar.is_ready()`

Kullanılan sürümün veritabanı seçili derleme ve güncel şema sürümü için hazırsa `True` döner.

```python
print(cv.installed_releases(), cv.is_ready())  # -> ['2026-11'] True
```

### Clinvar.latest_release

`Clinvar.latest_release()`

Sunucudaki en son aylık arşivi döndürür; ağ bağlantısı gerekir. Arşiv listesi okunamazsa `RuntimeError` verir.

### Clinvar.update

`Clinvar.update(force=False)`

Sürümün arşivini indirir ve seçili derleme için veritabanını hazırlar. Nesne bir sürümle kurulduysa o sürüm alınır ve hazırsa ağa çıkılmaz; kurulmadıysa sunucudaki en son arşiv alınır. Eski yılların arşivleri alt klasörlerde aranır. NCBI arşiv dosyaları için MD5 özeti yayımlamaz. Bu yüzden indirme dosya boyutuyla, içerik ise açılırken gzip CRC32 denetimiyle doğrulanır; dosyanın SHA-256 özeti manifest'e yazılır. Bozuk dosya silinir ve `RuntimeError` verilir. İndirilmiş dosya yeniden indirilmez, öteki derleme aynı dosyadan hazırlanır. `force=True` dosyayı yeniden indirir ve veritabanını yeniden kurar. Hazırlanan sürümü döndürür.

Hazırlıkta her varyant (AlleleID) bir kez yazılır, konum önce seçili derlemeden alınır. Metin sütunları (tip, klinik anlam, inceleme durumu) ve hastalıklar ayrı tablolarda bir kez tutulur, varyant satırı yalnızca kodlarını taşır.

<!-- calistirma -->
```python
Clinvar(release="2026-11", assembly="GRCh37").update()
```

### Clinvar.summary

`Clinvar.summary()`

Sürüm, derleme, varyant, gen ve hastalık sayıları, her klinik anlam sınıfındaki varyant sayısı ve hazırlıkta tanınmayan terim sayıları. Değerler hazırlık sırasında manifest'e yazılır, sorgu yapılmaz.

```python
ozet = cv.summary()
print(ozet["variants"], ozet["genes"], ozet["pathogenic"], ozet["likely_pathogenic"], ozet["unknown_significance"])  # -> 4707240 29360 265706 183241 0
```

## 5. Sorgular

### Clinvar.variants

`Clinvar.variants(gene=None, *, limit=None, **filters)`

Filtrelere uyan varyantlar, AlleleID sırasıyla. `gene` sembol ya da NCBI Gene ID'dir (liste olabilir), `limit` en çok satır sayısı, `**filters` [filtrelerdir](#2-filtreler). Satır sütunları: `allele_id`, `variation_id`, `gene_symbol` (çok genli varyantta `;` ile ayrılmış), `gene_id`, `n_genes`, `name`, `transcript`, `hgvs_c`, `hgvs_p`, `effect`, `type`, `significance`, `classes` (`|` ile ayrılmış), `stars`, `review_status`, `last_evaluated` (YYYY-AA-GG), `origin`, `rs`, `submitters`, `assembly`, `chromosome`, `start`, `stop`, `position_vcf`, `ref`, `alt` (VCF yazımı), `conditions` ve `condition_ids` (`|` ile ayrılmış). Gen bulunamazsa `KeyError`, filtre değeri geçersizse `ValueError` verir.

```python
m694v = [r for r in cv.variants("MEFV", significance="pathogenic", min_stars=1) if r["hgvs_p"] == "Met694Val"][0]
print(m694v["transcript"], m694v["hgvs_c"], m694v["effect"], m694v["classes"], m694v["stars"])  # -> NM_000243.3 c.2080A>G missense pathogenic|likely_pathogenic 2
print(len(cv.variants("MEFV")), len(cv.variants("MEFV", significance=["pathogenic", "likely_pathogenic"], min_stars=1)))  # -> 1425 23
print(len(cv.variants(effect="synonymous", significance="pathogenic")))  # -> 474
```

### Clinvar.genes

`Clinvar.genes(condition=None, **filters)`

Filtrelere uyan varyantların genleri, en çok varyantı olan önce. `condition` koşul filtresidir. Satır sütunları: `ncbi_gene_id`, `gene_symbol`, `variants`, `pathogenic` (patojen ya da olası patojen), `uncertain`, `benign` (benign ya da olası benign), `conflicting`, `max_stars`. Sayılar filtreden geçen varyantlar üzerindendir: `significance="pathogenic"` verilirse `uncertain` sütunu 0 olur. `gene` filtresi verilirse çok genli varyantların öteki genleri listelenmez.

```python
print(cv.genes("OMIM:249100", significance="pathogenic"))  # -> [{'ncbi_gene_id': '4210', 'gene_symbol': 'MEFV', 'variants': 12, 'pathogenic': 12, 'uncertain': 0, 'benign': 0, 'conflicting': 0, 'max_stars': 2}]
print([r["gene_symbol"] for r in cv.genes(significance=["pathogenic", "likely_pathogenic"], min_stars=2)[:5]])  # -> ['BRCA2', 'BRCA1', 'NF1', 'ATM', 'TTN']
```

### Clinvar.symbols

`Clinvar.symbols(condition=None, **filters)`

Filtrelere uyan genlerin sembolleri, alfabetik. `condition` koşul filtresidir, öteki filtreler `genes` ile aynıdır. Sonuç doğrudan `CUB` ya da `crn cub table -` girdisi olarak kullanılabilir.

```python
kardiyo = cv.symbols("cardiomyopathy", significance=["pathogenic", "likely_pathogenic"], min_stars=2, effect="nonsense")
print(len(kardiyo), kardiyo[:5])  # -> 36 ['ABCC9', 'ALMS1', 'ALPK3', 'BAG3', 'COQ4']
```

### Clinvar.conditions

`Clinvar.conditions(name=None, *, informative=True, limit=None, **filters)`

Hastalıklar: adında `name` geçenler ve/veya filtrelere uyanlar, en çok varyantı olan önce. `name` adın bir parçası, `re.compile` deseni ya da kimliktir. `informative` `True` ise `not provided`, `not specified` ve `see cases` çıkarılır. `limit` en çok satır sayısı, `**filters` varyant filtreleridir (ör. `gene="MEFV"`). Satır sütunları: `condition_name`, `condition_ids`, `variants`, `genes`; sayılar filtreden geçen varyantlar üzerindendir.

```python
print([(r["condition_name"], r["variants"], r["genes"]) for r in cv.conditions("mediterranean fever")])  # -> [('Familial Mediterranean fever', 1187, 1), ('Familial Mediterranean fever, autosomal dominant', 344, 3)]
print([r["condition_name"] for r in cv.conditions(gene="MEFV", limit=3)])  # -> ['Familial Mediterranean fever', 'Familial Mediterranean fever, autosomal dominant', 'Acute febrile neutrophilic dermatosis']
```

## 6. crn nesnelerine dönüşüm

### Clinvar.annotate

`Clinvar.annotate(obj, gene=None, **filters)`

`obj` nesnesine (`Gene`, `DNA`, `RNA` ya da `PRT`) genin ClinVar özetini yazar: `clinvar_condition` (bilgilendirici hastalık adları, tekrarsız), `clinvar_variants` (filtreden geçen varyant sayısı), `clinvar_pathogenic` (bunlardan patojen ya da olası patojen olanlar) ve `clinvar_release`. Aynı anahtarların eski değerleri silinir, değerler CRN dosyasında korunur. `gene` sembolü verir; verilmezse features'taki `gene_name`, o da yoksa `Gene` adı kullanılır. Sembol belirlenemezse `ValueError`, gen ClinVar'da yoksa `KeyError` verir. `**filters` varyant filtreleridir. Aynı nesneyi döndürür.

```python
gc = Gencode(verbose=False)
gen = cv.annotate(gc.canonical("MEFV"), min_stars=1)
print(gen.features.first("clinvar_release"), gen.features.first("clinvar_pathogenic") > 0, gen.features.get("clinvar_condition")[:2])  # -> 2026-11 True ['Inborn genetic diseases', 'Autoinflammatory syndrome']
```

### Clinvar.to_genes

`Clinvar.to_genes(condition=None, source=None, group=None, key='group', **filters)`

Filtrelere uyan genleri `Gene` listesi olarak çıkarır. Her gen bir DNA, bir transkript ve bir protein içerir, yani doğrudan `CUB`'a verilebilir. Genin ClinVar özeti features'a yazılır (bkz. `annotate`). Protein kodlamayan ya da proteini transkriptle eşleşmeyen genler atlanır ve `UserWarning` ile bildirilir. `condition` koşul filtresidir. `source` `Gencode` (ana transkript) ya da `Ncbi` (NM_ öncelikli ilk izoform) olabilir, `None` ise `sys/` içindeki GENCODE kullanılır. `group` verilirse her genin features'ına `key` anahtarıyla yazılır. `**filters` varyant filtreleridir.

```python
genler = cv.to_genes("OMIM:249100", source=gc, significance="pathogenic", group="fmf", key="grup")
print([g.name for g in genler], genler[0].features.first("grup"))  # -> ['MEFV_ENST00000219596.6'] fmf
```

### Clinvar.to_sequences

`Clinvar.to_sequences(kind, condition=None, source=None, **filters)`

Filtrelere uyan genlerin dizilerini verir, sonuç `write_fasta` ile yazılabilir. `kind` `"dna"` (genomik), `"rna"` (ana transkript) ya da `"protein"`dir, başka değer `ValueError` verir. `condition`, `source` ve `**filters` `to_genes` ile aynıdır.

```python
dizi = cv.to_sequences("protein", "OMIM:249100", source=gc, significance="pathogenic")
print(type(dizi[0]).__name__, len(dizi[0].seq), dizi[0].features.first("clinvar_release"))  # -> PRT 781 2026-11
```

### Clinvar.pack

Taban sınıftan ([`GeneSetSource.pack`](_geneset.md#genesetsourcepack)) gelir; `Clinvar`'da grup değeri sözlük değilse `condition` filtresine aktarılır. Bütün gruplara uygulanacak filtreler (`significance`, `min_stars` ...) anahtar kelime olarak verilir. `write_tsv` de aynı tabandan gelir.

```python
rapor = cv.pack({"fmf": "OMIM:249100", "gpx": {"gene": "GPX4"}}, "clinvar_gruplar.crn", source=gc, significance="pathogenic")
print(rapor["release"], rapor["groups"])  # -> 2026-11 {'fmf': 1, 'gpx': 1}
```

## 7. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_SIGNIFICANCE_TERMS` | ClinVar terimlerinden sınıflara eşleme |
| `_PATHOGENIC`, `_BENIGN` | `genes` sayımlarında birleştirilen sınıflar |
| `_FEATURE_KEYS` | `annotate`'in yazdığı anahtarlar |
| `Clinvar.SOURCE` | features'taki `source` değeri (`ClinVar`) |
| `Clinvar.FEATURE_KEYS`, `Clinvar.RELEASE_KEY`, `Clinvar.GROUP_FILTER` | Taban sınıfın kancaları (bkz. [_geneset.md](_geneset.md#1-alt-sınıfın-sağladıkları)) |
| `_ASSEMBLY_RANK` | Konumun hangi derleme satırından alınacağını belirleyen sıra |
| `_MONTHS` | Tarih dönüştürmede ay kısaltmaları; yerel ayardan bağımsızdır |
| `_NAME` | Varyant adından transkript, `c.` ve `p.` yazımını ayıran düzenli ifade |
| `_AA` | Üç harfli amino asit deseni |
| `_CONDITION_ID` | Koşul sorgusunun kimlik olup olmadığını ayıran desen |
| `_REQUIRED` | Arşivde bulunması gereken sütunlar; eksikse hazırlık durur |
| `_SCHEMA`, `_INDEXES` | Veritabanı tabloları ve dizinleri |
| `_UPSERT` | Varyant ekleme; aynı AlleleID'de daha uygun derlemenin konumunu yazar |
| `_SELECT` | `variants` sorgusunun kod tablolarını ve hastalıkları birleştiren gövdesi |
| `_ROW_ORDER` | `variants` satırlarının sütun sırası |
| `_as_list` | Tek değeri listeye sarar, `None`'ı boş liste yapar |
| `_mask`, `_unmask` | Sınıf listesi ile bit maskesi arasında dönüşüm |
| `_date` | `"Dec 17, 2024"` → `"2024-12-17"` |
| `_condition_ids` | Kimlikleri birleştirir; `MONDO:MONDO:x` yazımını `MONDO:x` yapar |
| `_conditions` | `PhenotypeIDS` ve `PhenotypeList`'ten (ad, kimlikler) çiftleri; koşul sayıları tutmazsa adlar kimliksiz döner |
| `_gene_symbols` | `;` ile ayrılmış sembolleri (CNV metninde `:` ile ayrılmış doz duyarlı genleri) tekrarsız listeye çevirir |
| `_CNV_GENES` | CNV'lerin `subset of N genes: ...` ve `covers N genes, ...` gen metni |
| `_gene_count` | Gen sayısı; CNV metninde ClinVar'ın toplamı, öbür durumlarda `LOC` lokusları sayılmaz |
| `_int`, `_allele`, `_dash` | Sayı, alel ve boş değer dönüştürücüleri |
| `Clinvar.__repr__` | `Clinvar('2026-11', 'GRCh38', '<sys>/clinvar')` |
| `Clinvar._normalize_release` | Sürüm biçimini denetler |
| `Clinvar._manifest` | Sürüm klasörünün manifest'ini okur |
| `Clinvar._prepared` | Klasörde seçili derlemenin güncel şemayla hazırlanıp hazırlanmadığı |
| `Clinvar._open` | `User-Agent` ile HTTP isteği açar |
| `Clinvar._locate` | Arşiv dosyasının adresini ve boyutunu bulur; 404'te yıl alt klasörünü dener |
| `Clinvar._download` | `.part` dosyasına indirir, boyutu denetler, SHA-256 hesaplar |
| `Clinvar._records` | Arşiv satırlarını okur; ilerlemeyi sıkıştırılmış baytlarla izler |
| `Clinvar._prepare` | Arşivi SQLite veritabanına çevirir, sınıf sayılarını ve tanınmayan terimleri döndürür |
| `Clinvar._ensure` | Veri hazır değilse indirir |
| `Clinvar._connection` | Veritabanını salt okunur açar, desen işlevini kaydeder |
| `Clinvar._close` | Bağlantıyı kapatır |
| `Clinvar._regex_match` | `re.compile` koşullarını SQL içinde uygulayan işlev |
| `Clinvar._symbols_for` | Sembol ya da Gene ID'leri ClinVar sembollerine çevirir |
| `Clinvar._condition_clause` | Koşul sorgusunun SQL parçası |
| `Clinvar._where` | Bütün filtrelerin SQL koşulu, parametreleri ve çözülen semboller |
| `Clinvar._describe` | `_build` kancası: kurulan gene ClinVar özetini yazar |

`Clinvar`, `GeneSetSource` tabanından türer; ortak yardımcılar [_geneset.md](_geneset.md) belgesindedir.
