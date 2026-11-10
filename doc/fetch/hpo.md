# `crn.fetch.hpo`: semptomdan gen grubuna

Kaynak: [`crn/fetch/hpo.py`](../../crn/fetch/hpo.py)

Human Phenotype Ontology, klinik fenotipleri bir ontoloji olarak düzenler ve onları hastalıklar üzerinden genlere bağlar. `Hpo` bir semptomdan gen listesi çıkarır, listeyi filtreler ve `Gene`, `DNA`, `RNA`, `PRT` nesnelerine dönüştürür. Ontoloji sayesinde bir terim alt terimleriyle birlikte sorgulanır: `Seizure` istendiğinde `Focal seizure` ile işaretlenmiş genler de gelir.

```text
sys/hpo/<sürüm>/                    ör. sys/hpo/v2026-09-01/
    raw/hp.obo                      ontoloji: terimler, eş anlamlılar, hiyerarşi
    raw/genes_to_phenotype.txt      gen, terim, sıklık, hastalık
    raw/genes_to_disease.txt        gen, hastalık, ilişki tipi
    raw/phenotype.hpoa              hastalık notları: kanıt, başlangıç, kalıtım
    manifest.json                   hazırlık tamamlandığında yazılır
```

<!-- toc -->
- [1. Yorumlarken bilinmesi gerekenler](#1-yorumlarken-bilinmesi-gerekenler)
- [2. Filtreler](#2-filtreler)
- [3. Sabitler](#3-sabitler): [`frequency_value`](#frequency_value)
- [4. Sürüm ve kurulum](#4-sürüm-ve-kurulum): [`Hpo`](#hpo) · [`hpo_dir`](#hpohpo_dir) · [`release`](#hporelease) · [`release_dir`](#hporelease_dir) · [`installed_releases`](#hpoinstalled_releases) · [`is_ready`](#hpois_ready) · [`latest_release`](#hpolatest_release) · [`update`](#hpoupdate) · [`summary`](#hposummary)
- [5. Ontoloji](#5-ontoloji): [`term`](#hpoterm) · [`ancestors`](#hpoancestors) · [`descendants`](#hpodescendants) · [`search`](#hposearch)
- [6. Sorgular](#6-sorgular): [`annotations`](#hpoannotations) · [`genes`](#hpogenes) · [`symbols`](#hposymbols) · [`phenotypes`](#hpophenotypes) · [`diseases`](#hpodiseases)
- [7. crn nesnelerine dönüşüm](#7-crn-nesnelerine-dönüşüm): [`annotate`](#hpoannotate) · [`to_genes`](#hpoto_genes) · [`to_sequences`](#hpoto_sequences) · [`pack`](#hpopack) · [`write_tsv`](#hpowrite_tsv)
- [8. İç yardımcılar](#8-iç-yardımcılar)
<!-- /toc -->

## 1. Yorumlarken bilinmesi gerekenler

HPO'da gen-fenotip ilişkisi hastalık üzerindendir: gen bir hastalığa, hastalık fenotipe bağlıdır. Bu, genin o fenotipi tek başına yaptığı anlamına gelmez; fenotip hastalığın başka bir bileşeninden kaynaklanıyor olabilir.

`genes_to_phenotype.txt` hastalıkta bulunmadığı açıkça yazılmış (`NOT`) notları ve hiç gözlenmemiş (`0/N`) satırları da içerir. Bunlar yüklenirken çıkarılır, sayıları `summary()` içindedir. Sıklık kategorilerinde alt sınır kullanılır: "Frequent" (%30-79) 0,30 sayılır, yani `min_frequency` geçilen satırın o sıklığa ulaştığını garanti eder.

İlişki etiketlerine (`association_type`) tek başına güvenilmemelidir. SCA2, LRRK2 ve GBA1 gibi bilinen monogenik ilişkiler veride `POLYGENIC` olarak geçer; `association="MENDELIAN"` filtresi bu genleri siler. PolyQ gibi bir gen kümesi kuruyorsanız sonucu bildiğiniz genlerle karşılaştırın.

## 2. Filtreler

`annotations` ve onu kullanan bütün metotlar (`genes`, `symbols`, `diseases`, `to_genes`, `to_sequences`, `pack`, `annotate`, `phenotypes`) aynı filtreleri alır. Hepsi birlikte sağlanmalıdır (VE).

| Filtre | Anlamı |
|---|---|
| `term` | Terim ya da terimler (birleşim); alt terimler dahildir |
| `exclude` | Bu terimlerle ilişkili genleri sonuçtan çıkarır |
| `gene` | Sembol ya da NCBI Gene ID |
| `disease` | OMIM/ORPHA kimliği, hastalığın tam adı ya da `re.compile` |
| `min_frequency` | En az bu sıklık (0-1); sıklığı bilinmeyen satırlar elenir |
| `association` | `MENDELIAN`, `POLYGENIC`, `UNKNOWN` |
| `evidence` | `PCS` (yayın), `TAS` (uzman), `IEA` (otomatik) |
| `onset` | Hastalığın başlangıç terimi, alt terimleriyle |
| `inheritance` | Hastalığın kalıtım terimi, alt terimleriyle |
| `descendants` | `False` ise yalnızca terimin kendisiyle işaretli satırlar |

## 3. Sabitler

`HPO_RELEASES` GitHub sürüm API adresidir. `PHENOTYPE_ROOT` (`HP:0000118`) fenotip dalının kökü, `ONSET_ROOT` (`HP:0003674`) başlangıç dalının köküdür. `ASSOCIATIONS`, `EVIDENCES` ve `SEQUENCE_KINDS` kabul edilen değerlerdir. `FREQUENCY_TERMS` sıklık terimlerini alt sınırlarına çevirir: Obligate 1,0; Very frequent 0,8; Frequent 0,3; Occasional 0,05; Very rare 0,01; Excluded 0. `TermQuery` ve `Filter` tip takma adlarıdır.

### frequency_value

`frequency_value(frequency)`

HPO sıklık yazımını (`frequency`) 0-1 arasına çevirir, bilinmiyorsa `None` verir. `"7/12"`, `"45%"` ve `"HP:0040282"` yazımlarını tanır.

```python
import re
from crn.fetch import Hpo, Gencode
from crn.fetch.hpo import frequency_value

print(round(frequency_value("7/12"), 3), frequency_value("45%"), frequency_value("HP:0040282"), frequency_value(""))  # -> 0.583 0.45 0.3 None
```

## 4. Sürüm ve kurulum

### Hpo

`Hpo(release=None, sys_dir=None, verbose=True)`

`release` `"v2026-09-01"` ya da `"2026-09-01"` biçimindedir; `None` ise `sys/` içindeki en yeni hazır sürüm kullanılır, hiç yoksa `update()` en son sürümü indirir. `sys_dir` veri klasörü, `verbose` ilerleme çubuklarıdır.

```python
hpo = Hpo(verbose=False)
print(hpo.SOURCE, hpo.release)  # -> HPO v2026-09-01
```

### Hpo.hpo_dir

`Hpo.hpo_dir  (özellik)`

Veri klasörü: `sys/hpo`.

### Hpo.release

`Hpo.release  (özellik)`

Kullanılan sürüm.

```python
print(hpo.release)  # -> v2026-09-01
```

### Hpo.release_dir

`Hpo.release_dir  (özellik)`

Sürümün klasörü; sürüm belirli değilse `ValueError` verir.

### Hpo.installed_releases

`Hpo.installed_releases()`

Hazırlığı tamamlanmış sürümleri eskiden yeniye verir.

```python
print(hpo.installed_releases())  # -> ['v2026-09-01']
```

### Hpo.is_ready

`Hpo.is_ready()`

Sürümün bütün dosyaları (`_FILES`) yerinde mi söyler.

```python
print(hpo.is_ready())  # -> True
```

### Hpo.latest_release

`Hpo.latest_release()`

Sunucudaki en son sürümü verir, ağ gerekir.

<!-- calistirma -->
```python
print(hpo.latest_release())
```

### Hpo.update

`Hpo.update(force=False)`

Sürümün dosyalarını indirir. Nesne bir sürümle oluşturulduysa o sürüm (hazırsa ağa çıkılmaz), oluşturulmadıysa en son sürüm alınır. Her dosya GitHub'ın SHA-256 özetiyle doğrulanır. `force=True` yeniden indirir. Hazırlanan sürümü döndürür.

<!-- calistirma -->
```python
Hpo(release="2026-09-01").update(force=True)
```

### Hpo.summary

`Hpo.summary()`

Sürüm, kayıt sayıları ve yükleme sırasında çıkarılan satırları verir.

```python
ozet = hpo.summary()
print(sorted(ozet))  # -> ['annotations', 'diseases', 'genes', 'release', 'removed_not', 'removed_zero_frequency', 'terms']
```

## 5. Ontoloji

### Hpo.term

`Hpo.term(term)`

Terimin kaydını verir: `id, name, definition, synonyms, layperson, parents, children`. `term` HP kimliği, eski kimlik, terim adı ya da eş anlamlıdır. Bulunamazsa, belirsizse ya da kullanımdan kalkmış ve yerine geçen yoksa `ValueError` verir; hata mesajı yakın terimleri önerir.

```python
t = hpo.term("Seizure")
print(t["id"], t["name"], len(t["children"]) > 3)  # -> HP:0001250 Seizure True
```

### Hpo.ancestors

`Hpo.ancestors(term)`

`term` teriminin bütün üst terimlerini HP kimliği olarak verir, kendisi hariç.

```python
print(hpo.ancestors("Seizure")[:3], "HP:0000118" in hpo.ancestors("Seizure"))  # -> ['HP:0000001', 'HP:0000118', 'HP:0000707'] True
```

### Hpo.descendants

`Hpo.descendants(term)`

`term` teriminin bütün alt terimlerini HP kimliği olarak verir, kendisi hariç.

```python
print(len(hpo.descendants("Seizure")) > 50)  # -> True
```

### Hpo.search

`Hpo.search(text, limit=None)`

Adında ya da eş anlamlılarında `text` geçen terimleri verir. Satır `id, name, matched, genes` içerir; `genes` alt terimler dahil gen sayısıdır. Tam eşleşmeler önce, sonra gen sayısına göre sıralanır. `text` büyük/küçük harf duyarsız metin ya da `re.compile` desenidir, `limit` satır sınırıdır.

```python
print(hpo.search("seizure", limit=1)[0]["name"])  # -> Seizure
print(len(hpo.search(re.compile(r"^Ataxia"), limit=3)))  # -> 3
```

## 6. Sorgular

### Hpo.annotations

`Hpo.annotations(term=None, *, exclude=None, gene=None, disease=None, min_frequency=None, association=None, evidence=None, onset=None, inheritance=None, descendants=True)`

Filtrelere uyan gen-terim-hastalık satırlarını verir. Satır `ncbi_gene_id, gene_symbol, hpo_id, hpo_name, frequency, frequency_value, disease_id, disease_name, association_type, evidence` içerir. Parametreler [Filtreler](#2-filtreler) bölümündedir: `term`, `exclude`, `gene`, `disease`, `min_frequency`, `association`, `evidence`, `onset`, `inheritance`, `descendants`. Terim bulunamazsa ya da değer geçersizse `ValueError` verir.

```python
satir = hpo.annotations("Seizure", gene="SCN1A", min_frequency=0.3, evidence=["PCS", "TAS"], association="MENDELIAN")
print(len(satir) > 0, sorted(satir[0])[:4])  # -> True ['association_type', 'disease_id', 'disease_name', 'evidence']
print(len(hpo.annotations("Seizure", exclude="Autism", onset="Infantile onset", inheritance="Autosomal dominant inheritance", descendants=True, disease=re.compile("epilep", re.I))) > 0)  # -> True
```

### Hpo.genes

`Hpo.genes(term=None, **filters)`

Filtrelere uyan genlerin özetini verir, en çok hastalıkla ilişkili olan önce. Satır `ncbi_gene_id, gene_symbol, diseases, terms, max_frequency` içerir; sayılar filtreden geçen satırlardandır. `term` ve `**filters` `annotations` ile aynıdır.

```python
print(sorted(hpo.genes("Ataxia", min_frequency=0.8)[0]))  # -> ['diseases', 'gene_symbol', 'max_frequency', 'ncbi_gene_id', 'terms']
```

### Hpo.symbols

`Hpo.symbols(term=None, **filters)`

Filtrelere uyan genlerin sembollerini alfabetik verir. `term` ve `**filters` `annotations` ile aynıdır.

```python
print(len(hpo.symbols("HP:0001250", min_frequency=0.3, association="MENDELIAN")) > 100)  # -> True
```

### Hpo.phenotypes

`Hpo.phenotypes(gene, **filters)`

`gene` geninin fenotip terimlerini verir; kalıtım ve başlangıç terimleri dahil değildir. Satır `hpo_id, hpo_name, diseases, max_frequency` içerir, en çok hastalıkta geçen önce gelir. `**filters` `annotations` filtreleridir. Gen HPO'da yoksa `KeyError`.

```python
print(sorted(hpo.phenotypes("SCN1A")[0]))  # -> ['diseases', 'hpo_id', 'hpo_name', 'max_frequency']
```

### Hpo.diseases

`Hpo.diseases(name=None, **filters)`

Hastalıkları verir: adında `name` geçenler ve/veya filtrelere (`**filters`) uyanlar. Satır `disease_id, disease_name, genes, terms, inheritance, onset` içerir; son ikisi `|` ile ayrılmış terim adlarıdır. En çok gen önce gelir. `name` metin ya da `re.compile` desenidir.

```python
h = hpo.diseases("parkinson", term="Tremor")
print(len(h) > 0, sorted(h[0])[:3])  # -> True ['disease_id', 'disease_name', 'genes']
```

## 7. crn nesnelerine dönüşüm

### Hpo.annotate

`Hpo.annotate(obj, gene=None, term=PHENOTYPE_ROOT, **filters)`

`obj` nesnesine (`Gene`, `DNA`, `RNA` ya da `PRT`) HPO terimlerini ve hastalıklarını yazar. features'a `hpo_id`, `hpo_name`, `disease_id`, `disease_name` listeleri ve `hpo_release` eklenir, aynı anahtarların eski değerleri silinir. Gen HPO'da yoksa listeler boş kalır. Bu değerler CRN dosyasında korunur. `gene` sembolü verir; verilmezse features'taki `gene_name`, o da yoksa `Gene` adı kullanılır, belirlenemezse `ValueError` verilir. `term` yazılacak terimleri sınırlar (varsayılan bütün fenotipler), `**filters` `annotations` filtreleridir. Aynı nesneyi döndürür.

```python
gc = Gencode(verbose=False)
gen = gc.canonical("SCN1A")
hpo.annotate(gen, min_frequency=0.8)
print(gen.features.first("hpo_release"), len(gen.features.get("hpo_id")) > 0)  # -> v2026-09-01 True
prt = gc.protein("ENST00000269305")
hpo.annotate(prt, gene="TP53", term="Neoplasm")
print(len(prt.features.get("disease_id")) > 0)  # -> True
```

### Hpo.to_genes

`Hpo.to_genes(term=None, source=None, group=None, key='group', **filters)`

Filtrelere uyan genleri `Gene` listesi olarak çıkarır. Her gen bir DNA, bir transkript ve bir protein içerir, yani doğrudan `CUB`'a verilebilir. Eşleşen HPO terimleri ve hastalıklar features'a yazılır. Protein kodlamayan ya da standart kodla çevrilemeyen genler (mitokondri gibi) atlanır ve `UserWarning` ile bildirilir.

`term` terim ya da terimlerdir. `source` `Gencode` ise MANE_Select, yoksa Ensembl_canonical transkript alınır; ORF'u proteine eşleşmezse ilk uyan transkript seçilir ve `canonical_tag` `other_transcript` olur. `source` `Ncbi` ise NM_ öncelikli ilk izoform alınır. `None` ise `sys/` içindeki GENCODE kullanılır. `group` verilirse her genin features'ına `key` anahtarıyla yazılır. `**filters` `annotations` filtreleridir.

```python
genler = hpo.to_genes(gene=["ATXN1", "ATXN3"], source=gc, group="sca", key="grup")
print([g.name.split("_")[0] for g in genler], genler[0].features.first("grup"))  # -> ['ATXN3', 'ATXN1'] sca
```

### Hpo.to_sequences

`Hpo.to_sequences(kind, term=None, source=None, **filters)`

Filtrelere uyan genlerin dizilerini verir, sonuç `write_fasta` ile yazılabilir. `kind` `"dna"` (genomik), `"rna"` (ana transkript) ya da `"protein"`dir, başka değer `ValueError` verir. `term`, `source` ve `**filters` `to_genes` ile aynıdır. Her dizinin features'ında gen adı, HPO terimleri ve hastalıklar vardır.

```python
dizi = hpo.to_sequences("protein", gene="ATXN1", source=gc)
print(type(dizi[0]).__name__, dizi[0].seq[:8])  # -> PRT MKSNQERS
```

### Hpo.pack

Taban sınıftan ([`GeneSetSource.pack`](_geneset.md#genesetsourcepack)) gelir; `Hpo`'da grup değeri sözlük değilse `term` filtresine aktarılır.

`Hpo.pack(groups, file_name, source=None, key='group', overlap='exclude', **filters)`

Gen gruplarını tek CRN listesine yazar; dosya `DifferentialCodonUsage.from_crn` ile doğrudan açılır. `groups` grup adından terime, terim listesine ya da gruba özel filtrelere giden sözlüktür (`{"term": "Ataxia", "min_frequency": 0.3}`). `file_name` yazılacak dosyadır. `source` `to_genes` ile aynıdır, `key` grup anahtarıdır. `overlap` birden çok gruba giren genler içindir: `"exclude"` hepsinden çıkarır (varsayılan, çünkü iki gruptaki gen farkı siler), `"first"` ilk grupta bırakır, `"error"` hata verir. `**filters` bütün gruplara uygulanır.

`file`, `release`, `groups` (grup → gen sayısı), `overlap` (çakışan semboller) ve `skipped` (atlanan gen → neden) içeren bir rapor döndürür. `overlap` geçersizse, gruplar çakışıyorsa (`"error"`) ya da bir grup boş kalırsa `ValueError` verir.

```python
rapor = hpo.pack({"a": {"gene": ["ATXN1", "ATXN3"]}, "b": {"gene": ["SCN1A", "ATXN3"]}}, "gruplar.crn", source=gc, key="group", overlap="exclude")
print(rapor["groups"], rapor["overlap"])  # -> {'a': 1, 'b': 1} ['ATXN3']
```

### Hpo.write_tsv

`statik: Hpo.write_tsv(rows, file_name)`

Sorgu sonucunu (`rows`) `file_name` dosyasına sekmeyle ayrılmış yazar.

```python
Hpo.write_tsv(hpo.genes(gene="ATXN1"), "atxn1.tsv")
print(open("atxn1.tsv").readline().split("\t")[1])  # -> gene_symbol
```

## 8. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_OBO`, `_GENE_PHENOTYPE`, `_GENE_DISEASE`, `_HPOA` | Ham dosya adları |
| `_FILES` | Bu dört dosyanın demeti |
| `_SYNONYM` | OBO `synonym` satırının düzenli ifadesi |
| `_DEFINITION` | OBO `def` satırının düzenli ifadesi |
| `_NO_SYMBOL` | Sembolsüz gen işareti (`-`) |
| `_FEATURE_KEYS` | `annotate`'in sıfırladığı features anahtarları |
| `Hpo.FEATURE_KEYS`, `Hpo.RELEASE_KEY`, `Hpo.GROUP_FILTER` | Taban sınıfın kancaları (bkz. [_geneset.md](_geneset.md#1-alt-sınıfın-sağladıkları)) |
| `_as_list` | Tek değeri listeye sarar, `None`'ı boş liste yapar |
| `Hpo.__repr__` | `Hpo('v2026-09-01', '<sys>/hpo')` |
| `Hpo._normalize_release` | `"2026-09-01"` → `"v2026-09-01"` |
| `Hpo._release_info` | GitHub'dan sürüm bilgisini alır; `User-Agent` gönderir, 404'te `ValueError`, öteki hatalarda `RuntimeError` verir |
| `Hpo._download` | `.part` dosyasına indirir ve SHA-256 doğrular |
| `Hpo._sha256` | Dosyanın özetini ilerleme çubuğuyla hesaplar |
| `Hpo._ensure` | Veri hazır değilse indirir; yerel sürüm varsa ağa çıkmaz |
| `Hpo._rows` | Yorum satırlarını atarak TSV okur |
| `Hpo._parse_ontology` | `hp.obo`'yu terimler, eş anlamlılar, alt/üst bağlantılar ve eski kimlikler olarak okur |
| `Hpo._load` | Dosyaları bir kez belleğe alır, `NOT` ve `0/N` satırlarını çıkarır, gen ve terim dizinlerini kurar |
| `Hpo._closure` | Tohum terimlerden bağlantıları izleyerek ulaşılan bütün terimleri bulur |
| `Hpo._resolve` | Kimlik, eski kimlik, ad ya da eş anlamlıdan güncel HP kimliği |
| `Hpo._terms` | Sorguyu terim kümesine çevirir, alt terimleri ekler |
| `Hpo._gene_ids` | Sembol ya da Gene ID listesini kimlik kümesine çevirir |
| `Hpo._disease_ids` | Hastalık sorgusunu kimlik kümesine çevirir |
| `Hpo._write_features` | Satırlardaki terim ve hastalıkları tekrarsız yazar |
| `Hpo._describe` | `_build` kancası: kurulan gene süzgeçten geçen HPO terim ve hastalıklarını yazar |

`Hpo`, `GeneSetSource` tabanından türer. `pack` ve `write_tsv` ile `_check_choices`, `_gencode`, `_canonical`, `_isoform`, `_build`, `_to_genes`, `_to_sequences` yardımcıları bu tabanda tanımlıdır ve `Clinvar` ile paylaşılır; ayrıntı [_geneset.md](_geneset.md) belgesindedir.
