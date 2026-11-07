# `crn.fetch._store`: ortak yerel depo

Kaynak: [`crn/fetch/_store.py`](../../crn/fetch/_store.py)

Veri kaynaklarının paylaştığı küçük yardımcılar. Modül içe aktarılarak değil, `Ncbi`, `Gencode`, `Hgnc` ve `Hpo` üzerinden kullanılır; ancak yeni bir veri kaynağı eklemek isteyen burada başlar.

<!-- toc -->
- [1. Depo konumu](#1-depo-konumu)
- [2. Dosya yardımcıları](#2-dosya-yardımcıları): [`atomic_write`](#atomic_write) · [`split_fasta`](#split_fasta) · [`read_single_fasta`](#read_single_fasta)
- [3. Taban sınıf](#3-taban-sınıf): [`LocalSource`](#localsource)
- [4. İç yardımcılar](#4-iç-yardımcılar)
<!-- /toc -->

## 1. Depo konumu

`DEFAULT_SYS_DIR` varsayılan veri klasörüdür: `CRN_SYS_DIR` ortam değişkeni tanımlıysa o, değilse proje kökündeki `sys/`. Paket içe aktarılırken bir kez okunur, bu yüzden değişken Python başlamadan önce ayarlanmalıdır.

## 2. Dosya yardımcıları

### atomic_write

`atomic_write(path, text)`

`text` metnini önce `path` + `.tmp` olarak yazar, sonra yerine taşır. İndirme yarıda kesilse bile yarım dosya kalmaz; `sys/` içindeki bir dosyanın varlığı tamamlanmış olduğunun kanıtıdır.

```python
from pathlib import Path
from crn.fetch._store import atomic_write, split_fasta, read_single_fasta, LocalSource

atomic_write(Path("a.fa"), ">kayit bir\nATGC\nGG\n")
print(Path("a.fa").read_text().count(">"))  # -> 1
```

### split_fasta

`split_fasta(text)`

Çok kayıtlı FASTA metnini (`text`) (ilk kelime, kaydın metni) çiftlerine böler. NCBI yanıtlarını kayıt başına dosyalara dağıtmak için kullanılır.

```python
print([k for k, _ in split_fasta(">a x\nAT\n>b\nGC\n")])  # -> ['a', 'b']
```

### read_single_fasta

`read_single_fasta(path)`

`path` içindeki tek kayıtlı FASTA'dan (başlık, dizi) verir, satırları birleştirir.

```python
print(read_single_fasta(Path("a.fa")))  # -> ('kayit bir', 'ATGCGG')
```

## 3. Taban sınıf

### LocalSource

`LocalSource(sys_dir=None, verbose=True)`

Yerel klasörle çalışan veri kaynağının tabanıdır. `sys_dir` verilmezse `DEFAULT_SYS_DIR`, `verbose` ilerleme çubuklarını açar. Alt sınıf `SOURCE` adını tanımlar; çıkarılan her nesnenin features'ına `source` olarak yazılır.

```python
class Kaynak(LocalSource):
    SOURCE = "ornek"

k = Kaynak(sys_dir="depo", verbose=False)
print(k.sys_dir, k.SOURCE)  # -> depo ornek
```

## 4. İç yardımcılar

| Ad | Görevi |
|---|---|
| `LocalSource._log` | `verbose` açıksa ilerleme çubuklarını bozmadan yazar |
| `LocalSource._bar` | tqdm çubuğu; `unit="B"` ise bayt ölçeğiyle, `verbose=False` ise gizli |
| `LocalSource._progress` | Bir öğe listesini tqdm ile sarar |
| `LocalSource._features` | `source` ile boş olmayan değerleri nesnenin features'ına yazar; `empty` hangi değerlerin boş sayılacağını belirler |
