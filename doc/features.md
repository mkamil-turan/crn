# `crn.features`: ek bilgi deposu

Kaynak: [`crn/features.py`](../crn/features.py)

Kaynak, gen adı, CDS koordinatı gibi bilgiler `Sequence` ve `Gene` sınıflarına yeni alan olarak eklenmez, `Features` içine yazılır. Böylece veri kaynakları sınıfları değiştirmeden istedikleri bilgiyi taşır ve bu bilgi CRN dosyasına kayıpsız geçer. Her anahtar bir değer listesi tutar, çünkü bir gene birden çok hastalık ya da HPO terimi bağlanabilir.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`Features`](#features) · [`add`](#featuresadd) · [`set`](#featuresset) · [`get`](#featuresget) · [`first`](#featuresfirst) · [`remove`](#featuresremove) · [`clear`](#featuresclear) · [`keys`](#featureskeys) · [`items`](#featuresitems) · [`to_dict`](#featuresto_dict) · [`copy`](#featurescopy)
- [2. Python protokolleri](#2-python-protokolleri)
<!-- /toc -->

## 1. Sınıf

### Features

`Features()`

Boş bir depo oluşturur; anahtarlar eklenme sırasını korur. Dışarı verilen listeler kopyadır, yani dönen listeyi değiştirmek depoyu değiştirmez. Anahtar boş olmayan bir `str` olmalıdır, değilse `TypeError` ya da `ValueError` verilir.

```python
from crn import Features

f = Features()
f.add("gene", "BRCA1")
f.add("gene", "BRCA2")
print(f.get("gene"))  # -> ['BRCA1', 'BRCA2']
print(f)  # -> Keys:gene
```

### Features.add

`Features.add(key, value)`

`key` anahtarının listesinin sonuna `value` ekler. Anahtar yoksa oluşturur.

```python
f.add("cds", (368, 2173))
print(f.first("cds"))  # -> (368, 2173)
```

### Features.set

`Features.set(key, value)`

`key` anahtarının bütün değerlerini silip yerine tek `value` koyar. `add` listeyi büyütür, `set` değiştirir; gen grubu gibi tek değerli bilgiler için `set` kullanılır.

```python
f.set("group", "hasta")
f.set("group", "kontrol")
print(f.get("group"))  # -> ['kontrol']
```

### Features.get

`Features.get(key, default=None)`

`key` anahtarının değer listesinin kopyasını verir; anahtar yoksa `default` döner.

```python
print(f.get("yok", default=[]))  # -> []
```

### Features.first

`Features.first(key, default=None)`

Listenin ilk değerini verir; anahtar yoksa `default` döner. Tek değerli anahtarlarda en sık kullanılan okuma yoludur. `key` aranan anahtardır.

```python
print(f.first("gene"), f.first("yok", default="-"))  # -> BRCA1 -
```

### Features.remove

`Features.remove(key)`

Anahtarı bütün değerleriyle siler. `key` varsa `True`, yoksa `False` döner.

```python
print(f.remove("cds"), f.remove("cds"))  # -> True False
```

### Features.clear

`Features.clear()`

Bütün anahtarları siler.

```python
g = Features(); g.add("a", 1); g.clear()
print(len(g))  # -> 0
```

### Features.keys

`Features.keys()`

Anahtarları eklenme sırasıyla liste olarak verir.

```python
print(f.keys())  # -> ['gene', 'group']
```

### Features.items

`Features.items()`

`(anahtar, değer listesi)` çiftlerini üretir. Listeler kopyadır.

```python
print(list(f.items()))  # -> [('gene', ['BRCA1', 'BRCA2']), ('group', ['kontrol'])]
```

### Features.to_dict

`Features.to_dict()`

Depoyu `{anahtar: liste}` sözlüğü olarak verir. Sözlük de listeler de kopyadır.

```python
print(f.to_dict())  # -> {'gene': ['BRCA1', 'BRCA2'], 'group': ['kontrol']}
```

### Features.copy

`Features.copy()`

Derin kopya döndürür. Dilimleme ve dönüşüm metotları yeni nesnelere bu kopyayı verir, bu yüzden türetilmiş dizinin features'ını değiştirmek özgün diziyi etkilemez.

```python
kopya = f.copy()
kopya.add("gene", "TP53")
print(len(f.get("gene")), len(kopya.get("gene")))  # -> 2 3
```

## 2. Python protokolleri

| Ad | Davranış |
|---|---|
| `Features.__contains__` | `"gene" in f` |
| `Features.__len__` | Anahtar sayısı |
| `Features.__iter__` | Anahtarlar üzerinde gezinir; dolaşırken silme güvenlidir, çünkü liste kopyasını dolaşır |
| `Features.__repr__` | `Keys:gene,group` |
| `Features._check_key` | Anahtarın boş olmayan `str` olduğunu denetleyen iç yardımcı |
| `Features.__slots__` | Nesneyi tek bir `_data` alanıyla sınırlar, örnek başına `__dict__` tutulmaz |
