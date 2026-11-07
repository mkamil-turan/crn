# `crn.fetch.session`: NCBI HTTP oturumu

Kaynak: [`crn/fetch/session.py`](../../crn/fetch/session.py)

NCBI, saniyedeki istek sayısını sınırlar ve aşan istemciyi 429 ile geri çevirir. `NcbiSession` bu sınırı kendisi uygular, geçici sunucu hatalarında yeniden dener ve bağlantıları havuzda tutar. `Ncbi` kendi oturumunu açar; toplu indirmede API anahtarı vermek için oturumu siz kurup geçirirsiniz.

<!-- toc -->
- [1. Oturum](#1-oturum): [`NcbiSession`](#ncbisession) · [`get`](#ncbisessionget) · [`is_open`](#ncbisessionis_open) · [`close`](#ncbisessionclose)
- [2. Hız sınırlayıcı](#2-hız-sınırlayıcı): [`_RateLimiter`](#_ratelimiter) · [`wait`](#_ratelimiterwait) · [`interval`](#_ratelimiterinterval)
<!-- /toc -->

## 1. Oturum

### NcbiSession

`NcbiSession(api_key=None, e_mail=None, timeout=30, max_retries=5, pool_connections=10, pool_maxsize=20)`

`api_key` verilmezse NCBI'nın anonim sınırına (saniyede 3 istek, `_RATE_NO_KEY`) uyulur; verilirse sınır saniyede 10'dur (`_RATE_WITH_KEY`). Anahtar her isteğe parametre olarak eklenir. `e_mail` User-Agent başlığına yazılır: NCBI sorun çıkan istemciye ulaşabilmek için bunu ister. `timeout` istek zaman aşımıdır (saniye). `max_retries` 429, 500, 502, 503 ve 504 yanıtlarında yeniden deneme sayısıdır, bekleme üstel artar (`backoff_factor=1.2`). `pool_connections` ve `pool_maxsize` bağlantı havuzunun boyutlarıdır.

Anahtarı kaynak koda yazmak yerine ortam değişkeninden okuyun.

```python
import os
from crn.fetch import NcbiSession, Ncbi

oturum = NcbiSession(api_key=os.environ.get("NCBI_API_KEY", "this_is_an_api_key"), timeout=20, max_retries=3)
print(oturum)  # -> NcbiSession(api_key=***_key, hız=10 istek/sn, durum=Açık)
print(oturum.is_open())  # -> True
oturum.close()
```

### NcbiSession.get

`NcbiSession.get(url, **kwargs)`

Hız sınırını uygulayarak `url` adresine GET gönderir. `**kwargs` `requests.Session.get`'e gider (ör. `params`, `headers`); `timeout` verilmezse oturumunki kullanılır. HTTP hatasında `requests.exceptions.RequestException`, kapalı oturumda `RuntimeError` verir.

<!-- calistirma -->
```python
with NcbiSession() as session:
    yanit = session.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                        params={"db": "gene", "term": "BDNF[sym] AND human[orgn]", "retmode": "json"})
    print(yanit.json()["esearchresult"]["idlist"])
```

### NcbiSession.is_open

`NcbiSession.is_open()`

Oturum açıksa `True`, `close` çağrıldıysa `False` verir.

```python
print(oturum.is_open())  # -> False
```

### NcbiSession.close

`NcbiSession.close()`

Oturumu kapatır. Sonraki `get` çağrısı `RuntimeError` verir. `with` bloğu bunu kendisi yapar: `NcbiSession.__enter__` oturumu döndürür, `NcbiSession.__exit__` kapatır.

```python
with NcbiSession() as s:
    icinde = s.is_open()
print(icinde, s.is_open())  # -> True False
```

`NcbiSession.__repr__` anahtarın yalnızca son dört karakterini gösterir, böylece oturum nesnesi günlüğe düşse bile anahtar sızmaz.

Bir oturum birden çok `Ncbi` nesnesine verilebilir; hız sınırı istekleri ortak saydığı için farklı organizmalarla çalışırken de sınır aşılmaz.

```python
with NcbiSession() as s:
    insan = Ncbi("human", session=s, verbose=False)
    fare = Ncbi("mouse", session=s, verbose=False)
print(insan.session is fare.session)  # -> True
```

## 2. Hız sınırlayıcı

### _RateLimiter

`_RateLimiter(min_interval)`

İki istek arasına en az `min_interval` saniye koyan iş parçacığı güvenli sınırlayıcıdır. Bir kilit ve son istek zamanını tutar, bu yüzden aynı oturumu paylaşan iş parçacıkları sınırı birlikte aşamaz.

### _RateLimiter.wait

`_RateLimiter.wait()`

Gerekirse bir sonraki istek zamanına kadar bekler.

### _RateLimiter.interval

`_RateLimiter.interval  (özellik)`

Geçerli bekleme aralığını saniye olarak verir; atama kilit altında yapılır.

```python
from crn.fetch.session import _RateLimiter

sinir = _RateLimiter(0.001)
sinir.wait()
print(sinir.interval)  # -> 0.001
```
