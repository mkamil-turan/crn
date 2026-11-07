# `crn.fetch`: veri kaynakları

Kaynak: [`crn/fetch/__init__.py`](../../crn/fetch/__init__.py)

Dört veri kaynağı aynı mantıkla çalışır: veri yerel `sys/` klasöründe tutulur, bir şey istendiğinde önce oraya bakılır, yoksa indirilir. İkinci istekte ağa çıkılmaz. Klasör varsayılan olarak proje kökündedir; `sys_dir=` parametresiyle ya da `CRN_SYS_DIR` ortam değişkeniyle değiştirilir. İndirme ve hazırlık adımları tqdm ilerleme çubuğu gösterir, `verbose=False` bunları kapatır.

| Sınıf | Veri | İndirme | Belge |
|---|---|---|---|
| `Ncbi` | RefSeq: genomik DNA, tüm RNA'lar, proteinler | Gen başına, E-utilities | [ncbi.md](ncbi.md) |
| `NcbiSession` | HTTP oturumu: hız sınırı ve yeniden deneme | | [session.md](session.md) |
| `Gencode` | Ensembl/GENCODE: genom, transkriptler, proteinler, GTF | Sürüm başına, FTP (1,3 GB indirme, açılmış halde 5 GB) | [gencode.md](gencode.md) |
| `Hgnc` | Onaylı gen sembolleri, gen aileleri, 10 türün ortologları | Tablo olarak | [hgnc.md](hgnc.md) |
| `Hpo` | Semptomdan gene: ontoloji ve gen-fenotip notları | Sürüm başına, GitHub | [hpo.md](hpo.md) |

Ortak taban sınıf [`_store.md`](_store.md) içindedir. `Ncbi` ve `Gencode` aynı adlı metotlar sunar (`dna`, `rna`, `protein`, `gene`, `region`, `update`), bu yüzden kod bir kaynaktan ötekine az değişiklikle geçer. Hangisini seçeceğiniz sorunuza bağlıdır: tek bir gene bakıyorsanız `Ncbi` birkaç saniyede yeter, genom ölçeğinde çalışıyorsanız `Gencode` yerel dosyadan okuduğu için ağ beklemez ve NCBI hız sınırına takılmaz.

```python
from crn.fetch import Ncbi, NcbiSession, Gencode, Hgnc, Hpo

print(Gencode(verbose=False).release)  # -> v50
```
