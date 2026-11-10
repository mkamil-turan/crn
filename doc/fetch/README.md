# `crn.fetch`: veri kaynakları

Kaynak: [`crn/fetch/__init__.py`](../../crn/fetch/__init__.py)

Veri kaynakları aynı mantıkla çalışır: veri yerel `sys/` klasöründe tutulur, bir şey istendiğinde önce oraya bakılır, yoksa indirilir. İkinci istekte ağa çıkılmaz. Klasör varsayılan olarak proje kökündedir; `sys_dir=` parametresiyle ya da `CRN_SYS_DIR` ortam değişkeniyle değiştirilir. İndirme ve hazırlık adımları tqdm ilerleme çubuğu gösterir, `verbose=False` bunları kapatır.

| Sınıf | Veri | İndirme | Belge |
|---|---|---|---|
| `Ncbi` | RefSeq: genomik DNA, tüm RNA'lar, proteinler; herhangi bir organizma (ad ya da taxid) | Gen başına, E-utilities | [ncbi.md](ncbi.md) |
| `NcbiSession` | HTTP oturumu: hız sınırı ve yeniden deneme | | [session.md](session.md) |
| `Gencode` | Ensembl/GENCODE: genom, transkriptler, proteinler, GTF; yalnızca insan (9606) ve fare (10090) | Sürüm başına, FTP (1,3 GB indirme, açılmış halde 5 GB) | [gencode.md](gencode.md) |
| `Hgnc` | Onaylı insan gen sembolleri, gen aileleri, HCOP ortologları (19 tür; tür ad ya da taxid ile) | Tablo olarak | [hgnc.md](hgnc.md) |
| `Hpo` | Semptomdan gene: ontoloji ve gen-fenotip notları | Sürüm başına, GitHub | [hpo.md](hpo.md) |
| `Clinvar` | Klinik varyantlardan gene ve hastalığa | Arşiv ayı başına, FTP (450 MB indirme, GRCh38 için 1 GB SQLite) | [clinvar.md](clinvar.md) |

Ortak taban sınıf [`_store.md`](_store.md) içindedir; gen listesi veren `Hpo` ve `Clinvar` ayrıca [`_geneset.md`](_geneset.md) tabanını paylaşır. `Ncbi` ve `Gencode` aynı adlı metotlar sunar (`dna`, `rna`, `protein`, `gene`, `region`, `update`), bu yüzden kod bir kaynaktan ötekine az değişiklikle geçer. Hangisini seçeceğiniz sorunuza bağlıdır: tek bir gene bakıyorsanız `Ncbi` birkaç saniyede yeter, genom ölçeğinde çalışıyorsanız `Gencode` yerel dosyadan okuduğu için ağ beklemez ve NCBI hız sınırına takılmaz.

```python
from crn.fetch import Ncbi, NcbiSession, Gencode, Hgnc, Hpo, Clinvar

print(Gencode(verbose=False).release)  # -> v50
```
