# crn_fetch_ncbi

`crn/fetch/crn_fetch_ncbi.py` modülü, NCBI'dan (RefSeq) genomik DNA, transkript
(RNA), protein ve bunların hepsini bir arada tutan `Gene` nesnelerini çıkaran
**`Ncbi`** sınıfını tanımlar. İstekler `crn/fetch/crn_fetch_session.py`
içindeki **`NcbiSession`** üzerinden yapılır.

Örnek betik: [test/ncbi_test.py](../test/ncbi_test.py)
Genel kılavuz: [crn_kullanim_kilavuzu.md](crn_kullanim_kilavuzu.md) (bölüm 11 ve 13)

Bu belgedeki tüm çıktılar MEFV geni (NCBI Gene ID `4210`) üzerinde
çalıştırılarak alınmıştır.

## İçindekiler

1. [Hızlı başlangıç](#1-hızlı-başlangıç)
2. [Çalışma mantığı](#2-çalışma-mantığı)
3. [`sys/ncbi/` klasörü](#3-sysncbi-klasörü)
4. [`NcbiSession` — HTTP oturumu](#4-ncbisession--http-oturumu)
5. [`Ncbi` oluşturmak](#5-ncbi-oluşturmak)
6. [`search_gene` — sembolden Gene ID'ye](#6-search_gene--sembolden-gene-idye)
7. [`update` — indirme](#7-update--indirme)
8. [`gene` — tüm içerikle Gene](#8-gene--tüm-içerikle-gene)
9. [`dna` — genomik DNA](#9-dna--genomik-dna)
10. [`rna` — transkript](#10-rna--transkript)
11. [`protein` — protein](#11-protein--protein)
12. [`region` — herhangi bir genomik bölge](#12-region--herhangi-bir-genomik-bölge)
13. [Gene ile çalışmak](#13-gene-ile-çalışmak)
14. [Kaydetme: CRN ve FASTA](#14-kaydetme-crn-ve-fasta)
15. [Tarifler](#15-tarifler)
16. [Hatalar](#16-hatalar)
17. [Features anahtarları özeti](#17-features-anahtarları-özeti)
18. [Sık sorulanlar ve ipuçları](#18-sık-sorulanlar-ve-ipuçları)

---

## 1. Hızlı başlangıç

```python
from crn.fetch import NcbiSession, Ncbi

session = NcbiSession(api_key="...")      # anahtar isteğe bağlı (bkz. bölüm 4)
ncbi = Ncbi(session=session, organism="human")

ncbi.search_gene("MEFV")        # '4210'
mefv = ncbi.gene("MEFV")        # ya da ncbi.gene("4210")
print(mefv)                     # MEFV: DNA <14607>, 2 transkript, 2 protein
```

Betikleri projenin sanal ortamıyla çalıştırın (`tqdm` ve `requests` sistem
Python'unda kurulu olmayabilir):

```bash
.venv/bin/python test/ncbi_test.py
```

## 2. Çalışma mantığı

```text
"MEFV" ──search_gene──▶ "4210" ──update──▶ sys/ncbi/ ──gene/dna/rna/protein──▶ DNA, RNA, PRT, Gene
           (ağ, bir kez)            (ağ)        (yerel)            (yerel okuma)
```

1. Gen sembolü NCBI Gene ID'sine çevrilir. Eşleşme `sys/` içinde saklandığı
   için aynı sembol için bir daha ağa çıkılmaz.
2. Genin kaydı (`sys/ncbi/genes/<gene_id>.json`) ve bu kaydın andığı tüm
   dosyalar varsa **ağa çıkılmaz**; her şey diskten okunur.
3. Bir dosya eksikse `update()` kendiliğinden çağrılır.

`update()`'in indirdikleri:

- gen özeti (ad, açıklama, kromozom, sitogenetik konum, takma adlar);
- genin genomik bölgesi;
- gene bağlı tüm RefSeq RNA'lar (NM_, NR_, varsa XM_, XR_) ve özellik
  tabloları;
- bu RNA'ların kodladığı proteinler (NP_, XP_).

Yani `update()`'i elle çağırmak **zorunlu değildir**. İlk `gene()`, `dna()`
çağrısı onu kendisi yapar.

## 3. `sys/ncbi/` klasörü

```text
sys/ncbi/
    symbols/<organizma>.json        sembol -> Gene ID           (ör. {"MEFV": "4210"})
    genes/<gene_id>.json            özet + RNA/protein listesi  (update'in son adımı)
    nuccore/<erişim>.fa             transkript dizisi           (NM_000243.3.fa)
    nuccore/<erişim>.ft             transkriptin özellik tablosu (CDS, protein_id, GeneID)
    protein/<erişim>.fa             protein dizisi              (NP_000234.1.fa)
    region/<erişim>_<baş>_<son>_<plus|minus>.fa
                                    genomik bölge               (NC_000016.10_3242027_3256633_minus.fa)
```

MEFV için `genes/4210.json`:

```json
{
 "gene_id": "4210",
 "gene_name": "MEFV",
 "description": "MEFV innate immunity regulator, pyrin",
 "organism": "Homo sapiens",
 "chromosome": "16",
 "map_location": "16p13.3",
 "aliases": ["FMF", "MEF", "PAAND", "TRIM20"],
 "region": ["NC_000016.10", 3242027, 3256633, "-"],
 "transcripts": ["NM_000243.3", "NM_001198536.2"],
 "proteins": {"NP_000234.1": "NM_000243.3", "NP_001185465.2": "NM_001198536.2"},
 "updated": "..."
}
```

- Dosyalar önce `.tmp` olarak yazılıp sonra yerine konur; yarıda kesilen bir
  indirme bozuk dosya bırakmaz.
- `genes/<id>.json` en son yazılır. Bu dosya yoksa gen "hazır" sayılmaz ve bir
  sonraki istekte yeniden indirilir.
- Klasörü silmek güvenlidir; gerekenler ilk istekte yeniden indirilir.

## 4. `NcbiSession` — HTTP oturumu

```python
NcbiSession(api_key=None, e_mail=None, timeout=30, max_retries=5,
            pool_connections=10, pool_maxsize=20)
```

| Parametre | Açıklama |
|---|---|
| `api_key` | NCBI API anahtarı. Yoksa 3 istek/sn, varsa 10 istek/sn. |
| `e_mail` | NCBI'nın istediği iletişim adresi; `User-Agent` başlığına yazılır. |
| `timeout` | Her isteğin zaman aşımı (saniye). |
| `max_retries` | 429, 500, 502, 503 ve 504 yanıtlarında yeniden deneme sayısı (artan beklemeyle). |
| `pool_connections`, `pool_maxsize` | Bağlantı havuzu ayarları. |

```python
>>> NcbiSession()
NcbiSession(Anonim, hız=3 istek/sn, durum=Açık)
>>> NcbiSession(api_key="...")
NcbiSession(api_key=***abcd, hız=10 istek/sn, durum=Açık)
```

### Kullanım biçimleri

```python
# 1) Oturumu Ncbi kendisi açsın (anonim, 3 istek/sn)
ncbi = Ncbi()

# 2) Kendi oturumunuzu verin
session = NcbiSession(api_key="...", e_mail="ben@ornek.com")
ncbi = Ncbi(session=session)
...
session.close()

# 3) Context manager (önerilen): blok bitince oturum kapanır
with NcbiSession(api_key="...") as session:
    ncbi = Ncbi(session=session)
    mefv = ncbi.gene("MEFV")
```

> Blok bittikten sonra `ncbi` ağa çıkmak isterse `RuntimeError: Kapalı bir
> oturum üzerinden istek yapılamaz.` alırsınız. Diskte hazır olan kayıtlar
> yine okunabilir.

### API anahtarını koda yazmamak

Anahtarı betiğe gömmek yerine ortam değişkeninden okuyun:

```python
import os
session = NcbiSession(api_key=os.environ.get("NCBI_API_KEY"))
```

```bash
set -x NCBI_API_KEY "anahtarınız"     # fish
```

Anahtar [NCBI hesap ayarlarından](https://account.ncbi.nlm.nih.gov/settings/)
ücretsiz alınır.

## 5. `Ncbi` oluşturmak

```python
Ncbi(organism="human", session=None, sys_dir=None, verbose=True)
```

| Parametre | Açıklama |
|---|---|
| `organism` | `human`, `mouse`, `rat`, `zebrafish`, `fruit_fly`, `worm`. Başka bir değer `ValueError` verir. **Yalnızca sembol aramasını** etkiler; Gene ID ya da erişimle yapılan isteklerde önemi yoktur. |
| `session` | `NcbiSession`; verilmezse anonim bir oturum açılır. |
| `sys_dir` | Veri klasörü; varsayılan proje kökündeki `sys/`. Veriler `<sys_dir>/ncbi/` altına yazılır. |
| `verbose` | `True` ise indirme ve okuma adımları `tqdm` ilerleme çubuklarıyla gösterilir. |

```python
ncbi = Ncbi()                                   # insan, varsayılan klasör
ncbi = Ncbi(organism="mouse")                   # fare sembolleri
ncbi = Ncbi(sys_dir="/veri/crn_sys")            # başka bir klasör
ncbi = Ncbi(verbose=False)                      # sessiz (betikler, testler)
```

Özellikler:

| Özellik | Açıklama |
|---|---|
| `ncbi.organism` | Seçili organizma. |
| `ncbi.session` | Kullanılan `NcbiSession`. |
| `ncbi.sys_dir` | Veri kökü (`Path`). |
| `ncbi.ncbi_dir` | `sys_dir / "ncbi"`. |
| `ncbi.verbose` | İlerleme çubukları açık mı (sonradan değiştirilebilir). |

## 6. `search_gene` — sembolden Gene ID'ye

```python
search_gene(gene) -> str
```

Gen sembolünü NCBI Gene ID'sine çevirir. **Kayıt indirmez**, yalnızca kimliği
bulur.

- Rakamlardan oluşan değer zaten Gene ID kabul edilip olduğu gibi döner
  (`int` de verilebilir).
- Arama `<sembol>[Gene Name] AND <taxon> AND alive[prop]` biçimindedir:
  yalnızca seçili organizmadaki, yürürlükteki (silinmemiş) kayıtlar aranır.
- Bulunan eşleşme `symbols/<organizma>.json` dosyasına yazılır.

```python
>>> ncbi.search_gene("MEFV")
'4210'
>>> ncbi.search_gene(4210)
'4210'
>>> Ncbi(organism="mouse").search_gene("Mefv")
'54483'
```

> Diğer tüm metotlar içeride `search_gene`'i çağırdığından, onlara da sembol ya
> da Gene ID verebilirsiniz. Bu adımı ayrıca yapmanız gerekmez.

## 7. `update` — indirme

```python
update(gene, force=False) -> str    # Gene ID döner
```

Genin tüm kayıtlarını indirip `sys/` içine yazar. Dört adımda ilerler (ilerleme
çubuğunda görünür): arama ve özet, RNA'lar, proteinler, genomik bölge.

```python
ncbi.update("MEFV")               # eksikleri indir
ncbi.update("4210", force=True)   # hepsini NCBI'daki güncel hâliyle yenile
```

| Kayıt | `force=False` iken |
|---|---|
| Gen özeti ve RNA listesi | Her çağrıda NCBI'dan alınır. |
| RNA FASTA | Her çağrıda indirilir (NCBI sayısal UID döndürdüğü için hangisinin diskte olduğu önceden bilinemez), ama diskte olan dosyanın üzerine yazılmaz. |
| RNA özellik tablosu, protein, genomik bölge | Yalnızca eksik olanlar indirilir. |
| `genes/<id>.json` | Her çağrıda yeniden yazılır. |

**Ne zaman elle çağırmalı?**

- Ağ erişimi olan bir anda verileri önceden hazırlamak için (ör. birçok gen için
  toplu indirme, bkz. [bölüm 15](#15-tarifler)).
- NCBI'da yeni bir sürüm ya da yeni bir transkript çıktığını düşünüyorsanız
  (`force=True`).

Sadece okumak istiyorsanız `update` yerine doğrudan `gene()` kullanın; kayıt
hazırsa ağa hiç çıkmaz.

## 8. `gene` — tüm içerikle Gene

```python
gene(gene) -> Gene
```

Gene ait her şeyi tek bir `Gene` nesnesinde toplar:

```text
Gene (MEFV)
 ├── dna          : DNA            genomik dizi (genin zincirinde, 5'->3')
 ├── transcripts  : {id: RNA}      NM_000243.3, NM_001198536.2
 ├── proteins     : {id: PRT}      NP_000234.1, NP_001185465.2
 └── features     : Features       gene_id, gene_name, description, organism,
                                   chromosome, map_location, aliases, genomic_region
```

- Transkriptler önce doğrulanmış (NM_, NR_), sonra tahmini (XM_, XR_) olacak
  şekilde sıralanır. Kimlikleri sürümlü erişimleridir.
- Her RNA'nın `genomic_strand` alanı genin zincirinden doldurulur (MEFV için
  `-`).
- Proteinlerin kimlikleri erişimleridir. Her proteinin `transcript_id`
  özelliği, onu kodlayan transkripti gösterir.
- Genomik koordinatı olmayan genlerde `dna` `None` olur.

```python
>>> mefv = ncbi.gene("MEFV")
>>> mefv
MEFV: DNA <14607>, 2 transkript, 2 protein
>>> mefv.features.to_dict()
{'source': ['NCBI'], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
 'description': ['MEFV innate immunity regulator, pyrin'], 'organism': ['Homo sapiens'],
 'chromosome': ['16'], 'map_location': ['16p13.3'],
 'aliases': ['FMF', 'MEF', 'PAAND', 'TRIM20'],
 'genomic_region': [('NC_000016.10', 3242027, 3256633, '-')]}
```

Gene nesnesinin içinde gezinmek için bkz. [bölüm 13](#13-gene-ile-çalışmak).

## 9. `dna` — genomik DNA

```python
dna(gene) -> DNA
```

Genin genomik DNA'sını **genin kendi (kodlayan) zincirinde, 5'→3'**
döndürür. Gen (−) zincirdeyse dizi zaten ters tümleyicidir. Yani dizinin başı
genin 5' ucudur.

```python
>>> d = ncbi.dna("MEFV")
>>> d
MEFV: CTACCAGAAGCCAGACAGCTGGCTCGAGCCTCTCCTGCTC...AATTGACTGGGCCTCATGAATAAAAAAATTAGAAAATCTG <14607> [DNA] 5'->3'
>>> d.features.to_dict()
{'source': ['NCBI'], 'chromosome': ['NC_000016.10'], 'start': [3242027], 'end': [3256633],
 'chromosome_strand': ['-'], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
 'description': ['MEFV innate immunity regulator, pyrin']}
```

`start`/`end` referans kromozom üzerinde 1 tabanlı ve iki ucu dahil
koordinatlardır. Bunlar NCBI özetindeki 0 tabanlı değerlerden çevrilmiştir.

## 10. `rna` — transkript

```python
rna(accession, genomic_strand=GenomicStrand.Unknown) -> RNA
```

RefSeq transkriptini RNA olarak döndürür (`sys/` içinde yoksa indirir). Gene
bağlı olmak zorunda değildir; herhangi bir NM_/NR_/XM_/XR_ erişimi verilebilir.

| `accession` | Davranış |
|---|---|
| `"NM_000243.3"` | O sürüm. |
| `"NM_000243"` | Diskteki en yeni sürüm; diskte yoksa NCBI'daki güncel sürüm. |

```python
>>> r = ncbi.rna("NM_000243")
>>> r
NM_000243.3: CUACCAGAAGCCAGACAGCUGGCUCGAGCCUCUCCUGCUC...AAUUGACUGGGCCUCAUGAAUAAAAAAAUUAGAAAAUCUG <3506> [RNA] 5'->3' (?)
>>> r.features.to_dict()
{'source': ['NCBI'], 'transcript_id': ['NM_000243.3'],
 'title': ['Homo sapiens MEFV innate immunity regulator, pyrin (MEFV), transcript variant 1, mRNA'],
 'predicted': [False], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
 'protein_id': ['NP_000234.1'],
 'utr5': [(1, 46)], 'cds': [(47, 2392)], 'utr3': [(2393, 3506)]}
```

### Genomik zincir

Tek başına çağrıldığında RNA'nın hangi zincirden geldiği bilinmez ve `(?)`
görünür. Bunu biliyorsanız kendiniz verin, ya da RNA'yı `gene()` üzerinden alın:

```python
>>> from crn.io import GenomicStrand
>>> ncbi.rna("NM_000243.3", GenomicStrand.Minus)
NM_000243.3: CUACC...AUCUG <3506> [RNA] 5'->3' (-)
```

### CDS ve UTR'ler

Protein kodlayan transkriptlerde `utr5`, `cds`, `utr3` özellikleri
`(başlangıç, bitiş)` biçimindedir: transkript üzerinde **1 tabanlı, iki ucu
dahil**. Python dilimine çevirirken başlangıçtan 1 çıkarın:

```python
start, end = r.features.first("cds")    # (47, 2392)
cds = r[start - 1:end]                   # RNA, 2346 nt, AUG ile başlar
utr5 = r[:start - 1]                     # 46 nt

if "utr3" in r.features:
    s, e = r.features.first("utr3")
    utr3 = r[s - 1:e]
```

Kodlamayan RNA'larda (NR_, XR_) `cds`, `utr5`, `utr3` ve `protein_id`
bulunmaz.

`predicted` özelliği XM_/XR_ (tahmini, hesaplamalı) kayıtlarda `True` olur.

## 11. `protein` — protein

```python
protein(accession) -> PRT
```

| `accession` | Davranış |
|---|---|
| `"NP_000234.1"`, `"NP_000234"` | Protein doğrudan (sürümsüzde en yeni). |
| `"NM_000243"` (ya da XM_) | Önce transkript okunur, onun `protein_id`'si alınır. Proteinin features'ına `transcript_id` de yazılır. |
| `"NR_..."` (ya da XR_) | Kodlamayan RNA: `KeyError`. |

```python
>>> ncbi.protein("NP_000234").features.to_dict()
{'source': ['NCBI'], 'protein_id': ['NP_000234.1'], 'title': ['pyrin isoform 1 [Homo sapiens]']}

>>> p = ncbi.protein("NM_000243")
>>> p
NP_000234.1: MAKTPSDHLLSTLEELVPYDFEKFKFKLQNTSVQKEHSRI...TFASCSFSGPLQPIFSPGTRDGGKNTAPLTICPVGGQGPD <781> [PRT] N->C
>>> p.features.to_dict()
{'source': ['NCBI'], 'protein_id': ['NP_000234.1'], 'transcript_id': ['NM_000243.3'],
 'title': ['pyrin isoform 1 [Homo sapiens]']}

>>> ncbi.protein("NR_002819.5")
KeyError: 'NCBI: NR_002819.5 protein kodlamıyor'
```

## 12. `region` — herhangi bir genomik bölge

```python
region(accession, start, end, strand="+", name="") -> DNA
```

Bir kromozomdan (ya da contig'den) istenen bölgeyi indirir. Gen bilgisi
gerekmez.

| Parametre | Açıklama |
|---|---|
| `accession` | Kromozom erişimi, ör. `NC_000016.10` (GRCh38 kromozom 16). Kromozom adı (`chr16`) değil. |
| `start`, `end` | 1 tabanlı, iki ucu dahil. `1 <= start <= end` olmalı. |
| `strand` | `"+"` referans zincir; `"-"` ters tümleyici. İkisinde de dizi 5'→3' yazılır. |
| `name` | DNA'nın adı; boşsa `NC_000016.10:3242028-3242100(-)` biçiminde olur. |

```python
>>> b = ncbi.region("NC_000016.10", 3242028, 3242100, "-")
>>> b
NC_000016.10:3242028-3242100(-): ATGAGCCACCATGCCGGGCCATGATTCTTAAGAGAATTGACTGGGCCTCATGAATAAAAAAATTAGAAAATCT <73> [DNA] 5'->3'
>>> b.features.to_dict()
{'source': ['NCBI'], 'chromosome': ['NC_000016.10'], 'start': [3242028], 'end': [3242100],
 'chromosome_strand': ['-']}

>>> ncbi.region("NC_000016.10", 3242028, 3242100, "-", name="MEFV_3uc")   # özel ad
```

Her farklı `(erişim, başlangıç, bitiş, zincir)` dörtlüsü ayrı bir dosya olarak
saklanır. Aynı bölge ikinci kez istendiğinde diskten okunur.

**Genin çevresini almak** (ör. promotör için 2 kb yukarı akış):

```python
acc, start, end, strand = ncbi.gene("MEFV").features.first("genomic_region")
if strand == "-":           # (-) zincirde genin 5' ucu 'end' tarafındadır
    promoter = ncbi.region(acc, end + 1, end + 2000, "-", name="MEFV_promoter")
else:
    promoter = ncbi.region(acc, start - 2000, start - 1, "+", name="MEFV_promoter")
```

## 13. Gene ile çalışmak

`ncbi.gene()` sonucu sıradan bir `crn.io.Gene` nesnesidir.

```python
mefv = ncbi.gene("MEFV")

mefv.name                       # 'MEFV'
mefv.dna                        # DNA <14607>
mefv.features.first("map_location")      # '16p13.3'
mefv.features.get("aliases")             # ['FMF', 'MEF', 'PAAND', 'TRIM20']

# Transkriptler: {kimlik: RNA} (kopya sözlük)
for tid, rna in mefv.transcripts.items():
    start, end = rna.features.first("cds")
    print(tid, len(rna), (start, end), rna.features.first("protein_id"))
# NM_000243.3    3506 (47, 2392) NP_000234.1
# NM_001198536.2 3048 (47, 1384) NP_001185465.2

# Proteinler: {kimlik: PRT}
for pid, prt in mefv.proteins.items():
    print(pid, len(prt), prt.features.first("transcript_id"))
# NP_000234.1    781 NM_000243.3
# NP_001185465.2 445 NM_001198536.2

# Tek tek erişim (kimlikler SÜRÜMLÜ olmalı)
rna = mefv.transcript("NM_000243.3")
prt = mefv.protein("NP_000234.1")
```

### RNA ↔ protein eşleştirme

`Gene` sınıfı ikisi arasında bağ tutmaz. Bağı features sağlar:

```python
# Transkriptten proteine
prt = mefv.protein(rna.features.first("protein_id"))

# Proteinden transkripte
rna = mefv.transcript(prt.features.first("transcript_id"))

# Tüm çiftler
pairs = {tid: r.features.first("protein_id")
         for tid, r in mefv.transcripts.items()
         if "protein_id" in r.features}
```

### Gene'i düzenlemek

```python
mefv.remove_transcript("NM_001198536.2")     # True
mefv.remove_protein("NP_001185465.2")        # True
mefv.features.add("note", "yalnızca kanonik izoform")
mefv.name = "MEFV_canonical"
```

`add_transcript` / `add_protein` ile başka bir kaynaktan gelen dizileri de
ekleyebilirsiniz. Aynı kimlik ikinci kez eklenirse `ValueError` verir.

## 14. Kaydetme: CRN ve FASTA

### Tüm geni tek dosyada (CRN)

```python
from crn.io import Gene

mefv.to_crn("MEFV.crn")
mefv2 = Gene.from_crn("MEFV.crn")    # ağ ya da sys/ gerekmez
print(mefv2)                         # MEFV: DNA <14607>, 2 transkript, 2 protein
```

Ad, features, DNA, tüm transkript ve proteinler, kimlikleri ve features'larıyla
birlikte yazılır. Biçim pickle değildir; kendi ikili biçimidir.

### Tek bir diziyi FASTA'ya

```python
from crn.io import PRT

mefv.protein("NP_000234.1").to_fasta("NP_000234.1.fa")
mefv.dna.to_fasta("MEFV_genomic.fa")

prts = PRT.from_fasta_file("NP_000234.1.fa")     # liste döner
```

## 15. Tarifler

### CDS'i çevirip NCBI proteiniyle karşılaştırmak

```python
for tid, rna in mefv.transcripts.items():
    start, end = rna.features.first("cds")
    translated = rna[start - 1:end].translation()          # sonu '*' (stop)
    ncbi_prt = mefv.protein(rna.features.first("protein_id"))
    print(tid, translated.seq.rstrip("*") == ncbi_prt.seq)
# NM_000243.3 True
# NM_001198536.2 True
```

### Birçok geni toplu indirmek

```python
import os
from crn.fetch import NcbiSession, Ncbi

genes = ["MEFV", "NLRP3", "TNFRSF1A", "MVK"]

with NcbiSession(api_key=os.environ.get("NCBI_API_KEY")) as session:
    ncbi = Ncbi(session=session, verbose=False)
    for symbol in genes:
        try:
            gene_id = ncbi.update(symbol)
            print(symbol, gene_id, "tamam")
        except KeyError as e:
            print(symbol, "bulunamadı:", e)

# Daha sonra, ağ olmadan:
ncbi = Ncbi(verbose=False)
collection = {s: ncbi.gene(s) for s in genes}
```

### Yalnızca kanonik (ilk NM_) izoformu almak

```python
nm = [t for t in mefv.transcripts if t.startswith("NM_")]
canonical = mefv.transcript(nm[0])
```

> NCBI'nın "MANE Select" bilgisi bu sınıf tarafından çekilmez. Buradaki "ilk
> NM_" yalnızca sıralamaya dayanır; klinik olarak kanonik izoform olduğu
> garanti değildir.

### Başka bir organizma

```python
mouse = Ncbi(organism="mouse")
mouse.gene("Mefv")          # fare sembolleri büyük-küçük harf düzenine dikkat
```

Gene ID ya da erişimle çalışırken `organism` önemsizdir:
`Ncbi().gene("54483")` de fare Mefv'i getirir.

### Verileri farklı bir klasörde tutmak

```python
ncbi = Ncbi(sys_dir="/mnt/veri/crn_sys")     # /mnt/veri/crn_sys/ncbi/... altına yazar
```

## 16. Hatalar

| Hata | Ne zaman |
|---|---|
| `ValueError: Desteklenmeyen organizma: dog; seçenekler: ...` | `Ncbi(organism=...)` listede değilse. |
| `KeyError: "NCBI: 'X' geni human organizmasında bulunamadı"` | Sembol o organizmada yoksa (ya da yanlış organizma seçiliyse). |
| `KeyError: 'NCBI: Gene ID 123 bulunamadı'` | Gene ID geçersiz ya da silinmişse. |
| `KeyError: 'NCBI: Gene ID ... için genomik koordinat yok'` | `dna()`: genin kromozom konumu yoksa. |
| `KeyError: "NCBI: 'NM_...' transkripti bulunamadı"` | `rna()`: erişim NCBI'da yoksa. |
| `KeyError: 'NCBI: NR_... protein kodlamıyor'` | `protein()`'e kodlamayan RNA verildiyse. |
| `KeyError: "NCBI: '...' proteini bulunamadı"` | `protein()`: protein erişimi yoksa. |
| `ValueError: strand '+' ya da '-' olmalı` / `Geçersiz bölge` | `region()` parametreleri hatalıysa. |
| `KeyError: "MEFV: 'NM_000243' transkripti yok"` | `Gene.transcript()`'e sürümsüz kimlik verildiyse. |
| `RuntimeError: Kapalı bir oturum üzerinden istek yapılamaz.` | `with` bloğu kapandıktan sonra ağ gerektiren bir çağrı yapıldıysa. |
| `requests.HTTPError` / `ConnectionError` | Ağ hatası; 429/5xx yanıtları önce `max_retries` kez yeniden denenir. |

```python
try:
    g = ncbi.gene(symbol)
except KeyError as e:
    print("Atlandı:", e)
```

## 17. Features anahtarları özeti

| Nesne | Anahtarlar |
|---|---|
| `Gene` (`gene()`) | `source`, `gene_id`, `gene_name`, `description`, `organism`, `chromosome`, `map_location`, `aliases` (her takma ad ayrı değer), `genomic_region` = `(erişim, başlangıç, bitiş, zincir)` |
| `DNA` (`dna()`) | `source`, `chromosome` (erişim), `start`, `end`, `chromosome_strand`, `gene_id`, `gene_name`, `description` |
| `DNA` (`region()`) | `source`, `chromosome`, `start`, `end`, `chromosome_strand` |
| `RNA` (`rna()`) | `source`, `transcript_id`, `title`, `predicted`, `gene_id`, `gene_name`, `protein_id`¹, `utr5`¹, `cds`¹, `utr3`¹ |
| `PRT` (`protein()`) | `source`, `protein_id`, `title`, `transcript_id`² |

¹ Yalnızca protein kodlayan transkriptlerde; `utr5`/`utr3` ancak boş değilse.
² Transkript erişimiyle istendiğinde ya da `gene()` içinden geldiğinde.

`source` her zaman `"NCBI"`dır. Tüm konumlar 1 tabanlı ve iki ucu dahildir.

## 18. Sık sorulanlar ve ipuçları

**`search_gene` neden sadece bir sayı döndürüyor?**
O yalnızca kimlik çözücüdür. İndirmek ve nesne almak için `gene()`, `dna()`,
`rna()`, `protein()` kullanın. Bunlar sembolü de kabul eder.

**Her çalıştırmada "RNA dizileri" çubuğu görünüyor, neden?**
Betikte `update()`'i açıkça çağırıyorsanız RNA FASTA'ları her seferinde
indirilir (bkz. [bölüm 7](#7-update--indirme)). `update()` satırını kaldırın;
`gene()` kayıt hazırsa ağa çıkmaz.

**Verilerim eskidi mi?**
`genes/<id>.json` içindeki `updated` alanı son indirmenin UTC zamanıdır.
Yenilemek için `ncbi.update(gene, force=True)`.

**Gencode ile farkı ne?**
`Gencode` tüm GENCODE sürümünü bir kez indirip yerelde tutar. `Ncbi` ise
**yalnızca istenen genin** RefSeq kayıtlarını indirir. İkisinin `region`,
`dna`, `rna`, `protein`, `gene` metotları aynı mantıkla çalışır.

**Koordinatlar hangi genom sürümüne göre?**
NCBI özetindeki güncel referans (insan için GRCh38). Kromozom erişiminin sürümü
(`NC_000016.10`) bunu belirler.
