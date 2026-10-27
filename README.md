# crn

`crn`, moleküler biyolojinin santral dogmasını (DNA → RNA → protein) Python
nesneleriyle modelleyen bir biyoinformatik çerçevesidir. Üç iş yapar:

1. **Dizileri doğru modellemek.** DNA, RNA ve protein ayrı sınıflardır; her
   biri alfabesini denetler, yazılış yönünü (5'→3' / 3'→5') bilir ve ek
   bilgilerini (`features`) taşır. Dönüşümler biyolojik adlarıyla yapılır:
   `transcription`, `reverse_transcription`, `translation`, `find_orfs`,
   `match_orf`.
2. **Gerçek veriyi tek komutla almak.** NCBI RefSeq ve GENCODE'dan bir genin
   genomik DNA'sı, tüm transkriptleri ve proteinleri tek bir `Gene`
   nesnesinde toplanır. İndirilen her şey yerel `sys/` klasöründe saklanır;
   ikinci istekte ağa çıkılmaz.
3. **Kodon kullanımını ve proteinleri ölçmek.** `crn.cub`, kodon kullanım
   yanlılığının (CUB) literatürdeki temel ölçümlerini (A–H bölümleri),
   `ProteinProperties` ise protein özelliklerini (E-A … E-G) hesaplar. Sayımlar
   gen başına bir kez yapıldığından genom ölçeğinde de hızlıdır.

Ayrıntılı başvuru belgeleri:

| Belge | İçerik |
|---|---|
| [docs/kullanim_kilavuzu.md](docs/kullanim_kilavuzu.md) | Kullanım kılavuzu: tüm modüller (io, fetch, cub, protein), her fonksiyon örnekleriyle, uçtan uca senaryolar |
| [docs/kilavuz.md](docs/kilavuz.md) | Veri modeli ve veri kaynaklarının her metodu, çalıştırılmış örneklerle |
| [docs/features.md](docs/features.md) | `Features` deposu |
| [docs/ncbi.md](docs/ncbi.md) | `Ncbi` ve `NcbiSession` |
| [docs/dogrulanacak_sabitler.md](docs/dogrulanacak_sabitler.md) | Kaynağıyla doğrulanması gereken sabitler |

---

## İçindekiler

1. [Yapı](#1-yapı)
2. [Kurulum](#2-kurulum)
3. [Hızlı başlangıç](#3-hızlı-başlangıç)
4. [Veri modeli](#4-veri-modeli)
5. [Veri kaynakları](#5-veri-kaynakları)
6. [Kodon kullanımı — crn.cub](#6-kodon-kullanımı--crncub)
7. [Protein özellikleri — ProteinProperties](#7-protein-özellikleri--proteinproperties)
8. [İyi uygulamalar ve sık hatalar](#8-iyi-uygulamalar-ve-sık-hatalar)
9. [Sınırlılıklar](#9-sınırlılıklar)
10. [Eski API'den geçiş](#10-eski-apiden-geçiş)

---

## 1. Yapı

```
crn/
├── constants.py   alfabeler, kodonlar, tümleme tabloları, GeneticCode
├── tables.py      kütleler, pKa kümeleri, hidropati ölçekleri, tAI cezaları ...
├── features.py    Features: anahtar -> değer listesi deposu
├── seq.py         Sequence, NucleicAcid, DNA, RNA, PRT; Molecule, Strand, GenomicStrand
├── gene.py        Gene: genomik DNA + transkriptler + proteinler; isoforms()
├── fileio.py      FASTA ve CRN (pickle'sız ikili) dosya biçimleri
├── protein.py     ProteinProperties (E-A … E-G)
├── fetch/         Ncbi, Gencode, NcbiSession; yerel sys/ deposu
└── cub/           CUB (A) ve bölümler B–H; summary_table
```

Bağımlılıklar tek yönlüdür: `fetch` ve `cub` veri modelini kullanır; veri
modeli hiçbirine bağlı değildir. `cub` ise `fetch`'e bağlı değildir. Ona
hangi kaynaktan gelirse gelsin bir `Gene` listesi vermek yeterlidir.

```
   NCBI RefSeq ─┐                         ┌─ FASTA / .crn dosyası
   GENCODE ─────┤  Ncbi().gene("MEFV")    │  read_fasta, Gene.from_crn
                ▼  Gencode().gene("MEFV") ▼
          ┌─────────────────────────────────────────┐
          │ Gene                                    │
          │  ├── dna          genomik DNA (1)       │
          │  ├── transcripts  {id: RNA} (n)         │
          │  ├── proteins     {id: PRT} (m)         │
          │  └── features     gene_id, aliases ...  │
          └─────────────────────────────────────────┘
                │  gene.isoforms()
                ▼
          [Gene, Gene, ...]   her biri: 1 DNA + 1 transkript + 1 protein
                │  CUB(isoforms)        transkript → CDS (proteinle doğrulanmış ORF)
                ▼
          CUB ──► Composition, Synonymous, Adaptation, EffectiveCodons,
                  Divergence, CodonContext, ProteinProperties
                │
                ▼
          cub.table()  → gen başına tüm ölçümler (TSV)   ·   grafikler
```

**Tasarım ilkeleri**

1. **Ek bilgiler `features` içindedir.** Kimlik, kaynak, CDS koordinatı gibi
   bilgiler sınıf alanı değil, anahtar → değer listesi deposudur. Yeni bir
   bilgi taşımak için sınıfı değiştirmek gerekmez.
2. **Root alanları istisnadır.** Yalnızca dizinin biyolojik anlamını
   değiştirenler sınıftadır: `strand` (yazılış yönü), RNA'da `genomic_strand`.
3. **Diziler büyük harfe çevrilir ve doğrulanır.** `DNA("atgx")` →
   `ValueError: DNA dizisinde geçersiz karakter(ler): X`.
4. **Yeni nesne döndüren işlemler bağımsızdır.** Dilim, `complementary()`,
   `transcription()` ... `features`'ın kopyasını taşır. `CUB` da kendisine
   verilen genleri değiştirmez.
5. **Hatalar sessizce yutulmaz.** Belirsiz durumlar `ValueError`/`KeyError`
   ile durur ya da `UserWarning`/`CubWarning` ile bildirilir.
6. **Sabitler tek yerdedir** (`constants.py`, `tables.py`); fonksiyonlara
   parametre olarak girer ve değiştirilebilir.
7. **Ağır işler bir kez yapılır.** Kodon sayımları `CUB` kurulurken; optimal
   kodonlar, w değerleri, CPS skorları ilk kullanımda hesaplanıp saklanır.

---

## 2. Kurulum

Python ≥ 3.12; bağımlılıklar `tqdm`, `requests`, `matplotlib` (yalnızca
grafik çizilirken yüklenir).

```bash
python -m venv .venv
source .venv/bin/activate.fish        # bash/zsh: source .venv/bin/activate
pip install -e .
```

Betikleri proje kökünden, sanal ortamın Python'u ile çalıştırın:

```bash
.venv/bin/python test/cub_test.py
```

Veriler varsayılan olarak proje kökündeki `sys/` klasörüne yazılır;
`CRN_SYS_DIR` ortam değişkeniyle ya da `Ncbi(sys_dir=...)` /
`Gencode(sys_dir=...)` ile değiştirilebilir.

```python
from crn import DNA, RNA, PRT, Gene, Features, Strand, GenomicStrand, GeneticCode, ProteinProperties
from crn.fetch import Ncbi, Gencode, NcbiSession
from crn.cub import (CUB, CubWarning, Composition, Synonymous, Adaptation, EffectiveCodons,
                     Divergence, CodonContext, read_trna_counts)
from crn.constants import CODONS, AMINO_ACIDS, RNA_BASES
from crn.tables import HYDROPATHY_SCALES, PKA_SETS
```

---

## 3. Hızlı başlangıç

```python
from crn.fetch import Ncbi
from crn.cub import CUB, Adaptation, EffectiveCodons

ncbi = Ncbi()
mefv = ncbi.gene("MEFV")              # MEFV: DNA <14607>, 2 transkript, 2 protein
isoforms = mefv.isoforms()            # ['MEFV_NM_000243.3', 'MEFV_NM_001198536.2']
cub = CUB(isoforms)                   # CUB: 2 gen, L=1226

def first_nm(symbol):
    return next(g for g in ncbi.gene(symbol).isoforms()
                if g.features.first("transcript_id").startswith("NM_"))

high = CUB([first_nm(s) for s in ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
                                   "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")])

name = "MEFV_NM_000243.3"
Adaptation(cub, reference=high).cai(name)        # 0.7254
EffectiveCodons(cub).enc(name)                   # 45.44

rows = cub.table(reference=high, file_name="mefv.tsv")   # gen başına 30+ ölçüm
```

---

## 4. Veri modeli

Her metodun ayrıntısı: [docs/kilavuz.md](docs/kilavuz.md).

### 4.1. Diziler

| Sınıf | Alfabe | Root alanları | Başlıca işlemler |
|---|---|---|---|
| `Sequence` | denetlenmez | — | `seq`, `name`, `features`, dilim, `in`, `==`, FASTA/CRN |
| `DNA` | `ACGT` + IUPAC + `-` | `strand` | `complement`, `complementary`, `reverse_complement`, `transcription`, `find_orfs` |
| `RNA` | `ACGU` + IUPAC + `-` | `strand`, `genomic_strand` | `complement`, `complementary`, `reverse_complement`, `reverse_transcription`, `translation`, `match_orf` |
| `PRT` | 20 aa + `BZJXUO*-` | — | — |

`DNA` ve `RNA` ortak `NucleicAcid` tabanından türer. `strand` dizinin
**yazılış yönüdür**: `complement()`/`complementary()` tümleyiciyi bu dizinin
altına hizalı (3'→5') verir, `reverse_complement()` karşı zinciri aynı yönde
(5'→3') verir. `seq_5to3()` dizinin 5'→3' okunuşunu (boşluklar atılmış) döndürür.

```python
>>> d = DNA("ATGAAATTTGGGTAACCCATGCCCTAG", name="g")
>>> [(o.name, o.seq, o.features.first("orf_start")) for o in d.find_orfs()]
[('g_ORF1', 'ATGAAATTTGGGTAA', 1), ('g_ORF2', 'ATGCCCTAG', 19)]
>>> m = d.transcription(template=False)      # kodlayan zincir: yalnızca T->U
>>> m.translation()
g: MKFG*PMP* <9> [PRT] N->C
>>> m.match_orf("MP")                        # çevirisi proteinle aynı ORF
g: AUGCCCUAG <9> [RNA] 5'->3' (?)
```

### 4.2. Gene ve izoformlar

`Gene` bir genomik DNA, istenen sayıda transkript (`{id: RNA}`) ve protein
(`{id: PRT}`) tutar. Kodon kullanımı her izoformu ayrı ister:
`gene.isoforms()` her transkripti, çevirisi bir proteinle birebir aynı ORF'si
üzerinden proteiniyle eşler ve **bir DNA + bir transkript + bir protein**
içeren yeni genler döndürür. Önce kaynağın bağladığı protein
(`transcript_id`) denenir; aynı CDS'i veren kopyalardan biri (önce NM_, sonra
XM_) tutulur. Eşleşmeyenler `UserWarning` ile bildirilir.

### 4.3. Dosyalar

| Biçim | Okuma | Yazma | Kayıp |
|---|---|---|---|
| FASTA | `DNA.from_fasta_file(p)`, `read_fasta(p, cls)` (`.gz` de olur) | `seq.to_fasta(p)`, `write_fasta(seqs, p)` | tip, strand, features |
| CRN (dizi) | `Sequence.from_crn(p)` (doğru sınıfı seçer) | `seq.to_crn(p)` | yok |
| CRN (gen) | `Gene.from_crn(p)` | `gene.to_crn(p)` | yok |

CRN, pickle kullanmayan, CRC32 korumalı ikili biçimdir; tanımı
`crn/fileio.py`'nin başındadır. Eski (sürüm 3) dosyalar okunmaya devam eder.

---

## 5. Veri kaynakları

| | `Ncbi` | `Gencode` |
|---|---|---|
| Kaynak | NCBI E-utilities (RefSeq) | GENCODE FTP (EBI) |
| İndirme | gen başına, istek üzerine (KB'lar) | sürümün tamamı bir kez (~1,2 GB, açılmış ~6 GB) |
| Kimlikler | NM_/NR_/XM_/XR_, NP_/XP_ | ENST, ENSP, ENSG |
| Organizma | insan, fare, sıçan, zebra balığı, sirke sineği, solucan | insan, fare |
| Ağ | ilk istekte (hız sınırı: 3/sn, API anahtarıyla 10/sn) | yalnızca `update()` |

Her iki sınıf da aynı arayüzü verir: `gene(ad)`, `dna(ad)`, `rna(kimlik)`,
`protein(kimlik)`, `region(kromozom, başlangıç, bitiş, zincir)`, `update(...)`.
Dönen nesnelerin `features`'ına `source` ve kaynağa özgü kimlikler yazılır;
transkriptlerde `cds`, `utr5`, `utr3` koordinatları da bulunur.

```python
import os
from crn.fetch import Ncbi, NcbiSession

ncbi = Ncbi(session=NcbiSession(api_key=os.environ.get("NCBI_API_KEY")))
rna = ncbi.rna("NM_000243")                      # sürümsüz de olur
start, end = rna.features.first("cds")
cds = rna[start - 1:end]
```

Analizde tek kaynak kullanın: NCBI ve GENCODE aynı gen için farklı sınırlar,
transkript kümeleri ve kimlikler verir.

---

## 6. Kodon kullanımı — crn.cub

### 6.1. CUB

```python
cub = CUB(isoforms)                        # Gene listesi ya da .crn yolları
cub = CUB(mito_genes, genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL)
panel = CUB([isoforms[0]]) + high + low    # birleştirme; sayımlar yeniden yapılmaz
```

- Her gende **tam olarak bir DNA, bir transkript, bir protein** olmalıdır
  (`gene.isoforms()`); aksi halde hesap başlamadan `ValueError` verilir.
- Transkriptin çevirisi proteinle aynı olan ORF'si (CDS) hesaba girer. Verilen
  genler değiştirilmez; `cub.genes` CDS'e indirgenmiş kopyaları, `cub.cds(ad)`
  CDS dizisini verir.
- Kodon sayımları kurulumda gen başına bir kez yapılır; her ölçüm 64 kodonluk
  sayımdan hesaplanır. Kodon çiftleri ilk istekte sayılıp saklanır.
- Bütün metotlar `gene=None` ise tüm genleri birlikte, gen adı verilirse
  yalnızca o geni hesaplar. Adlar: `cub.names`.

### 6.2. Bölümler ve ölçümler

| Bölüm | Sınıf | Ölçümler (metotlar) |
|---|---|---|
| A. Ham sayımlar | `CUB` | `codon_counts`, `total_codons` (L), `amino_acid_counts`, `codon_freqs`, `codon_per_thousand`; kümeler: `sense_codons`, `codon_families`, `synonymous_codons`, `fourfold_codons` |
| B. Nükleotit bileşimi | `Composition(cub)` | `gc_percent`, `gc1`, `gc2`, `gc12`, `gc3`, `gc3s`, `base3s`, `at3s`, `summary`, `nucleotide_counts/freqs` |
| C. Sinonim tercih | `Synonymous(cub, reference, background)` | `codon_ratio`, `rscu`, `preferred_codons`, `underused_codons`, `optimal_codons` (`"rscu"`/`"chi2"`), `fop` |
| D. Referans indeksleri | `Adaptation(cub, reference, background)` | `relative_adaptiveness` (w), `cai`, `cbi` |
| E. Protein | `ProteinProperties(cub)` | bkz. [§7](#7-protein-özellikleri--proteinproperties) |
| F. ENC ve grafikler | `EffectiveCodons(cub)` | `homozygosity`, `enc`, `expected_enc`, `enc_ratio`, `enc_points`, `neutrality`, `pr2`, `pr2_points`; `plot_enc`, `plot_enc_ratio`, `plot_neutrality`, `plot_pr2` |
| G. Entropi ve sapma | `Divergence(cub, reference)` | `scuo`, `cdc`, `rcdi`, `rcbs`, `expected_codon_freqs` |
| H. tRNA ve bağlam | `CodonContext(cub, reference)` | `trna_weights`, `tai`, `dinucleotide_counts`, `rho`, `codon_pair_scores`, `cpb` |

Formüller ve kaynaklar her modülün başındaki tabloda numaralarıyla yazılıdır
(`crn/cub/*.py`).

### 6.3. Referans setinin seçimi

| Ölçüm | Parametre | Önerilen set |
|---|---|---|
| Fop, CAI, CBI | `reference` | Yüksek ifadeli genler (ribozomal proteinler, EEF1A1) |
| χ² optimal kodonlar | `background` | Düşük ifadeli genler (ör. koku reseptörleri) |
| RCDI | `Divergence(reference=...)` | Konak genomunun geniş ORF seti |
| CPS/CPB | `CodonContext(reference=...)` | Genomun geniş ORF seti |

`reference` verilmezse `cub`'ın kendisi kullanılır; tek gende bu, genin
kendini kendisiyle karşılaştırması demektir. `method="chi2"` için `reference`
ile `background` farklı setler olmalıdır.

### 6.4. Gen paneli: tablo ve grafikler

`cub.table(...)` her gen için bir satır döndürür. Verilen referanslara göre
sütunlar eklenir; `file_name` ile TSV olarak da yazılır (pandas:
`pd.DataFrame(rows)`).

```python
trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
panel = CUB([isoforms[0]]) + high + low
rows = panel.table(reference=high, background=low,       # Fop, CAI, CBI (+ χ²)
                   host=high + low,                       # RCDI
                   genome=high + low,                     # CPB
                   trna_counts=trna,                      # tAI
                   file_name="panel.tsv")
```

```
gen                      GC3s    ENC    CAI    tAI    pI
MEFV_NM_000243.3        0.695  45.44  0.725  0.339  8.30
RPL3_NM_000967.4        0.711  44.68  0.800  0.370 10.19
RPL4_NM_000968.4        0.460  54.76  0.711  0.357 11.07
...
OR1E1_NM_003553.3       0.667  47.02  0.716  0.340  7.98
```

Sütunlar: `gene, L, GC%, GC1, GC2, GC12, GC3, GC3s, ENC, ENC_expected,
ENC_ratio, SCUO, CDC, RCBS`; `reference` ile `Fop, CAI, CBI`; `background`
ile `Fop_chi2, CBI_chi2`; `host` ile `RCDI`; `genome` ile `CPB`;
`trna_counts` ile `tAI`; ve `ProteinProperties.summary` sütunları (`length,
MW, pI, net_charge_7, GRAVY, aromaticity, aliphatic_index, extinction_280,
IVYWREL, CvP, TOP_IDP`).

Dört F grafiği tek şekilde:

```python
import matplotlib.pyplot as plt

eff = EffectiveCodons(panel)
fig, axs = plt.subplots(2, 2, figsize=(11, 10))
eff.plot_enc(ax=axs[0, 0])
eff.plot_enc_ratio(ax=axs[0, 1])
eff.plot_neutrality(ax=axs[1, 0], labels=True)
eff.plot_pr2(ax=axs[1, 1])
fig.savefig("f_grafikleri.png", dpi=300)
```

Bu panelde (25 gen) nötrallik eğimi 0,129 (R² = 0,07): GC3'teki değişimin
ancak ~%13'ü mutasyon baskısıyla açıklanıyor, seçilim baskın. Ham veriler:
`eff.enc_points()`, `eff.neutrality()`, `eff.pr2_points()`.

### 6.5. tRNA havuzu ve kodon bağlamı

```python
ctx = CodonContext(cub, reference=high + low)
w = ctx.trna_weights(trna)                  # {kodon: w}; en büyük UGC, en küçük CUA
ctx.tai(name, trna_counts=trna)             # 0.3386
ctx.tai(name, trna_counts=bakteri_trna, kingdom="prokaryote")   # AUA için lizidin kuralı
ctx.rho(name)["CG"]                         # 0.467 (CpG az temsil)
ctx.rho(name, frame="31")                   # kodon sınırındaki dinükleotitler
ctx.cpb(name)                               # 0.1649
```

### 6.6. Kendi verinizi kullanmak

`CUB` yalnızca `Gene` ister:

```python
cds = "AUGGCUGCCGCAGCGUAA"
g = Gene("sentetik", DNA(cds.replace("U", "T")))
g.add_transcript(RNA(cds), "t1")
g.add_protein(PRT("MAAAA"), "p1")
CUB([g])                                     # CUB: 1 gen, L=5
```

FASTA'dan: `RNA.from_fasta_file("cds.fa")` ve `PRT.from_fasta_file("prot.fa")`
ile okuyup her çift için bir `Gene` kurun.

---

## 7. Protein özellikleri — ProteinProperties

`ProteinProperties` bir `CUB` (proteinler CDS'lerden), bir `PRT` ya da `PRT`
listesiyle çalışır; iki yoldan aynı protein için aynı sonuç çıkar.

```python
props = ProteinProperties(cub)
props.length(name), props.molecular_weight(name)    # 781, 86443.77
props.formula_string(name)                          # 'C3779H6038N1104O1167S26'
props.isoelectric_point(name)                       # 8.30 (Bjellqvist)
props.isoelectric_point(name, pka="EMBOSS")         # 8.11
props.gravy(name)                                   # -0.645
props.summary(name)                                 # başlıca skaler özellikler

ProteinProperties(PRT("MKWVTFISLLLLFSSAYS", name="sinyal")).gravy()   # 1.289
```

| Bölüm | Metotlar |
|---|---|
| E-A. Sayım ve kütle | `length`, `composition`, `group_composition`, `molecular_weight`, `monoisotopic_mass`, `mean_residue_mass`, `residue_counts`, `disulfide_pairs`, `scale_average` |
| E-B. Atom bileşimi | `formula`, `formula_string`, `atom_count`, `atoms_per_residue`, `carbon_oxidation_state` |
| E-C. Yük | `charged_counts`, `simple_net_charge`, `fcr`, `ncpr`, `net_charge` (pH, pKa kümesi), `isoelectric_point`, `titration_curve` |
| E-D. Hidropati | `gravy`, `hydropathy` (Kyte–Doolittle, Hopp–Woods, Eisenberg, Wimley–White, Janin), `aromaticity`, `aliphatic_index`, `boman_index` |
| E-E. Optik | `extinction_coefficient` (`reduced`), `absorbance_01`, `partial_specific_volume`, `mean_residue_volume`, `mean_residue_surface` |
| E-F. Stabilite | `ivywrel`, `cvp_bias`, `solubility_cv`, `solubility_probability`, `composition_entropy` |
| E-G. Düzensizlik | `uversky`, `top_idp`, `composition_features`, `linear_score` |

Birden çok protein birlikte hesaplanırken her biri ayrı zincirdir (kütleye
zincir başına bir H2O, uç grupları ve disülfit eşleri zincir başına).

---

## 8. İyi uygulamalar ve sık hatalar

1. **Her zaman izoformlara ayırın:** `CUB(gene.isoforms())`. `unique_cds=True`'yu
   kapatmayın; aynı CDS'li kopyalar referansları ve regresyonu çarpıtır.
2. **İndirmeyi ve ORF aramasını bir kez yapın**, izoformları `.crn` olarak
   saklayın; `CUB` `.crn` yollarını doğrudan kabul eder.
3. **Ölçüm sınıflarını bir kez kurun**, döngüde yalnızca metot çağırın; ya da
   doğrudan `cub.table()` kullanın.
4. **Referansları açıkça verin**; varsayılan (`cub`'ın kendisi) tek gende
   anlamsızdır.
5. **Uyarıları okuyun.** `CubWarning` referansta olmayan kodonları (RCDI) ve
   çiftleri (CPB) bildirir; susturmak için
   `warnings.filterwarnings("ignore", category=CubWarning)`.
6. **API anahtarını ortam değişkeninde tutun** (`NCBI_API_KEY`).

| Belirti | Neden | Çözüm |
|---|---|---|
| `ValueError: Her gende bir genomik DNA, bir transkript ve bir protein olmalı` | Gen izoformlara ayrılmamış | `gene.isoforms()` |
| `TypeError: genes bir liste olmalı` | `CUB(gene)` | `CUB([gene])` |
| `ValueError: … örtüşen ORF yok` | Protein transkriptten çevrilemiyor (yanlış genetik kod, kodlamayan RNA) | Genetik kodu ve kaydı kontrol edin |
| `ValueError: DNA dizisinde geçersiz karakter(ler)` | Alfabede olmayan harf | Diziyi temizleyin ya da doğru sınıfı kullanın (RNA'da U) |
| Nötrallik `nan` | İkiden az gen | Çok genli panel |
| ENC `nan` | Bir dejenerelik sınıfı için veri yok (çok kısa gen) | Daha uzun gen ya da gen seti |
| `KeyError: CUB: '…' geni yok` | Gen adı yanlış | `cub.names` |
| `ValueError: trna_counts verilmeli` | tAI için tRNA verisi yok | `read_trna_counts(...)` |

---

## 9. Sınırlılıklar

- **ORF eşleme** yalnızca transkriptin kendi zincirinde, AUG ile başlayan ve
  tam eşleşen ORF'lerde yapılır. Selenoproteinler (UGA = Sec), okunan stop
  kodonları ve AUG dışı başlangıçlar eşleşmez.
- **Genetik kodlar:** standart (NCBI 1) ve omurgalı mitokondriyal (NCBI 2).
  Yenisi `GeneticCode`'a 64 harflik NCBI dizisiyle eklenir.
- **tAI** gen kopya sayısını tRNA bolluğunun vekili kabul eder; s değerleri
  maya için optimize edilmiştir (tür özgü stAI yoktur).
- **CDC** Zhang ve ark.'nın bootstrap anlamlılık testini içermez.
- **Wilkinson–Harrison** sabitleri Davis ve ark. (1999) yeniden
  kalibrasyonundandır; özgün makaleyle doğrulanmalıdır
  ([docs/dogrulanacak_sabitler.md](docs/dogrulanacak_sabitler.md)).
- **Protein-Sol gibi skorların** model ağırlıkları pakette yoktur;
  `linear_score`'a dışarıdan verilir.
- Veri modeli ve kaynaklara özgü diğer sınırlılıklar:
  [docs/kilavuz.md §16](docs/kilavuz.md#16-sınırlılıklar).

---

## 10. Eski API'den geçiş

| Eski | Yeni |
|---|---|
| `from crn.io import DNA, ...` | `from crn import DNA, ...` |
| `crn.io.crn_io_const` (alfabeler, `GeneticCode`) | `crn.constants` |
| `crn.io.crn_io_const` (kütleler, ölçekler, pKa, tAI) | `crn.tables` |
| `Maps.*` tabloları, `maps=` parametreleri | `crn.constants.DNA_COMPLEMENT`, `RNA_COMPLEMENT`, `DNA_TO_RNA`, `RNA_TO_DNA`; parametre kaldırıldı |
| `seq._type` | `seq.molecule` |
| `d.complementary().seq[::-1]` | `d.reverse_complement()` |
| `gene.match_transcripts_to_proteins()` | `gene.isoforms()` |
| `CUB(genes=[...])` (gen yerinde kırpılır) | `CUB([...])` (genler değişmez; `cub.genes` kopyalar) |
| `CUB(genes=a.genes + b.genes)` | `a + b` |
| `ProteinMass`, `ProteinAtoms`, `ProteinCharge`, `ProteinHydropathy`, `ProteinOptics`, `ProteinStability`, `ProteinDisorder` | `ProteinProperties` (aynı metot adları) |
| `ProteinDisorder.features()` | `ProteinProperties.composition_features()` |
| Ölçüm metotlarındaki `codons=` parametresi | kaldırıldı (`codon_counts`, `nucleotide_counts`, `pr2` dışında) |
| `RNA.match_orf` sonucunun adı `<ad>_ORF<n>` | RNA'nın kendi adı |
| `crn.fetch.crn_fetch_ncbi` vb. | `crn.fetch.ncbi`, `crn.fetch.gencode`, `crn.fetch.session` |
