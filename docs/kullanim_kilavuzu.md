# crn Kullanım Kılavuzu

Bu kılavuz `crn` çerçevesini baştan sona anlatır: ne işe yaradığını, hangi
modüllerden oluştuğunu, her modülün her fonksiyonunu ve bunların birlikte nasıl
kullanıldığını. Bütün örnekler bu sürümde çalıştırılmıştır. Çıktılar,
yorumlarda (`# ...`) ya da `>>>` satırlarının altında gerçekte üretilen
değerlerdir. Ağ bağlantısı gerektiren ve çalıştırılmamış birkaç örnek ayrıca
belirtilmiştir.

---

## İçindekiler

1. [crn nedir, ne işe yarar?](#1-crn-nedir-ne-işe-yarar)
2. [Modüller: io, fetch, cub](#2-modüller-io-fetch-cub)
3. [Modüllerin fonksiyonlarına genel bakış](#3-modüllerin-fonksiyonlarına-genel-bakış)
4. [io — veri modeli (ayrıntılı başvuru)](#4-io--veri-modeli-ayrıntılı-başvuru)
   - [4.1 crn.constants](#41-crnconstants)
   - [4.2 crn.tables](#42-crntables)
   - [4.3 crn.features — Features](#43-crnfeatures--features)
   - [4.4 crn.seq — Sequence, DNA, RNA, PRT](#44-crnseq--sequence-dna-rna-prt)
   - [4.5 crn.gene — Gene](#45-crngene--gene)
   - [4.6 crn.fileio — FASTA ve CRN](#46-crnfileio--fasta-ve-crn)
5. [fetch — veri kaynakları (ayrıntılı başvuru)](#5-fetch--veri-kaynakları-ayrıntılı-başvuru)
6. [cub — kodon kullanım yanlılığı (ayrıntılı başvuru)](#6-cub--kodon-kullanım-yanlılığı-ayrıntılı-başvuru)
7. [ProteinProperties — protein özellikleri (ayrıntılı başvuru)](#7-proteinproperties--protein-özellikleri-ayrıntılı-başvuru)
8. [Hatalar ve uyarılar](#8-hatalar-ve-uyarılar)
9. [Uçtan uca örnekler](#9-uçtan-uca-örnekler)

---

## 1. crn nedir, ne işe yarar?

`crn`, moleküler biyolojinin **santral dogmasını** (DNA → RNA → protein) Python
nesneleriyle modelleyen ve bu model üzerinde **kodon kullanım yanlılığı
(codon usage bias, CUB)** ile **protein fizikokimyası** analizleri yapan bir
biyoinformatik çerçevesidir.

Çerçevenin çözdüğü üç sorun vardır:

| Sorun | crn'in çözümü |
|---|---|
| Diziler yalnızca metin (`str`) olarak tutulunca yönü, tipi ve kaynağı kaybolur; DNA'ya yanlışlıkla `U` yazılabilir. | `DNA`, `RNA`, `PRT` ayrı sınıflardır. Her biri alfabesini denetler, yazılış yönünü (5'→3' / 3'→5') bilir ve ek bilgilerini `features` içinde taşır. |
| Bir genin genomik DNA'sını, tüm transkriptlerini ve proteinlerini farklı veritabanlarından toplayıp eşleştirmek zahmetlidir. | `fetch` paketi NCBI RefSeq ya da GENCODE'dan tek çağrıyla eksiksiz bir `Gene` kurar. İndirilenler yerel `sys/` klasöründe saklanır, ikinci istekte ağa çıkılmaz. |
| CUB ölçümleri (GC3s, RSCU, CAI, ENC, tAI, CPB...) farklı araçlara dağılmıştır, girdileri tutarsızdır. | `cub` paketi bu ölçümleri tek bir `CUB` nesnesi üzerinden, proteinle doğrulanmış CDS'lerden hesaplar. Gen başına tüm ölçümleri tek tabloda verir. |

Tipik bir akış şöyledir:

```
NCBI / GENCODE / FASTA  ──►  Gene  ──►  gene.isoforms()  ──►  CUB  ──►  ölçümler, tablo, grafik
                              │                                │
                              └── DNA, RNA, PRT nesneleri      └── ProteinProperties
```

**Kimler için?** Belirli bir genin izoformlarının ya da bir gen panelinin
kodon kullanımını karşılaştırmak isteyenler için. Örneğin yüksek ve düşük
ifadeli genleri, konak uyumunu (heterolog ifade) ya da mutasyon ve seçilim
baskısını incelemek isteyen araştırmacılar ve öğrenciler.

**Kurulum** (Python ≥ 3.12):

```bash
python -m venv .venv
source .venv/bin/activate          # fish: source .venv/bin/activate.fish
pip install -e .
```

Bağımlılıklar `requests`, `tqdm` ve `matplotlib`'dir. matplotlib yalnızca
grafik çizilirken yüklenir. Betikleri proje kökünden çalıştırın; örneklerdeki
`sys/...` yolları köke göredir.

---

## 2. Modüller: io, fetch, cub

`crn` üç katmandan oluşur. Bağımlılık tek yönlüdür: `fetch` ve `cub` veri
modelini kullanır, veri modeli hiçbirine bağlı değildir. `cub` de `fetch`'e
bağlı değildir.

```
crn/
├── constants.py   ┐
├── tables.py      │
├── features.py    │  io katmanı: veri modeli
├── seq.py         │  (from crn import ...)
├── gene.py        │
├── fileio.py      │
├── protein.py     ┘  (cub ile ortak kullanılan protein özellikleri)
├── fetch/            veri kaynakları: Ncbi, Gencode, NcbiSession
└── cub/              kodon kullanım yanlılığı: CUB + B–H bölümleri
```

### 2.1 io — veri modeli

Eski sürümlerde `crn.io` adlı ayrı bir paketti. Bu sürümde doğrudan `crn`
paketinin modülleridir ve `from crn import ...` ile kullanılır. Görevi,
biyolojik dizileri ve genleri **doğru, denetlenmiş ve kayıpsız saklanabilir**
nesneler olarak modellemektir.

| Modül | Ne işe yarar |
|---|---|
| `crn.constants` | Alfabeler, 64 kodon, 20 amino asit, tümleme tabloları, `GeneticCode` (standart ve omurgalı mitokondriyal kod). |
| `crn.tables` | Sayısal sabitler: kalıntı kütleleri, atom formülleri, pKa kümeleri, hidropati ölçekleri, ekstinksiyon katsayıları, tAI wobble cezaları. |
| `crn.features` | `Features`: diziye ya da gene ait ek bilgileri (kaynak, kimlik, CDS koordinatı...) anahtar → değer listesi olarak tutar. |
| `crn.seq` | `Sequence`, `DNA`, `RNA`, `PRT` sınıfları ve santral dogma işlemleri: transkripsiyon, ters transkripsiyon, translasyon, ORF arama. |
| `crn.gene` | `Gene`: bir genin genomik DNA'sı, transkriptleri ve proteinleri; `isoforms()` ile izoformlara ayırma. |
| `crn.fileio` | FASTA okuma/yazma ve pickle kullanmayan, CRC32 korumalı ikili **CRN** biçimi. |
| `crn.protein` | `ProteinProperties`: kütle, formül, yük, pI, hidropati, optik, stabilite ve düzensizlik ölçümleri. |

### 2.2 fetch — veri kaynakları

Gerçek genomik veriyi indirip `crn` nesnelerine çevirir. İki kaynak aynı
arayüzü sunar: `gene`, `dna`, `rna`, `protein`, `region`, `update`. Böylece
analiz kodu kaynaktan bağımsız kalır.

| Sınıf | Ne işe yarar |
|---|---|
| `Ncbi` | NCBI E-utilities üzerinden RefSeq kayıtları (NM_/XM_/NR_, NP_/XP_). Gen başına, istek üzerine indirir. İnsan, fare, sıçan, zebra balığı, sirke sineği ve solucan desteklenir. |
| `Gencode` | GENCODE sürümünü (insan, fare) bir kez indirip yerelde indeksler. Sonrası tamamen çevrimdışıdır. |
| `NcbiSession` | NCBI için hız sınırlı (3 istek/sn, API anahtarıyla 10/sn), yeniden denemeli HTTP oturumu. |

### 2.3 cub — kodon kullanım yanlılığı

İzoformlara ayrılmış genlerin CDS'lerinden kodon kullanım ölçümlerini
hesaplar. Merkezde `CUB` nesnesi vardır. Kodonlar kurulumda gen başına bir kez
sayılır, bütün bölümler bu sayımları paylaşır.

| Bölüm | Sınıf | Ne ölçer |
|---|---|---|
| A | `CUB` | Ham sayımlar: kodon/amino asit sayısı, frekanslar, kodon çiftleri; panel tablosu. |
| B | `Composition` | Nükleotit bileşimi: GC%, GC1/2/3, GC3s, A3s/U3s/G3s/C3s. |
| C | `Synonymous` | Eş anlamlı kodon tercihi: RSCU, tercih edilen/az kullanılan/optimal kodonlar, Fop. |
| D | `Adaptation` | Referans sete uyum: w, CAI, CBI. |
| E | `ProteinProperties` | Kodlanan proteinlerin özellikleri (bkz. §7). |
| F | `EffectiveCodons` | ENC, beklenen ENC, nötrallik ve PR2 analizi, dört grafik. |
| G | `Divergence` | Entropi ve sapma: SCUO, CDC, RCDI, RCBS. |
| H | `CodonContext` | tRNA uyumu (tAI), dinükleotit bolluğu (ρ), kodon çifti yanlılığı (CPS/CPB). |

---

## 3. Modüllerin fonksiyonlarına genel bakış

Bu bölüm her fonksiyonu tek satırla tanıtır. Ayrıntı ve örnekler §4–§7'dedir.
"→" dönüş tipini gösterir.

### 3.1 io katmanı

**crn.constants**

| Ad | Ne işe yarar |
|---|---|
| `DNA_BASES`, `RNA_BASES` | `"TCAG"`, `"UCAG"`: kodon sırasını belirleyen baz sıraları. |
| `CODONS` | 64 RNA kodonu, `UUU`'dan `GGG`'ye standart tablo sırasıyla. |
| `AMINO_ACIDS` | 20 standart amino asit, tek harfli: `ACDEFGHIKLMNPQRSTVWY`. |
| `GAP`, `STOP_SYMBOL` | Boşluk (`-`) ve stop (`*`) simgeleri. |
| `DNA_COMPLEMENT`, `RNA_COMPLEMENT` | `str.translate` tabloları: baz → tümleyicisi (IUPAC dahil). |
| `DNA_TO_RNA`, `RNA_TO_DNA` | `T↔U` dönüşüm tabloları. |
| `TRANSCRIPT_PRIORITY` | Aynı CDS'li izoformlardan hangisinin tutulacağı: `('NM_', 'XM_')`. |
| `GeneticCode` | Genetik kodlar: `STANDARD`, `VERTEBRATE_MITOCHONDRIAL`. Özellikleri: `ncbi_id`, `table`, `starts`, `stops`, `sense_codons`, `families`, `synonymous_codons`, `fourfold_codons`. |

**crn.tables**: fonksiyon içermez, yalnızca sabitler vardır (bkz. §4.2).

**crn.features — `Features`**

| Metot | Ne işe yarar |
|---|---|
| `add(key, value)` | Anahtarın değer listesine ekler. |
| `set(key, value)` | Anahtarın değerlerini tek değerle değiştirir. |
| `get(key, default=None)` → list | Değer listesinin kopyası. |
| `first(key, default=None)` | İlk değer. |
| `remove(key)` → bool | Anahtarı siler. |
| `clear()` | Hepsini siler. |
| `keys()`, `items()`, `to_dict()` | Anahtarlar, (anahtar, liste) çiftleri, sözlük kopyası. |
| `copy()` → Features | Derin kopya. |
| `in`, `len`, `for` | Anahtar var mı, anahtar sayısı, anahtarlar üzerinde gezinme. |

**crn.seq**

| Ad | Ne işe yarar |
|---|---|
| `Molecule` | Dizi tipi: `DNA`, `RNA`, `PRT`, `SEQ`. |
| `Strand` | Yazılış yönü: `FivePrime` (5'→3'), `ThreePrime` (3'→5'); `.flipped` tersini verir. |
| `GenomicStrand` | RNA'nın genomdaki zinciri: `Plus`, `Minus`, `Unstranded`, `Unknown`. |
| `Sequence(seq, name)` | Genel dizi (alfabe denetimi yok). `seq`, `name`, `features`, `molecule`, `root_fields()`, dilim, `in`, `==`. |
| `Sequence.from_sequence(s)` | Başka bir diziden bu sınıfta kopya. |
| `.from_fasta_file(p)`, `.to_fasta(p)` | FASTA okuma/yazma. |
| `.from_crn(p)`, `.to_crn(p)` | CRN okuma/yazma (kayıpsız). |
| `DNA` / `RNA`: `strand`, `seq_5to3()` | Yazılış yönü ve 5'→3' okunuş. |
| `complement()` | Yerinde tümleme (3'→5' yazılır). |
| `complementary()` | Tümleyici, yeni nesne. |
| `reverse_complement()` | Karşı zincir 5'→3', yeni nesne. |
| `DNA.transcription(template, genomic_strand)` → RNA | Transkripsiyon. |
| `DNA.find_orfs(min_len, genetic_code, alternative_starts)` → [DNA] | Altı çerçevede ORF arama. |
| `RNA.genomic_strand` | Genomik zincir. |
| `RNA.reverse_transcription()` → DNA | cDNA (ilk zincir). |
| `RNA.translation(genetic_code)` → PRT | Translasyon. |
| `RNA.match_orf(protein, genetic_code)` → RNA \| None | Çevirisi proteinle aynı ORF'yi bulur. |
| `PRT(seq, name)` | Protein dizisi (N→C). |

**crn.gene — `Gene`**

| Metot | Ne işe yarar |
|---|---|
| `Gene(name, dna=None)` | Gen oluşturur; `dna`, `features`. |
| `add_transcript(rna, id=None)`, `add_protein(prt, id=None)` | Ekler (kimlik yoksa nesnenin adı). |
| `transcript(id)`, `protein(id)` | Kimlikle getirir (`KeyError`). |
| `remove_transcript(id)`, `remove_protein(id)` → bool | Siler. |
| `transcripts`, `proteins` | `{kimlik: nesne}` kopyaları. |
| `isoforms(genetic_code, unique_cds, priority)` → [Gene] | Her biri 1 DNA + 1 transkript + 1 protein olan izoform genleri. |
| `to_crn(p)`, `Gene.from_crn(p)` | Geni kayıpsız saklar ve okur. |

**crn.fileio**

| Fonksiyon | Ne işe yarar |
|---|---|
| `read_fasta(p, cls=Sequence)` → list | FASTA'yı (`.gz` dahil) okur. |
| `write_fasta(seqs, p)` | FASTA yazar. |
| `write_crn(seq, p)`, `read_crn(p, cls=Sequence)` | Tek diziyi CRN olarak yazar ve okur. |
| `write_gene_crn(gene, p)`, `read_gene_crn(p)` | Geni CRN olarak yazar ve okur. |

**crn.protein — `ProteinProperties`**: 40'tan fazla metodun özeti §7'nin başındaki tablodadır.

### 3.2 fetch

| Sınıf / metot | Ne işe yarar |
|---|---|
| `NcbiSession(api_key, e_mail, timeout, max_retries, ...)` | Hız sınırlı HTTP oturumu: `get(url)`, `is_open()`, `close()`. |
| `Ncbi(organism, session, sys_dir, verbose)` | NCBI kaynağı. |
| `Ncbi.search_gene(sembol)` → str | Gen sembolü → NCBI Gene ID. |
| `Ncbi.update(gen, force)` → str | Genin tüm kayıtlarını indirir. |
| `Ncbi.gene(gen)` → Gene | DNA + tüm RefSeq transkriptleri ve proteinleri. |
| `Ncbi.dna(gen)`, `rna(acc)`, `protein(acc)`, `region(acc, başlangıç, bitiş, zincir)` | Tek tek nesneler. |
| `Gencode(species, release, sys_dir, base_url, verbose)` | GENCODE kaynağı. |
| `Gencode.update(force)` → str | Sürümü indirip hazırlar (~1,2 GB). |
| `Gencode.latest_release()`, `installed_releases()`, `is_ready()`, `release` | Sürüm bilgisi. |
| `Gencode.gene(gen)`, `dna`, `rna`, `protein`, `region(kromozom, ...)` | Nesneler. |

### 3.3 cub

| Sınıf / fonksiyon | Metotlar |
|---|---|
| `CUB(genes, genetic_code)` | `genes`, `names`, `genetic_code`, `cds`, `codons`, `sense_codons`, `codon_families`, `synonymous_codons`, `fourfold_codons`, `codon_counts`, `total_codons`, `amino_acid_counts`, `codon_freqs`, `codon_per_thousand`, `termini`, `chain_count`, `codon_pair_counts`, `table`, `+` |
| `Composition(cub)` | `nucleotide_counts`, `nucleotide_freqs`, `gc_percent`, `gc1`, `gc2`, `gc12`, `gc3`, `gc3s`, `base3s`, `at3s`, `summary` |
| `Synonymous(cub, reference, background)` | `codon_ratio`, `rscu`, `preferred_codons`, `underused_codons`, `optimal_codons`, `fop` |
| `Adaptation(cub, reference, background)` | `relative_adaptiveness`, `cai`, `cbi` |
| `EffectiveCodons(cub)` | `homozygosity`, `enc`, `expected_enc`, `enc_ratio`, `enc_points`, `neutrality`, `pr2`, `pr2_points`, `plot_enc`, `plot_enc_ratio`, `plot_neutrality`, `plot_pr2` |
| `Divergence(cub, reference)` | `expected_codon_freqs`, `scuo`, `cdc`, `rcdi`, `rcbs` |
| `CodonContext(cub, reference)` | `trna_weights`, `tai`, `dinucleotide_counts`, `rho`, `codon_pair_scores`, `cpb` |
| `read_trna_counts(p)` | GtRNAdb FASTA → `{antikodon: kopya sayısı}` |
| `summary_table(cub, ...)` | `cub.table()`'ın kendisi; tüm ölçümler gen başına |
| `CubWarning` | Hesaptan çıkarılan veriyi bildiren uyarı sınıfı |

---

## 4. io — veri modeli (ayrıntılı başvuru)

```python
from crn import (Sequence, DNA, RNA, PRT, Gene, Features, Strand, GenomicStrand,
                 GeneticCode, ProteinProperties, read_fasta, write_fasta)
from crn.seq import Molecule
from crn import constants, tables
```

### 4.1 crn.constants

Çerçevenin bütün biyolojik sabitleri buradadır. Fonksiyonlar bunları
varsayılan parametre olarak kullanır, dolayısıyla isterseniz kendi değerinizi
verebilirsiniz.

**Alfabeler ve kodonlar**

```python
>>> from crn.constants import CODONS, AMINO_ACIDS, DNA_BASES, RNA_BASES, GAP, STOP_SYMBOL
>>> len(CODONS), CODONS[:4]
(64, ('UUU', 'UUC', 'UUA', 'UUG'))
>>> AMINO_ACIDS, DNA_BASES, RNA_BASES, GAP, STOP_SYMBOL
('ACDEFGHIKLMNPQRSTVWY', 'TCAG', 'UCAG', '-', '*')
```

**Ne işe yarar?** `CODONS` bütün kodon sözlüklerinin anahtar sırasını belirler.
Bu yüzden bir ölçümün çıktısı her zaman aynı sırada gelir ve doğrudan tabloya
dökülebilir.

**Tümleme ve dönüşüm tabloları.** Bunlar `str.translate` için hazırlanmış
tablolardır. Hem büyük hem küçük harfleri ve IUPAC belirsiz bazlarını
(`R↔Y`, `K↔M`, `N↔N`...) kapsar.

```python
>>> from crn.constants import DNA_COMPLEMENT, RNA_COMPLEMENT, DNA_TO_RNA, RNA_TO_DNA
>>> "ATGC".translate(DNA_COMPLEMENT), "AUGC".translate(RNA_COMPLEMENT)
('TACG', 'UACG')
>>> "ATGC".translate(DNA_TO_RNA), "AUGC".translate(RNA_TO_DNA)
('AUGC', 'ATGC')
```

**`TRANSCRIPT_PRIORITY = ('NM_', 'XM_')`**: `Gene.isoforms()` aynı CDS'i veren
transkriptlerden birini tutarken önce elle doğrulanmış RefSeq (NM_), sonra
tahmini (XM_) kaydı seçer.

#### GeneticCode

Translasyon, ORF araması ve kodon sayımı hangi genetik kodla yapılacaksa o
`GeneticCode` üyesi verilir. Her üye NCBI tablo numarasını ve 64 harflik
amino asit dizisini taşır. Aşağıdaki özellikler bunlardan türetilip
önbelleğe alınır.

| Özellik | Ne verir |
|---|---|
| `ncbi_id` | NCBI tablo no (1, 2) |
| `table` | `{kodon: amino asit}`; stop `*` |
| `starts` | Başlangıç kodonları (frozenset) |
| `stops` | Stop kodonları |
| `sense_codons` | Amino asit kodlayan kodonlar |
| `families` | `{amino asit: [kodonlar]}`: eş anlamlı kodon aileleri |
| `synonymous_codons` | Birden çok kodonu olan amino asitlerin kodonları |
| `fourfold_codons` | 3. konumu tümüyle eş anlamlı (xyN) kutuların kodonları |

```python
>>> gc = GeneticCode.STANDARD
>>> gc.ncbi_id, gc.table["AUG"], gc.table["UGA"]
(1, 'M', '*')
>>> gc.starts, gc.stops
(frozenset({'UUG', 'CUG', 'GUG', 'AUG'}), frozenset({'UGA', 'UAG', 'UAA'}))
>>> len(gc.sense_codons), gc.families["L"]
(61, ('UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'))
>>> len(gc.synonymous_codons), len(gc.fourfold_codons)
(59, 32)

>>> mt = GeneticCode.VERTEBRATE_MITOCHONDRIAL
>>> mt.ncbi_id, mt.table["UGA"], mt.table["AUA"]     # UGA = Trp, AUA = Met
(2, 'W', 'M')
>>> mt.stops
frozenset({'AGG', 'UAG', 'AGA', 'UAA'})
```

**Nasıl kullanılır?** `genetic_code=` parametresi alan her yere verilir:
`RNA.translation`, `DNA.find_orfs`, `RNA.match_orf`, `Gene.isoforms`, `CUB`.
`starts`/`stops` küme olduğundan yazdırıldığında sırası değişebilir.

### 4.2 crn.tables

Protein ve tAI hesaplarının sayısal sabitleridir. Her sabit ilgili metodun
varsayılan parametresidir. Başka bir kaynağın değerleriyle hesap yapmak
için sözlüğü kopyalayıp değiştirin ve metoda verin.

| Sabit | İçerik | Kullanan metot |
|---|---|---|
| `RESIDUE_AVERAGE_MASS`, `RESIDUE_MONOISOTOPIC_MASS`, `WATER_AVERAGE_MASS`, `WATER_MONOISOTOPIC_MASS` | Kalıntı ve su kütleleri (Da) | `molecular_weight`, `monoisotopic_mass` |
| `AMINO_ACID_FORMULAS`, `WATER_FORMULA`, `ELEMENTS` | (C, H, N, O, S) atom sayıları | `formula` |
| `PKA_SETS` | `'EMBOSS'`, `'Lehninger'`, `'Bjellqvist'` pKa kümeleri | `net_charge`, `isoelectric_point` |
| `POSITIVE_RESIDUES`, `NEGATIVE_RESIDUES`, `IONIZABLE_POSITIVE`, `IONIZABLE_NEGATIVE` | Yüklü kalıntı kümeleri | yük metotları |
| `HYDROPATHY_SCALES` | `'Kyte-Doolittle'`, `'Hopp-Woods'`, `'Eisenberg'`, `'Wimley-White'`, `'Janin'` | `hydropathy`, `gravy` |
| `BOMAN_SCALE` | Kalıntı çözünme serbest enerjileri | `boman_index` |
| `EXTINCTION_COEFFICIENTS_280` | W, Y, sistin ε (280 nm) | `extinction_coefficient` |
| `RESIDUE_SPECIFIC_VOLUMES`, `RESIDUE_VOLUMES`, `RESIDUE_MAX_ASA` | Kısmi hacim, kalıntı hacmi, en büyük yüzey | optik metotlar |
| `AMINO_ACID_GROUPS` | hydrophobic, polar, aliphatic, aromatic, small, basic, acidic | `group_composition` |
| `IVYWREL_RESIDUES`, `CVP_CHARGED`, `CVP_POLAR`, `TURN_FORMING_RESIDUES` | Stabilite grupları | `ivywrel`, `cvp_bias`, `solubility_cv` |
| `TOP_IDP_SCALE` | Düzensizlik eğilim ölçeği | `top_idp` |
| `TAI_S_VALUES`, `TAI_WOBBLE_PAIRS` | tAI wobble cezaları ve eşleşme kuralları | `trna_weights`, `tai` |

```python
>>> from crn.tables import HYDROPATHY_SCALES, PKA_SETS
>>> list(HYDROPATHY_SCALES)
['Kyte-Doolittle', 'Hopp-Woods', 'Eisenberg', 'Wimley-White', 'Janin']
>>> list(PKA_SETS)
['EMBOSS', 'Lehninger', 'Bjellqvist']
```

Kaynağıyla doğrulanması gereken sabitler
[dogrulanacak_sabitler.md](dogrulanacak_sabitler.md) dosyasında listelenmiştir.

### 4.3 crn.features — Features

`Features`, bir dizi ya da gen hakkındaki **ek bilgileri** tutar. Her anahtarın
değeri bir **listedir**, çünkü bir genin birden çok takma adı ya da bir
transkriptin birden çok exon koordinatı olabilir. Yeni bir bilgi taşımak için
sınıfa yeni alan eklemek gerekmez, yalnızca `features`'a yazmak yeter.

`fetch`'in yazdığı başlıca anahtarlar: `source`, `gene_id`, `gene_name`,
`transcript_id`, `protein_id`, `cds`, `utr5`, `utr3`, `chromosome`, `start`,
`end`, `chromosome_strand`, `aliases`.

| Metot | Ne yapar | Örnek |
|---|---|---|
| `add(key, value)` | Listeye ekler | `f.add("aliases", "FMF")` |
| `set(key, value)` | Listeyi `[value]` yapar | `f.set("source", "NCBI")` |
| `get(key, default=None)` | Listenin **kopyası** | `f.get("aliases")` |
| `first(key, default=None)` | İlk değer | `f.first("cds")` |
| `remove(key)` | Siler; varsa `True` | `f.remove("source")` |
| `clear()` | Hepsini siler | |
| `keys()`, `items()` | Anahtarlar, (anahtar, liste) | |
| `to_dict()` | `{anahtar: liste}` kopyası | |
| `copy()` | Bağımsız derin kopya | |

```python
>>> f = Features()
>>> f.add("aliases", "FMF"); f.add("aliases", "MEF"); f.set("source", "NCBI")
>>> f
Keys:aliases,source
>>> f.get("aliases"), f.first("aliases"), f.first("yok", "-")
(['FMF', 'MEF'], 'FMF', '-')
>>> "source" in f, len(f), f.keys()
(True, 2, ['aliases', 'source'])
>>> f.to_dict()
{'aliases': ['FMF', 'MEF'], 'source': ['NCBI']}

>>> g = f.copy(); g.add("aliases", "X")      # kopya bağımsızdır
>>> f.get("aliases"), g.get("aliases")
(['FMF', 'MEF'], ['FMF', 'MEF', 'X'])
>>> f.remove("source"), f.remove("source"), list(f)
(True, False, ['aliases'])
```

**İpucu:** `get` bir kopya döndürür. `f.get("x").append(...)` deposu
değiştirmez, ekleme için `add` kullanın.

### 4.4 crn.seq — Sequence, DNA, RNA, PRT

#### Ortak davranış (Sequence)

| Üye | Ne işe yarar |
|---|---|
| `seq` | Dizi metni (her zaman büyük harf). Atanabilir, atanırken doğrulanır. |
| `name` | Ad. |
| `features` | `Features` deposu. |
| `molecule` | `Molecule.DNA/RNA/PRT/SEQ`: sınıfın tipi. |
| `root_fields()` | Sınıfa ait alanlar ve değerleri, ör. `{"strand": Strand.FivePrime}`. |
| `len(s)`, `s[i:j]`, `"ATG" in s`, `s1 == s2` | Uzunluk; dilim (aynı sınıfta, features kopyalı yeni nesne); alt dizi arama; tip, ad, dizi ve root alanlarına göre eşitlik. |
| `from_sequence(s)` | Başka bir diziden bu sınıfta kopya oluşturur. |
| `from_fasta_file(p)`, `to_fasta(p)` | FASTA (bkz. §4.6). |
| `from_crn(p)`, `to_crn(p)` | CRN (bkz. §4.6). |

`Sequence` alfabe denetimi yapmaz. Dizinin tipi bilinmediğinde (ör. FASTA'dan
okurken) kullanılır.

```python
>>> s = Sequence("hello world?", name="x")
>>> s, s.molecule, s.root_fields()
(x: HELLO WORLD? <12> [SEQ], <Molecule.SEQ: 'SEQ'>, {})
```

`repr` biçimi şudur: `ad: DİZİ <uzunluk> [tip] yön (genomik zincir)`. Uzun
dizilerin ortası `...` ile kısaltılır.

#### DNA

`DNA(seq, name="", strand=Strand.FivePrime)`. Alfabe `ACGT`, IUPAC belirsiz
bazlar (`RYKMSWBDHVN`) ve `-`'dir. Küçük harf büyüğe çevrilir, geçersiz
karakter `ValueError` verir.

```python
>>> d = DNA("atgcgtacgtta", name="d1")
>>> d, d.molecule, d.strand
(d1: ATGCGTACGTTA <12> [DNA] 5'->3', <Molecule.DNA: 'DNA'>, <Strand.FivePrime: ...>)
>>> len(d), d[0:3], "CGT" in d
(12, d1: ATG <3> [DNA] 5'->3', True)
>>> d == DNA("ATGCGTACGTTA", name="d1")
True
>>> DNA("ATGU")
ValueError: DNA dizisinde geçersiz karakter(ler): U
```

**Strand: yazılış yönü.** `strand` dizinin hangi uçtan okunduğunu söyler.
Biyolojik zinciri değil, metnin yönünü belirtir. `Strand.FivePrime.flipped`
değeri `Strand.ThreePrime` olur.

##### seq_5to3()

Dizinin 5'→3' okunuşunu verir. Boşluklar (`-`) atılır, `ThreePrime` dizi ters
çevrilir. Yönden bağımsız işlem yapmanız gerektiğinde (motif arama gibi)
kullanın.

```python
>>> DNA("ATG-CA", name="g").seq_5to3()
'ATGCA'
```

##### complement(), complementary(), reverse_complement()

Üçü de tümleyici zinciri verir. Aralarındaki fark yön ve yerinde değişimdir:

| Metot | Sonuç | Yön | Yeni nesne mi? |
|---|---|---|---|
| `complement()` | Tümleyici, bu dizinin **altına hizalı** | 3'→5' | Hayır, yerinde değiştirir ve `self` döner |
| `complementary()` | Aynısı | 3'→5' | Evet |
| `reverse_complement()` | Karşı zincir, **aynı yönde** okunur | 5'→3' | Evet |

```python
>>> c = d.complementary()
>>> c, c.seq_5to3()
(d1: TACGCATGCAAT <12> [DNA] 3'->5', 'TAACGTACGCAT')
>>> d.reverse_complement()                      # primer tasarımı için olağan biçim
d1: TAACGTACGCAT <12> [DNA] 5'->3'

>>> z = DNA("ATGC", name="z"); z.complement() is z
True
>>> z
z: TACG <4> [DNA] 3'->5'
```

**Ne zaman hangisi?** Çift zinciri alt alta göstermek için
`complementary()`, karşı zincirde 5'→3' çalışmak için (primer, ters zincir
ORF'si) `reverse_complement()` kullanın.

##### transcription(template=True, genomic_strand=GenomicStrand.Unknown) → RNA

RNA polimerazın yaptığını modeller. `template=True` ise bu DNA **kalıp
zincir** kabul edilir: RNA tümleyicidir, antiparaleldir (3'→5' yazılır).
`template=False` ise bu DNA **kodlayan zincirdir**: RNA aynı dizidir, yalnızca
T→U değişir. `genomic_strand` RNA'nın genomdaki zincirini kaydeder.

```python
>>> d.transcription()
d1: UACGCAUGCAAU <12> [RNA] 3'->5' (?)
>>> d.transcription(template=False, genomic_strand=GenomicStrand.Plus)
d1: AUGCGUACGUUA <12> [RNA] 5'->3' (+)
```

**İpucu:** NCBI ve GENCODE'dan gelen genomik DNA kodlayan zincirde, 5'→3'
yazılır. mRNA dizisini elde etmek için `template=False` kullanın.

##### find_orfs(minimum_length_of_orf=-1, genetic_code=STANDARD, alternative_starts=False) → [DNA]

Her iki zincirde ve üç çerçevede (toplam altı çerçeve) başlangıç kodonundan
aynı çerçevedeki ilk stop'a kadar olan **açık okuma çerçevelerini** bulur.
Stop'suz kısmi ORF'ler alınmaz. Aynı stop'u paylaşan iç içe ORF'lerin hepsi
döner.

| Parametre | Anlamı |
|---|---|
| `minimum_length_of_orf` | Stop dahil en küçük uzunluk (nt); `-1` sınırsız |
| `genetic_code` | Başlangıç/stop kodonlarının alınacağı kod |
| `alternative_starts` | `True` ise kodun tüm başlangıçları (CUG, GUG, UUG), `False` ise yalnızca ATG |

Dönen her ORF kodlayan yönde 5'→3' yazılmış bir `DNA`'dır. Adı
`<ad>_ORF<n>`'dir. `features`'a `orf_start`, `orf_end` (bu dizide 1 tabanlı,
her zaman start ≤ end), `orf_strand` (`+`/`-`) ve `orf_frame` (0, 1, 2)
yazılır.

```python
>>> g = DNA("ATGAAATTTGGGTAACCCATGCCCTAG", name="g")
>>> for orf in g.find_orfs():
...     print(orf, orf.features.to_dict())
g_ORF1: ATGAAATTTGGGTAA <15> [DNA] 5'->3' {'orf_start': [1], 'orf_end': [15], 'orf_strand': ['+'], 'orf_frame': [0]}
g_ORF2: ATGCCCTAG <9> [DNA] 5'->3' {'orf_start': [19], 'orf_end': [27], 'orf_strand': ['+'], 'orf_frame': [0]}

>>> [o.seq for o in g.find_orfs(minimum_length_of_orf=12)]
['ATGAAATTTGGGTAA']

>>> [o.seq for o in DNA("CTACAT", name="rc").find_orfs()]   # ters zincirde ORF
['ATGTAG']
```

#### RNA

`RNA(seq, name="", strand=Strand.FivePrime, genomic_strand=GenomicStrand.Unknown)`.
Alfabe `ACGU`, IUPAC ve `-`'dir. DNA'daki `strand`, `seq_5to3`, `complement`,
`complementary` ve `reverse_complement` aynen geçerlidir. `genomic_strand`
tümlemeyle değişmez.

##### reverse_transcription() → DNA

Ters transkriptazı modeller ve ilk zincir cDNA'yı üretir. Sonuç RNA'nın
tümleyicisidir ve antiparaleldir.

```python
>>> m = d.transcription(template=False, genomic_strand=GenomicStrand.Plus)
>>> m.reverse_transcription()
d1: TACGCATGCAAT <12> [DNA] 3'->5'
```

##### translation(genetic_code=STANDARD) → PRT

İlk nükleotitten son tam kodona kadar üçerli okur. Başlangıç ya da stop
aranmaz: stoplar `*` yazılır ve okuma sürer. Belirli bir ORF'yi çevirmek için
önce ORF'yi kesin (`find_orfs`, `match_orf` ya da `features["cds"]` ile
dilim).

```python
>>> RNA("AUGGCCUAA", name="r").translation()
r: MA* <3> [PRT] N->C
>>> RNA("AUGUGAUAA", name="mt").translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL)
mt: MW* <3> [PRT] N->C
>>> DNA("ATGAAATTTGGGTAACCCATGCCCTAG", name="g").transcription(template=False).translation()
g: MKFG*PMP* <9> [PRT] N->C
```

##### match_orf(protein, genetic_code=STANDARD) → RNA | None

Transkript içinde, çevirisi verilen proteinle **birebir aynı** olan AUG...stop
ORF'sini bulur. Transkript–protein eşleşmesini **dizi düzeyinde doğrular**. Bu
yüzden `Gene.isoforms()` ve `CUB` bunun üzerine kuruludur. `protein` bir `PRT`
ya da metin olabilir. Sondaki `*` olsa da olmasa da eşleşir. Birden çok eşleşme
varsa en 5' uçtaki seçilir.

Dönen RNA transkriptin adını taşır ve yalnızca ORF dizisidir (stop dahil).
`features`'ta `cds = (1, uzunluk)`, transkriptteki konumu gösteren `orf_start`
ve `orf_end` bulunur. Artık geçersiz olan `utr5` ve `utr3` silinir.

```python
>>> m = DNA("ATGAAATTTGGGTAACCCATGCCCTAG", name="g").transcription(template=False)
>>> orf = m.match_orf("MP")
>>> orf, orf.features.to_dict()
(g: AUGCCCUAG <9> [RNA] 5'->3' (?), {'cds': [(1, 9)], 'orf_start': [19], 'orf_end': [27]})
>>> m.match_orf("MKF")          # "MKFG" ORF'si var ama "MKF" birebir değil
None
```

#### PRT

`PRT(seq, name="")`. Alfabe 20 amino asit, `BZJXUO` (belirsiz ve nadir) ve `*`
(stop) ile `-`'dir. N-ucundan C-ucuna yazılır. Hesaplamalar
`ProteinProperties` ile yapılır (§7).

```python
>>> PRT("mkwvtf*", name="p")
p: MKWVTF* <7> [PRT] N->C
```

#### from_sequence(s)

Bir diziyi başka bir sınıfa çevirir. Ad, dizi ve `features` kopyalanır. Tipi
belirsiz okunan diziyi DNA'ya çevirmek için kullanılır (alfabe bu sırada
denetlenir).

```python
>>> DNA.from_sequence(Sequence("ATGC", name="s"))
s: ATGC <4> [DNA] 5'->3'
```

### 4.5 crn.gene — Gene

`Gene`, bir genin dizilerini bir arada tutan kapsayıcıdır: **bir** genomik DNA,
istenen sayıda transkript (`{kimlik: RNA}`) ve protein (`{kimlik: PRT}`).
Kimlikler ekleme sırasıyla tutulur.

```python
>>> g = Gene("demo", DNA("ATGGCTGCCTAAGG", name="demo"))
>>> g.add_transcript(RNA("AUGGCUGCCUAAGG", name="t1"))          # kimlik = ad
>>> g.add_transcript(RNA("GGAUGGCUGCCUAA", name="t2"), "T2")    # açık kimlik
>>> g.add_protein(PRT("MAA", name="p1"))
>>> g.features.add("gene_id", "42")
>>> g
demo: DNA <14>, 2 transkript, 1 protein
>>> list(g.transcripts), list(g.proteins)
(['t1', 'T2'], ['p1'])
>>> g.transcript("T2"), g.protein("p1")
(t2: GGAUGGCUGCCUAA <14> [RNA] 5'->3' (?), p1: MAA <3> [PRT] N->C)
```

| Metot | Ne yapar | Hata |
|---|---|---|
| `add_transcript(rna, transcript_id=None)` | Ekler, RNA'yı döndürür | Aynı kimlik varsa `ValueError` |
| `add_protein(prt, protein_id=None)` | Ekler, PRT'yi döndürür | Aynı kimlik varsa `ValueError` |
| `transcript(id)`, `protein(id)` | Getirir | Yoksa `KeyError` |
| `remove_transcript(id)`, `remove_protein(id)` | Siler; varsa `True` | |
| `transcripts`, `proteins` | Sözlüğün kopyası | |

```python
>>> g.add_protein(PRT("MAA", name="p1"))
ValueError: demo: 'p1' proteini zaten var
>>> g.transcript("yok")
KeyError: "demo: 'yok' transkripti yok"
```

#### isoforms(genetic_code=STANDARD, unique_cds=True, priority=('NM_', 'XM_')) → [Gene]

Kodon kullanımı her izoformu ayrı ister, `CUB` de bu yüzden her gende tam
**bir DNA + bir transkript + bir protein** bekler. `isoforms()` bu ayrımı
yapar:

1. Her transkript için çevirisi bir proteinle birebir aynı ORF aranır
   (`match_orf`). Önce kaynağın transkripte bağladığı protein
   (`features["protein_id"]`) denenir.
2. Eşleşen her çift için `<gen>_<transkript>` adlı yeni bir `Gene` kurulur.
   Genin features'ı kopyalanır, `isoform_of`, `transcript_id` ve `protein_id`
   eklenir.
3. `unique_cds=True` ise aynı CDS'i veren izoformlardan yalnızca biri tutulur
   (`priority` sırasıyla, önce NM_). Aynı dizi iki kez sayılırsa ölçümler ve
   regresyon çarpılır, bu yüzden varsayılanı değiştirmeyin.
4. Proteini olmayan (kodlamayan) ya da eşleşmeyen transkriptler `UserWarning`
   ile bildirilir.

```python
>>> iso = g.isoforms()
UserWarning: demo: aynı CDS'i veren izoformlar çıkarıldı (T2 -> t1)
>>> iso
[demo_t1: DNA <14>, 1 transkript, 1 protein]
>>> iso[0].features.to_dict()
{'gene_id': ['42'], 'isoform_of': ['demo'], 'transcript_id': ['t1'], 'protein_id': ['p1']}
>>> [i.name for i in g.isoforms(unique_cds=False)]
['demo_t1', 'demo_T2']
```

Gerçek bir genle (NCBI'dan MEFV):

```python
>>> Ncbi().gene("MEFV").isoforms()
[MEFV_NM_000243.3: ..., MEFV_NM_001198536.2: ...]
```

#### to_crn(p), Gene.from_crn(p)

Geni bütün içeriğiyle (ad, features, DNA, transkriptler, proteinler) tek bir
CRN dosyasına yazar ve geri okur. İndirme ve izoform ayırma işini bir kez
yapıp sonucu saklamak için idealdir. `CUB` `.crn` yollarını doğrudan kabul
eder.

```python
>>> g.to_crn("demo.crn")
>>> Gene.from_crn("demo.crn")
demo: DNA <14>, 2 transkript, 1 protein
```

### 4.6 crn.fileio — FASTA ve CRN

| Biçim | Ne saklar | Ne kaybeder | Ne zaman |
|---|---|---|---|
| FASTA | ad + dizi | tip, strand, features | Başka araçlarla veri alışverişi |
| CRN | her şey | hiçbir şey | crn içinde saklama, önbellek |

#### read_fasta(file_name, cls=Sequence) → list

FASTA dosyasındaki bütün kayıtları `cls` tipinde okur. `.gz` dosyalar
doğrudan açılır. Dizinin adı başlığın ilk kelimesidir. Kısa yolu
`DNA.from_fasta_file(p)`'dir (aynı şey, `cls=DNA`).

#### write_fasta(sequences, file_name)

Dizileri FASTA'ya yazar, dosya varsa üzerine yazılır. Kısa yolu
`seq.to_fasta(p)`'dir.

```python
>>> write_fasta([d, DNA("ATGAAATTTGGGTAACCCATGCCCTAG", name="g")], "ex.fa")
```
```
>d1 this sequence contains 0 features but not listed
ATGCGTACGTTA
>g this sequence contains 0 features but not listed
ATGAAATTTGGGTAACCCATGCCCTAG
```
```python
>>> DNA.from_fasta_file("ex.fa")
[d1: ATGCGTACGTTA <12> [DNA] 5'->3', g: ATGAAATTTGGGTAACCCATGCCCTAG <27> [DNA] 5'->3']
>>> [type(x).__name__ for x in read_fasta("ex.fa")]     # cls verilmezse Sequence
['Sequence', 'Sequence']
```

#### write_crn(seq, p), read_crn(p, cls=Sequence)

Tek bir diziyi kayıpsız yazar ve okur. Kısa yolları `seq.to_crn(p)` ve
`Sequence.from_crn(p)`'dir. `Sequence.from_crn` dosyadaki tipe uygun sınıfı
döndürür. `DNA.from_crn` gibi belirli bir sınıfla okunursa tip tutmalıdır.

```python
>>> d.features.add("source", "demo")
>>> d.to_crn("d.crn")
>>> back = Sequence.from_crn("d.crn")
>>> back, type(back).__name__, back.features.to_dict()
(d1: ATGCGTACGTTA <12> [DNA] 5'->3', 'DNA', {'source': ['demo']})
>>> RNA.from_crn("d.crn")
ValueError: d.crn: dosyadaki dizi tipi DNA, RNA sınıfı RNA tipi bekliyor
```

#### write_gene_crn(gene, p), read_gene_crn(p)

`Gene.to_crn` ve `Gene.from_crn`'in kendisidir (bkz. §4.5).

**CRN biçimi.** Sihirli bayt + sürüm + kayıt + CRC32'den oluşur. pickle
kullanılmadığından güvenilmeyen dosyadan kod çalıştırılamaz, bozuk dosya
CRC ile yakalanır. Ayrıntılı tanım `crn/fileio.py`'nin başındadır. Eski
sürümün dosyaları okunmaya devam eder.

---

## 5. fetch — veri kaynakları (ayrıntılı başvuru)

```python
from crn.fetch import Ncbi, Gencode, NcbiSession
```

### 5.1 Yerel depo: sys/

Bütün kaynaklar indirdiklerini yerel bir klasöre yazar. Sıra şöyledir:
`sys_dir=` parametresi, yoksa `CRN_SYS_DIR` ortam değişkeni, o da yoksa proje
kökündeki `sys/` klasörü. Dosyalar önce `.tmp` olarak yazılıp sonra yerine
konur, böylece yarım kalan indirme bozuk dosya bırakmaz.

```
sys/
├── ncbi/      genes/ nuccore/ protein/ region/ symbols/   (gen başına, istek üzerine)
├── gencode/   human/v50/ ...                               (sürüm başına, bir kez)
└── gtrnadb/   hg38-tRNAs.fa                                (tAI için, elle indirilir)
```

Dönen her nesnenin `features`'ına `source` (`'NCBI'` / `'GENCODE'`) ve kaynağa
özgü kimlikler yazılır. Transkriptlerde `utr5`, `cds` ve `utr3` 1 tabanlı,
iki ucu dahil `(başlangıç, bitiş)` koordinatlarıdır.

### 5.2 NcbiSession

`NcbiSession(api_key=None, e_mail=None, timeout=30, max_retries=5, pool_connections=10, pool_maxsize=20)`

NCBI E-utilities'e giden bütün istekleri taşır. NCBI'ın hız sınırına uyar:
anahtarsız saniyede 3, API anahtarıyla 10 istek. Geçici hatalarda yeniden
dener. `Ncbi` kendi oturumunu kendisi açar. Oturumu yalnızca API anahtarı ya
da e-posta vermek istediğinizde kurmanız gerekir.

| Metot | Ne yapar |
|---|---|
| `get(url, **kwargs)` → `requests.Response` | Hız sınırlı GET |
| `is_open()` → bool | Oturum açık mı |
| `close()` | Oturumu kapatır |

```python
>>> import os
>>> s = NcbiSession(api_key=os.environ.get("NCBI_API_KEY"))
>>> s.is_open()
True
>>> s.close(); s.is_open()
False
```

API anahtarını koda yazmayın, ortam değişkeninde tutun.

### 5.3 Ncbi

`Ncbi(organism="human", session=None, sys_dir=None, verbose=True)`

`organism` değeri `human`, `mouse`, `rat`, `zebrafish`, `fruit_fly` ya da `worm`
olabilir. `verbose=False` ilerleme çubuklarını kapatır. Bir kayıt `sys/`
içinde yoksa ilk istekte indirilir, sonraki isteklerde yerelden okunur.

```python
>>> ncbi = Ncbi(verbose=False)
>>> ncbi.ncbi_dir
PosixPath('/mnt/ssd/Programs/crn/sys/ncbi')
```

#### search_gene(gene) → str

Gen sembolünü NCBI Gene ID'ye çevirir. Rakamlardan oluşan değer zaten Gene ID
kabul edilir. Sonuç saklanır.

```python
>>> ncbi.search_gene("MEFV"), ncbi.search_gene("4210")
('4210', '4210')
```

#### gene(gene) → Gene

Bir genin genomik DNA'sını, gene bağlı **tüm** RefSeq transkriptlerini ve
proteinlerini tek bir `Gene` içinde toplar. En sık kullanılan metottur.

```python
>>> mefv = ncbi.gene("MEFV")
>>> mefv
MEFV: DNA <14607>, 2 transkript, 2 protein
>>> list(mefv.transcripts), list(mefv.proteins)
(['NM_000243.3', 'NM_001198536.2'], ['NP_000234.1', 'NP_001185465.2'])
>>> mefv.features.to_dict()
{'source': ['NCBI'], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
 'description': ['MEFV innate immunity regulator, pyrin'], 'organism': ['Homo sapiens'],
 'chromosome': ['16'], 'map_location': ['16p13.3'],
 'aliases': ['FMF', 'MEF', 'PAAND', 'TRIM20'],
 'genomic_region': [('NC_000016.10', 3242027, 3256633, '-')]}
```

#### dna(gene) → DNA

Genin genomik DNA'sını **kodlayan (sense) zincirde**, 5'→3' döndürür. Gen (−)
zincirdeyse referansın ters tümleyicisi alınır.

```python
>>> d = ncbi.dna("MEFV")
>>> d.features.to_dict()
{'source': ['NCBI'], 'chromosome': ['NC_000016.10'], 'start': [3242027], 'end': [3256633],
 'chromosome_strand': ['-'], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
 'description': ['MEFV innate immunity regulator, pyrin']}
```

#### rna(accession, genomic_strand=Unknown) → RNA

Bir RefSeq transkriptini olgun mRNA dizisi olarak döndürür. Sürüm eki
yazılmasa da olur (`NM_000243` → `NM_000243.3`).

```python
>>> r = ncbi.rna("NM_000243")
>>> r.features.to_dict()
{'source': ['NCBI'], 'transcript_id': ['NM_000243.3'], 'title': ['Homo sapiens MEFV ... mRNA'],
 'predicted': [False], 'gene_id': ['4210'], 'gene_name': ['MEFV'], 'protein_id': ['NP_000234.1'],
 'utr5': [(1, 46)], 'cds': [(47, 2392)], 'utr3': [(2393, 3506)]}
>>> start, end = r.features.first("cds")
>>> cds = r[start - 1:end]                       # 1 tabanlı -> Python dilimi
>>> len(cds), cds.seq[:12], cds.seq[-3:]
(2346, 'AUGGCUAAGACC', 'UGA')
```

#### protein(accession) → PRT

```python
>>> ncbi.protein("NP_000234").features.to_dict()
{'source': ['NCBI'], 'protein_id': ['NP_000234.1'], 'title': ['pyrin isoform 1 [Homo sapiens]']}
```

#### region(accession, start, end, strand="+", name="") → DNA

Kromozomdan ya da contig'den 1 tabanlı, iki ucu dahil bir bölge çıkarır.
`strand="-"` ise ters tümleyici alınır. Sonuç her durumda 5'→3' yazılır.

```python
>>> ncbi.region("NC_000016.10", 3_243_000, 3_243_020, "-", name="bolge")
bolge: TACTGTGCCCACGGAGCCTAC <21> [DNA] 5'->3'
```

#### update(gene, force=False) → str

Genin bütün kayıtlarını (gen özeti, genomik bölge, transkriptler, proteinler)
indirip `sys/` içine yazar ve Gene ID'yi döndürür. `gene()` bunu gerektiğinde
kendisi çağırır. Elle çağırmak yalnızca `force=True` ile kayıtları tazelemek
için gerekir. *(Ağ gerekir.)*

```python
ncbi.update("MEFV", force=True)      # '4210'
```

### 5.4 Gencode

`Gencode(species="human", release=None, sys_dir=None, base_url="https://ftp.ebi.ac.uk/pub/databases/gencode", verbose=True)`

`species` değeri `human` ya da `mouse` olabilir. `release=None` ise `sys/`
içindeki en yeni hazır sürüm kullanılır. GENCODE bir sürümün tamamını bir kez
indirir (~1,2 GB, açılmış hali ~6 GB) ve indeksler. Sonrasında bütün işlemler
çevrimdışıdır ve hızlıdır.

#### Sürüm yönetimi

| Üye | Ne verir |
|---|---|
| `release` | Kullanılan sürüm, ör. `'v50'` |
| `installed_releases()` | `sys/` içinde hazır sürümler |
| `is_ready()` | Kullanılan sürüm hazır mı |
| `latest_release()` | Sunucudaki en son sürüm *(ağ gerekir)* |
| `update(force=False)` | Sürümü indirip hazırlar, sürüm adını döndürür *(ağ gerekir)* |
| `species_dir`, `release_dir` | Yerel klasörler |

```python
>>> gc = Gencode(verbose=False)
>>> gc.release, gc.installed_releases(), gc.is_ready()
('v50', ['v50'], True)

# İlk kurulum (bir kez, ağ gerekir):
Gencode(release="v50").update()      # 'v50'
```

#### gene(gene) → Gene

```python
>>> mefv = gc.gene("MEFV")
>>> mefv
MEFV: DNA <16563>, 17 transkript, 15 protein
>>> mefv.features.to_dict()
{'source': ['GENCODE'], 'release': ['v50'], 'gene_id': ['ENSG00000103313.16'], 'gene_name': ['MEFV'],
 'gene_type': ['protein_coding'], 'chromosome': ['chr16'], 'start': [3240392], 'end': [3256954],
 'chromosome_strand': ['-']}
```

Aynı gen için NCBI 2, GENCODE 17 transkript verir. Kaynaklar farklı gen
sınırları ve transkript kümeleri kullanır. **Bir analizde tek kaynak
kullanın.**

#### dna(gene), rna(transcript), protein(identifier), region(chromosome, start, end, strand, name)

`Ncbi`'deki karşılıklarıyla aynı şekilde çalışır. Kimlikler ENSG, ENST ve
ENSP'dir. `protein()` ENSP ya da ENST kabul eder. `region()` kromozom adı
alır (`"chr16"`).

```python
>>> gc.rna("ENST00000339854.9").features.to_dict()
{'source': ['GENCODE'], 'release': ['v50'], 'transcript_id': ['ENST00000339854.9'],
 'transcript_name': ['MEFV-202'], 'transcript_type': ['protein_coding'],
 'gene_id': ['ENSG00000103313.16'], 'gene_name': ['MEFV'], 'chromosome': ['chr16'],
 'start': [3240392], 'end': [3256954], 'utr5': [(1, 367)], 'cds': [(368, 2173)],
 'utr3': [(2174, 4922)]}
>>> gc.protein("ENSP00000339639.4").features.first("transcript_id")
'ENST00000339854.9'
>>> gc.region("chr16", 3_243_000, 3_243_020, "+", name="b")
b: GTAGGCTCCGTGGGCACAGTA <21> [DNA] 5'->3'
```

(NCBI'daki `-` zincir bölgesi bunun ters tümleyicisidir:
`TACTGTGCCCACGGAGCCTAC`.)

---

## 6. cub — kodon kullanım yanlılığı (ayrıntılı başvuru)

```python
from crn.cub import (CUB, CubWarning, Composition, Synonymous, Adaptation,
                     EffectiveCodons, Divergence, CodonContext, read_trna_counts)
```

Bu bölümdeki örnekler aşağıdaki hazırlığı kullanır. `high` yüksek ifadeli
(ribozomal + EEF1A1), `low` düşük ifadeli (koku reseptörü) genlerdir.

```python
ncbi = Ncbi(verbose=False)
cub = CUB(ncbi.gene("MEFV").isoforms())
n = "MEFV_NM_000243.3"

def first_nm(symbol):
    return next(g for g in ncbi.gene(symbol).isoforms()
                if g.features.first("transcript_id").startswith("NM_"))

high = CUB([first_nm(s) for s in ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
                                   "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")])
low = CUB([first_nm(s) for s in ("OR1A1", "OR2W1", "OR5A1", "OR7A5", "OR1D2", "OR2J2",
                                  "OR4D1", "OR6A2", "OR10A2", "OR51E2", "OR2T4", "OR1E1")])
```

**Ortak kural:** Bütün ölçüm metotları `gene` parametresi alır. `gene=None`
(varsayılan) ise tüm genler **birlikte** (havuzlanmış sayımlar), gen adı
verilirse yalnızca o gen hesaplanır. Geçerli adlar `cub.names` içindedir.

### 6.1 A — CUB

`CUB(genes, genetic_code=GeneticCode.STANDARD)`

`genes` bir `Gene` listesi ya da `.crn` dosya yollarıdır (karışık da olabilir).

Kurulumda şunlar yapılır:

1. Her genin tam olarak 1 DNA + 1 transkript + 1 protein içerdiği denetlenir.
   Değilse hiçbir hesap yapılmadan `ValueError` verilir.
2. Transkriptte çevirisi proteinle aynı olan ORF bulunur (`match_orf`) ve
   hesaba yalnızca bu **CDS** girer.
3. Kodonlar gen başına bir kez sayılır. Bütün bölümler bu sayımı paylaşır.

Verilen genler değiştirilmez. `cub.genes` CDS'e indirgenmiş **kopyaları**
verir.

```python
>>> cub
CUB: 2 gen, L=1226
>>> cub.names, cub.genetic_code
(['MEFV_NM_000243.3', 'MEFV_NM_001198536.2'], <GeneticCode.STANDARD: NCBI 1, Standard>)
```

#### Birleştirme: `a + b`

İki CUB'u yeniden saymadan birleştirir. Genetik kodlar aynı olmalı ve iki tarafta aynı adlı gen bulunmamalıdır (yoksa `ValueError`). Panel
kurmak ve referans setlerini birleştirmek için kullanılır.

```python
>>> high, low, high + low
(CUB: 12 gen, L=3520, CUB: 12 gen, L=3781, CUB: 24 gen, L=7301)
```

#### Diziler ve kodon kümeleri

| Metot | Ne verir |
|---|---|
| `cds(gene)` | Genin hesaba giren CDS'i (RNA, 5'→3', stop dahil) |
| `codons(gene=None)` | Okuma sırasıyla kodon listesi (stop dahil) |
| `sense_codons(codons=CODONS)` | Amino asit kodlayanlar (61) |
| `codon_families(amino_acids=AMINO_ACIDS)` | `{aa: [kodonlar]}` |
| `synonymous_codons(codons=CODONS)` | Eş anlamlı kodonlar (59; Met, Trp ve stop hariç) |
| `fourfold_codons(codons=CODONS)` | Dört katlı dejenere kodonlar (32) |

```python
>>> cub.cds(n)[:15], len(cub.cds(n))
('AUGGCUAAGACCCCU', 2346)
>>> cub.codons(n)[:5], len(cub.codons(n)), len(cub.codons())
(['AUG', 'GCU', 'AAG', 'ACC', 'CCU'], 782, 1228)
>>> cub.codon_families()["S"]
['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC']
>>> len(cub.sense_codons()), len(cub.synonymous_codons()), len(cub.fourfold_codons())
(61, 59, 32)
```

#### Sayımlar ve frekanslar (ölçüm 1–5)

| Metot | Formül / anlam |
|---|---|
| `codon_counts(gene, codons=CODONS)` | 1. Her kodonun sayısı (stop dahil) |
| `total_codons(gene)` | 2. L: stop hariç kodon sayısı |
| `amino_acid_counts(gene, amino_acids)` | 3. Amino asit başına kodon sayısı |
| `codon_freqs(gene)` | 4. N_kodon / L (anlamlı kodonlar, toplam 1) |
| `codon_per_thousand(gene)` | 5. Frekans × 1000 (kodon kullanım tablolarının olağan birimi) |
| `termini(gene)` | Her proteinin (ilk, son) kalıntısı |
| `chain_count(gene)` | Hesaba giren protein zinciri sayısı |
| `codon_pair_counts(gene)` | Komşu kodon çiftleri, `"GCCGAA"` biçiminde. Stop/belirsiz içerenler hariç, önbelleğe alınır |

```python
>>> cc = cub.codon_counts(n)
>>> cc["CUG"], cc["UGA"], cub.codon_counts(n, codons=("AUG", "UGG"))
(52, 1, {'AUG': 9, 'UGG': 6})
>>> cub.total_codons(n), cub.total_codons()
(781, 1226)
>>> cub.amino_acid_counts(n, "LKW")
{'L': 71, 'K': 49, 'W': 6}
>>> round(cub.codon_freqs(n)["CUG"], 4), round(cub.codon_per_thousand(n)["CUG"], 2)
(0.0666, 66.58)
>>> cub.termini(), cub.chain_count(), cub.chain_count(n)
([('M', 'D'), ('M', 'C')], 2, 1)
>>> pairs = cub.codon_pair_counts(n)
>>> len(pairs), sorted(pairs.items(), key=lambda x: -x[1])[:3]
(604, [('CUGGAG', 6), ('CUGAAG', 6), ('GAGGAG', 5)])
```

#### table(**kwargs) → list[dict]

Her gen için bütün skaler ölçümleri tek satırda verir. Parametreleri
`summary_table` ile aynıdır.

| Parametre | Eklenen sütunlar |
|---|---|
| (her zaman) | `gene, L, GC%, GC1, GC2, GC12, GC3, GC3s, ENC, ENC_expected, ENC_ratio, SCUO, CDC, RCBS` |
| `reference=` | `Fop, CAI, CBI` |
| `background=` | `Fop_chi2, CBI_chi2` (χ² optimal kodonlarla) |
| `host=` | `RCDI` |
| `genome=` | `CPB` |
| `trna_counts=`, `kingdom=` | `tAI` |
| `protein=True` | `length, MW, pI, net_charge_7, GRAVY, aromaticity, aliphatic_index, extinction_280, IVYWREL, CvP, TOP_IDP` |
| `file_name=` | Sonucu ayrıca TSV olarak yazar |

`CubWarning`'ler tablo kurulurken susturulur.

```python
>>> rows = cub.table()
>>> list(rows[0])
['gene', 'L', 'GC%', 'GC1', 'GC2', 'GC12', 'GC3', 'GC3s', 'ENC', 'ENC_expected', 'ENC_ratio',
 'SCUO', 'CDC', 'RCBS', 'length', 'MW', 'pI', 'net_charge_7', 'GRAVY', 'aromaticity',
 'aliphatic_index', 'extinction_280', 'IVYWREL', 'CvP', 'TOP_IDP']
```

Tam bir panel tablosu §9.2'dedir. pandas ile `pd.DataFrame(rows)` doğrudan
çalışır.

### 6.2 B — Composition

`Composition(cub)`: nükleotit bileşimi. Mutasyon baskısını (GC zenginliği)
ölçer. Üçüncü konum (GC3, GC3s) sinonim seçilime en açık konumdur.

| Metot | Ölçüm |
|---|---|
| `nucleotide_counts(gene, position=None, codons=CODONS, bases="UCAG")` | Anlamlı kodonlarda baz sayıları; `position` 1/2/3 ya da tümü |
| `nucleotide_freqs(...)` | Aynısının oranları (toplam 1) |
| `gc_percent(gene)` | 6. GC% |
| `gc1`, `gc2`, `gc3` | 7, 8, 10. Konuma göre G+C oranı |
| `gc12(gene)` | 9. (GC1 + GC2) / 2 |
| `gc3s(gene)` | 11. Eş anlamlı kodonların 3. konumunda G+C |
| `base3s(gene)` | 12. A3s, U3s, G3s, C3s |
| `at3s(gene)` | 13. 1 − GC3s |
| `summary(gene)` | 6–13 tek sözlükte |

```python
>>> comp = Composition(cub)
>>> comp.nucleotide_counts(n, position=3)
{'U': 121, 'C': 248, 'A': 113, 'G': 299}
>>> comp.gc3s(n), comp.base3s(n)
(0.6945, {'U': 0.158, 'C': 0.3238, 'A': 0.1475, 'G': 0.3708})
>>> comp.summary(n)
{'GC%': 59.0269, 'GC1': 0.6133, 'GC2': 0.4571, 'GC12': 0.5352, 'GC3': 0.7004, 'GC3s': 0.6945,
 'A3s': 0.1475, 'U3s': 0.158, 'G3s': 0.3708, 'C3s': 0.3238, 'AT3s': 0.3055}
```

(Değerler okunabilirlik için 4 basamağa yuvarlanmıştır.)

### 6.3 C — Synonymous

`Synonymous(cub, reference=None, background=None)`: eş anlamlı kodon tercihi.
`reference` optimal kodonların belirlendiği set, `background` χ² yöntemi için
karşılaştırma setidir. `reference` verilmezse `cub`'ın kendisi kullanılır.

| Metot | Ölçüm |
|---|---|
| `codon_ratio(gene)` | 14. Kodonun kendi amino asidi içindeki payı (aile toplamı 1) |
| `rscu(gene)` | 15. RSCU = pay × aile büyüklüğü. 1 tarafsız, >1 tercih |
| `preferred_codons(gene, threshold=1.6)` | 16. RSCU > 1,6 olan kodonlar |
| `underused_codons(gene, threshold=0.6)` | 16. RSCU < 0,6 olanlar (hiç kullanılmayan amino asitler hariç) |
| `optimal_codons(method="rscu", threshold=1.6, alpha=0.01)` | Referansın optimal kodonları. `"rscu"`: referans RSCU > eşik. `"chi2"`: referans ve background arasında 2×2 χ², p < alpha. Önbelleğe alınır |
| `fop(gene, optimal_codons=None, method, denominator="synonymous", ...)` | 17. Fop = optimal kodon sayısı / payda. `denominator="optimal"`: yalnızca optimal kodonu olan amino asitler (CodonW) |

```python
>>> syn = Synonymous(cub, reference=high, background=low)
>>> round(syn.codon_ratio(n)["CUG"], 4)
0.7324
>>> {c: round(syn.rscu(n)[c], 3) for c in cub.codon_families()["L"]}
{'UUA': 0.0, 'UUG': 0.085, 'CUU': 0.338, 'CUC': 0.845, 'CUA': 0.338, 'CUG': 4.394}
>>> syn.preferred_codons(n)
['CUG', 'CCC', 'CAG', 'AUC', 'AAG', 'AGC', 'AGG', 'GUG', 'GGG']
>>> syn.underused_codons(n)
['UUU', 'UUA', 'UUG', 'UCG', 'UGU', 'CUU', 'CUA', 'CAA', 'CGU', 'CGA', 'AUA', 'ACG', 'AAA', 'GUA', 'GGU']
>>> syn.optimal_codons()                       # high setinin RSCU > 1,6 kodonları
['UCU', 'CUG', 'ACC', 'GUG', 'GCU', 'GGC']
>>> syn.optimal_codons(method="chi2")          # high'da low'a göre anlamlı fazla
['GCU', 'GAU']
>>> syn.fop(n), syn.fop(n, method="chi2"), syn.fop(n, denominator="optimal")
(0.1632, 0.0274, 0.3776)
```

**Yorum:** MEFV'de Leu'nun %73'ü CUG ile kodlanıyor (RSCU 4,39). Bu,
GC3 zenginliğiyle uyumlu güçlü bir tercihtir.

### 6.4 D — Adaptation

`Adaptation(cub, reference=None, background=None)`: genin yüksek ifadeli
referans sete ne kadar benzediği. Translasyon verimliliğinin vekil ölçüsüdür.

| Metot | Ölçüm |
|---|---|
| `relative_adaptiveness(absent_count=0.5)` | 18. w: referansta kodon sayısı / ailesindeki en büyük sayı. Referansta olmayan kodon için sayı `absent_count` alınır. Önbelleğe alınır |
| `cai(gene, w=None, absent_count=0.5)` | 19. CAI: w'lerin geometrik ortalaması (0–1). `w` dışarıdan da verilebilir |
| `cbi(gene, optimal_codons=None, method="rscu", ...)` | 20. CBI: (N_opt − N_rastgele) / (N_toplam − N_rastgele). 0 rastgele, 1 tam tercih, negatif ters tercih |

```python
>>> ad = Adaptation(cub, reference=high, background=low)
>>> w = ad.relative_adaptiveness()
>>> len(w), {c: round(w[c], 3) for c in cub.codon_families()["L"]}
(59, {'UUA': 0.068, 'UUG': 0.447, 'CUU': 0.427, 'CUC': 0.398, 'CUA': 0.117, 'CUG': 1.0})
>>> ad.cai(n), ad.cbi(n), ad.cbi(n, method="chi2")
(0.7254, 0.2069, -0.1454)
>>> Adaptation(cub).cai(n)            # referans verilmezse gen kendisiyle kıyaslanır
0.6962
```

### 6.5 F — EffectiveCodons

`EffectiveCodons(cub)`: etkin kodon sayısı (ENC/Nc) ve CUB'un mutasyon mu
seçilim mi kaynaklı olduğunu ayırt eden grafik analizleri. Nötrallik ve PR2
çok genli panel ister.

| Metot | Ölçüm |
|---|---|
| `homozygosity(gene)` | 23. Her amino asidin F değeri: (nΣp² − 1)/(n − 1) |
| `enc(gene)` | 23. Nc (Wright 1990): 20 aşırı yanlılık, 61 yanlılık yok. Bir dejenerelik sınıfı boşsa `nan` |
| `EffectiveCodons.expected_enc(gc3s)` | 24. Statik: yalnızca GC3s'ten beklenen Nc = 2 + s + 29/(s² + (1−s)²) |
| `enc_ratio(gene)` | 25. (beklenen − gözlenen) / beklenen. 0'dan uzaklık seçilim izi |
| `enc_points()` | 24–25. Grafik verisi: gen başına gc3s, enc, expected_enc, enc_ratio |
| `neutrality()` | 26. GC12 ~ GC3 regresyonu: `{'points', 'slope', 'intercept', 'r2'}`. İki genden az ise `nan` |
| `pr2(gene, codons=None)` | 27. x = G3/(G3+C3), y = A3/(A3+U3). (0,5; 0,5) tarafsız |
| `pr2_points(codons=None)` | 27. Gen başına (x, y) |
| `plot_enc(ax=None, labels=False, file_name=None)` | ENC–GC3s grafiği ve Wright eğrisi |
| `plot_enc_ratio(ax=None, bin_width=0.05, file_name=None)` | ENC oranı histogramı |
| `plot_neutrality(ax=None, labels=False, file_name=None)` | Nötrallik grafiği, regresyon doğrusu ve y = x |
| `plot_pr2(codons=None, ax=None, labels=False, file_name=None)` | PR2 grafiği |

Grafik metotları bir matplotlib `Axes` döndürür. `ax=` ile var olan bir
şekle çizer, `file_name=` ile dosyaya kaydeder.

```python
>>> panel = CUB([cub.genes[0]]) + high + low          # 25 gen
>>> eff = EffectiveCodons(panel)
>>> round(eff.homozygosity(n)["L"], 4), round(eff.enc(n), 2)
(0.5565, 45.44)
>>> round(EffectiveCodons.expected_enc(0.6945), 2), round(eff.enc_ratio(n), 4)
(53.07, 0.1438)
>>> eff.enc_points()[0]
{'gene': 'MEFV_NM_000243.3', 'gc3s': 0.6945, 'enc': 45.4365, 'expected_enc': 53.0703, 'enc_ratio': 0.1438}
>>> ne = eff.neutrality()
>>> ne["slope"], ne["intercept"], ne["r2"]
(0.1288, 0.3757, 0.0709)
>>> eff.pr2(n)
{'x': 0.5303, 'y': 0.5}
>>> eff.plot_enc(file_name="enc.png")
<Axes: ...>
```

**Yorum:** MEFV'nin ENC'si (45,4) GC3s'ten beklenenin (53,1) altında kalıyor
(oran 0,14). Yanlılığın bir kısmı yalnızca GC3 ile açıklanamıyor. Panelde
nötrallik eğimi 0,13'tür: GC3 değişiminin ancak ~%13'ü mutasyon baskısıyla
açıklanıyor, seçilim baskın.

### 6.6 G — Divergence

`Divergence(cub, reference=None)`: kodon kullanımının düzensizliği ve
beklenenden ya da referanstan sapması.

| Metot | Ölçüm |
|---|---|
| `expected_codon_freqs(gene)` | Konum baz oranlarından beklenen kodon oranları f1·f2·f3 |
| `scuo(gene)` | 28. SCUO: entropiye dayalı düzensizlik (0 tarafsız, 1 tam yanlı) |
| `cdc(gene)` | 29. CDC: gözlenen ve beklenen kodon oranları arasındaki kosinüs uzaklığı (0–1) |
| `rcdi(gene)` | 30. RCDI: referansa (konak) göre sapma. 1 referansla aynı, büyüdükçe uyumsuz. Referansta olmayan kodonlar `CubWarning` ile çıkarılır |
| `rcbs(gene)` | 31. RCBS: beklenene göre göreli yanlılık |

```python
>>> dv = Divergence(cub, reference=high + low)
>>> e = dv.expected_codon_freqs(n)
>>> round(e["CUG"], 4), round(sum(e.values()), 4)
(0.0249, 1.0)
>>> dv.scuo(n), dv.cdc(n), dv.rcdi(n), dv.rcbs(n)
(0.1713, 0.1175, 1.2799, 0.1668)
```

### 6.7 H — CodonContext

`CodonContext(cub, reference=None)`: tRNA havuzuna uyum ve kodon bağlamı.
`reference` CPS skorlarının hesaplandığı geniş set olmalıdır (genom ya da
büyük bir gen paneli).

#### read_trna_counts(file_name) → dict

GtRNAdb FASTA dosyasından antikodon başına tRNA gen kopya sayısını (tGCN)
okur. Antikodonlar 5'→3' RNA alfabesiyle döner. Psödogen ve SeC gibi kodon
okumayan tRNA'lar atlanır. Dosya
[gtrnadb.ucsc.edu](https://gtrnadb.ucsc.edu)'dan (`hg38-tRNAs.fa`) elle
indirilir.

```python
>>> trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
>>> len(trna), sum(trna.values())
(47, 431)
>>> {k: trna[k] for k in list(trna)[:4]}
{'AGC': 26, 'CGC': 4, 'UGC': 8, 'ACG': 7}
```

| Metot | Ölçüm |
|---|---|
| `trna_weights(trna_counts, s_values=TAI_S_VALUES, kingdom="eukaryote")` | 32. Her kodonun tAI ağırlığı w = W/W_max, W = Σ(1 − s)·tGCN (wobble eşleşmeleri dahil). `kingdom="prokaryote"` AUA için lizidin kuralını uygular |
| `tai(gene, trna_counts, s_values, kingdom)` | 32. tAI: w'lerin geometrik ortalaması (0–1) |
| `dinucleotide_counts(gene, frame=None, bases="UCAG")` | 33. Dinükleotit sayıları. `frame`: `None` (tüm dizi), `"12"`, `"23"`, `"31"` (kodon sınırı: 3. baz + sonraki 1. baz) |
| `rho(gene, frame=None, bases="UCAG")` | 33. ρ(XY) = f(XY)/(f(X)·f(Y)). <0,78 az temsil, >1,23 fazla temsil |
| `codon_pair_scores()` | 34. Referanstaki her çiftin CPS'si (log oran). Önbelleğe alınır |
| `cpb(gene)` | 34. CPB: genin çiftlerinin CPS ortalaması. Referansta olmayan çiftler `CubWarning` ile çıkarılır |

```python
>>> ctx = CodonContext(cub, reference=high + low)
>>> w = ctx.trna_weights(trna)
>>> max(w, key=w.get), min(w, key=w.get), round(w["CUG"], 4)
('UGC', 'CUA', 0.3434)
>>> ctx.tai(n, trna_counts=trna), ctx.tai(n, trna_counts=trna, kingdom="prokaryote")
(0.3386, 0.3392)
>>> dc = ctx.dinucleotide_counts(n)
>>> dc["CG"], dc["GC"], ctx.dinucleotide_counts(n, frame="31")["CG"]
(95, 201, 37)
>>> round(ctx.rho(n)["CG"], 3), round(ctx.rho(n, frame="31")["CG"], 3)
(0.467, 0.456)
>>> cps = ctx.codon_pair_scores()
>>> len(cps), sorted(cps.items(), key=lambda x: x[1])[:2]
(2514, [('GCCGCC', -2.1414), ('CCUAUG', -2.026)])
>>> ctx.cpb(n)
CubWarning: CPB: referansta hiç geçmeyen 119 kodon çifti hesaptan çıkarıldı
0.1649
```

**Yorum:** ρ(CG) = 0,47 değeri, omurgalılarda metilasyona bağlı tipik CpG
baskılanmasını gösteriyor. CPB'nin pozitif olması, MEFV'nin referansta
fazla temsil edilen kodon çiftlerini kullandığını söylüyor.

### 6.8 summary_table

`summary_table(cub, reference=None, background=None, host=None, genome=None, trna_counts=None, kingdom="eukaryote", protein=True, file_name=None)`

`cub.table(...)`'ın arkasındaki fonksiyondur (`crn.cub.report`). Parametreler
§6.1'deki tablodaki gibidir.

---

## 7. ProteinProperties — protein özellikleri (ayrıntılı başvuru)

`ProteinProperties(source)`, proteinlerin amino asit bileşiminden hesaplanan
fizikokimyasal özelliklerini verir. Literatürdeki ProtParam, EMBOSS pepstats
ve düzensizlik araçlarının ölçümlerini tek sınıfta toplar.

`source` üç biçimde olabilir:

| Kaynak | Proteinler | Gen adları |
|---|---|---|
| `CUB` | CDS'lerin çevirisi (CUB ile aynı genler) | `cub.names` |
| `PRT` | O protein | `[prt.name]` |
| `[PRT, ...]` | Liste | adları |

İki yoldan aynı protein için aynı sonuç çıkar. `gene=None` ise bütün
proteinler **ayrı zincirler** olarak birlikte hesaplanır: kütleye zincir
başına bir H2O eklenir, uç grupları ve disülfit eşleri zincir başına sayılır.
Stop (`*`) hesaba girmez.

```python
>>> props = ProteinProperties(cub)                       # §6'daki MEFV CUB'u
>>> props.names
['MEFV_NM_000243.3', 'MEFV_NM_001198536.2']
>>> sp = ProteinProperties(PRT("MKWVTFISLLLLFSSAYS", name="sinyal"))
>>> sp.names, round(sp.gravy(), 3), round(sp.molecular_weight(), 2)
(['sinyal'], 1.289, 2106.55)
>>> lst = ProteinProperties([PRT("MKWVTF", name="a"), PRT("MDEEK", name="b")])
>>> lst.chain_count(), lst.length(), lst.length("b")
(2, 11, 5)
```

Aşağıdaki bütün örnekler `n = "MEFV_NM_000243.3"` (pyrin, 781 aa) içindir.

### 7.1 E-A — Sayım ve kütle

| Metot | Ölçüm |
|---|---|
| `chain_count(gene)` | Zincir sayısı |
| `length(gene, amino_acids)` | 1. N: kalıntı sayısı |
| `composition(gene)` | 2. Her amino asidin yüzdesi |
| `group_composition(gene, groups=AMINO_ACID_GROUPS)` | 3. Grup oranları (0–1) |
| `scale_average(scale, gene)` | Herhangi bir kalıntı ölçeğinin bileşim ağırlıklı ortalaması (genel araç) |
| `molecular_weight(gene, masses, water)` | 4. Ortalama MW (Da) |
| `monoisotopic_mass(gene, ...)` | 5. Monoizotopik kütle (kütle spektrometrisi için) |
| `mean_residue_mass(gene)` | 6. MW / N |
| `residue_counts(gene, residues="CWY")` | 7. Seçili kalıntı sayıları |
| `disulfide_pairs(gene)` | 7. Olası disülfit bağı: ⌊nCys/2⌋ |

```python
>>> props.length(n), round(props.composition(n)["L"], 2)
(781, 9.09)
>>> props.group_composition(n)
{'hydrophobic': 0.501, 'polar': 0.499, 'aliphatic': 0.254, 'aromatic': 0.059,
 'small': 0.443, 'basic': 0.157, 'acidic': 0.127}
>>> props.scale_average(HYDROPATHY_SCALES["Kyte-Doolittle"], n)    # = GRAVY
-0.645
>>> props.molecular_weight(n), props.monoisotopic_mass(n), props.mean_residue_mass(n)
(86443.77, 86389.98, 110.68)
>>> props.residue_counts(n), props.disulfide_pairs(n)
({'C': 17, 'W': 6, 'Y': 19}, 8)
```

### 7.2 E-B — Atom bileşimi

| Metot | Ölçüm |
|---|---|
| `formula(gene, ...)` | 8. `{C, H, N, O, S}` atom sayıları |
| `formula_string(gene)` | 8. Metin olarak |
| `atom_count(gene)` | 9. Toplam atom |
| `atoms_per_residue(gene)` | 10. Kalıntı başına atom |
| `carbon_oxidation_state(gene, charge=0.0)` | 11. Ortalama karbon oksidasyon durumu (Dick) |
| `termini(gene)` | Zincir uçları (ilk, son kalıntı) |

```python
>>> props.formula(n)
{'C': 3779, 'H': 6038, 'N': 1104, 'O': 1167, 'S': 26}
>>> props.formula_string(n), props.atom_count(n)
('C3779H6038N1104O1167S26', 12114)
>>> round(props.atoms_per_residue(n), 3), round(props.carbon_oxidation_state(n), 3), props.termini(n)
(15.511, -0.09, [('M', 'D')])
```

### 7.3 E-C — Yük

| Metot | Ölçüm |
|---|---|
| `charged_counts(gene)` | 12. `{positive: R+K, negative: D+E}` |
| `simple_net_charge(gene)` | 13. Kaba net yük (R+K) − (D+E) |
| `fcr(gene)` | 14. Yüklü kalıntı oranı |
| `ncpr(gene)` | 15. Kalıntı başına net yük |
| `net_charge(gene, ph=7.0, pka="Bjellqvist")` | 16. Henderson–Hasselbalch ile Q(pH) (uçlar, His, Cys, Tyr dahil) |
| `isoelectric_point(gene, pka="Bjellqvist", tolerance, ph_range)` | 17. Q = 0 olan pH |
| `titration_curve(gene, pka, step=0.1, ph_range)` | 18. (pH, Q) çiftleri |

`pka` bir `PKA_SETS` adı (`'Bjellqvist'` ExPASy ile uyumludur, ayrıca
`'EMBOSS'` ve `'Lehninger'`) ya da kendi sözlüğünüz olabilir.

```python
>>> props.charged_counts(n), props.simple_net_charge(n)
({'positive': 104, 'negative': 99}, 5)
>>> round(props.fcr(n), 3), round(props.ncpr(n), 3)
(0.26, 0.006)
>>> round(props.net_charge(n), 2), round(props.net_charge(n, ph=5.0, pka="EMBOSS"), 2)
(6.15, 33.33)
>>> [round(props.isoelectric_point(n, pka=k), 2) for k in ("Bjellqvist", "EMBOSS", "Lehninger")]
[8.3, 8.11, 7.94]
>>> props.titration_curve(n, step=3.5)
[(0.0, 123.99), (3.5, 110.18), (7.0, 6.15), (10.5, -65.83), (14.0, -135.45)]
```

### 7.4 E-D — Hidropati

| Metot | Ölçüm |
|---|---|
| `gravy(gene)` | 19. Kyte–Doolittle ortalaması. Negatif hidrofilik |
| `hydropathy(gene, scale="Kyte-Doolittle")` | 20. Herhangi bir `HYDROPATHY_SCALES` ölçeği ya da sözlük |
| `aromaticity(gene)` | 21. (F+W+Y) / N |
| `aliphatic_index(gene)` | 22. Ikai alifatik indeksi (termostabilite göstergesi) |
| `boman_index(gene)` | 23. Boman indeksi. >2,48 yüksek bağlanma potansiyeli |

```python
>>> props.gravy(n)
-0.645
>>> {s: round(props.hydropathy(n, scale=s), 3) for s in HYDROPATHY_SCALES}
{'Kyte-Doolittle': -0.645, 'Hopp-Woods': 0.256, 'Eisenberg': -0.113, 'Wimley-White': 0.322, 'Janin': -0.213}
>>> props.aromaticity(n), props.aliphatic_index(n), props.boman_index(n)
(0.059, 73.2, 2.144)
```

### 7.5 E-E — Optik ve hacim

| Metot | Ölçüm |
|---|---|
| `extinction_coefficient(gene, reduced=False)` | 24. ε280 (M⁻¹cm⁻¹). `reduced=True` sistinleri saymaz |
| `absorbance_01(gene, reduced=False)` | 25. 1 mg/mL çözeltinin A280'i |
| `partial_specific_volume(gene)` | 26. v̄ (cm³/g), santrifüj ve SAXS için |
| `mean_residue_volume(gene)` | 27. Ortalama kalıntı hacmi (Å³) |
| `mean_residue_surface(gene)` | 27. Ortalama en büyük yüzey alanı (Å²) |

```python
>>> props.extinction_coefficient(n), props.extinction_coefficient(n, reduced=True)
(62310, 61310)
>>> round(props.absorbance_01(n), 3), round(props.partial_specific_volume(n), 4)
(0.721, 0.7312)
>>> round(props.mean_residue_volume(n), 2), round(props.mean_residue_surface(n), 2)
(131.73, 191.76)
```

### 7.6 E-F — Stabilite ve çözünürlük

| Metot | Ölçüm |
|---|---|
| `ivywrel(gene)` | 28. IVYWREL oranı (termofilik proteinlerde yüksek) |
| `cvp_bias(gene)` | 29. Yüklü − polar yanlılığı (%) |
| `solubility_cv(gene)` | 30. Wilkinson–Harrison kanonik değişkeni |
| `solubility_probability(gene, discriminant=1.71)` | 30. `{cv, soluble, probability}`: E. coli'de çözünür ifade tahmini |
| `composition_entropy(gene)` | 31. Bileşim entropisi (bit), en çok 4,32 |

```python
>>> props.ivywrel(n), props.cvp_bias(n), props.solubility_cv(n)
(0.373, 4.48, 3.392)
>>> props.solubility_probability(n)
{'cv': 3.392, 'soluble': False, 'probability': 0.847}
>>> props.composition_entropy(n)
4.104
```

`soluble=False, probability=0.847`: model, %84,7 olasılıkla **çözünmez**
(inclusion body) sınıfını öngörüyor.

### 7.7 E-G — Düzensizlik ve bileşim modelleri

| Metot | Ölçüm |
|---|---|
| `uversky(gene)` | 32. Uversky grafiği: ortalama net yük, ortalama hidropati, sınır, `disordered` |
| `top_idp(gene)` | 33. TOP-IDP ortalaması. Pozitif değer düzensizlik eğilimi |
| `composition_features(gene)` | 34. Bileşim tabanlı modellerin özellikleri: 20 aa oranı + `K-R`, `D-E`, `K+R`, `D+E`, `K+R-D-E`, `K+R+D+E`, `F+W+Y`, `pI`, `abs_charge`, `gravy`, `entropy`, `top_idp`, `fcr`, `ncpr` |
| `linear_score(gene, weights, intercept=0.0)` | 34. sabit + Σ ağırlık × özellik. Protein-Sol gibi doğrusal modelleri kendi ağırlıklarınızla uygular |
| `summary(gene)` | Başlıca skaler özellikler tek sözlükte (`cub.table()` sütunları) |

```python
>>> props.uversky(n)
{'mean_net_charge': 0.006, 'mean_hydropathy': 0.428, 'boundary': 0.416, 'disordered': False}
>>> props.top_idp(n), round(props.composition_features(n)["K-R"], 3)
(0.175, -0.008)
>>> round(props.linear_score(n, weights={"K": 1.0, "D-E": -0.5}, intercept=0.1), 4)
0.1864
>>> props.linear_score(n)
ValueError: weights verilmeli: composition_features() anahtarı -> regresyon katsayısı
>>> props.summary(n)
{'length': 781, 'MW': 86443.77, 'pI': 8.298, 'net_charge_7': 6.149, 'GRAVY': -0.645,
 'aromaticity': 0.059, 'aliphatic_index': 73.201, 'extinction_280': 62310,
 'IVYWREL': 0.373, 'CvP': 4.481, 'TOP_IDP': 0.175}
```

---

## 8. Hatalar ve uyarılar

crn belirsiz durumları sessizce geçmez. Ya hesabı durdurur (`ValueError`,
`KeyError`, `TypeError`) ya da veriyi dışarıda bıraktığını bildirir
(`UserWarning`, `CubWarning`).

| Mesaj | Neden | Çözüm |
|---|---|---|
| `DNA dizisinde geçersiz karakter(ler): U` | Alfabede olmayan harf | Doğru sınıfı kullanın (RNA'da U) ya da diziyi temizleyin |
| `genes bir liste olmalı, Gene verildi; tek gen için genes=[gene] yazın` | `CUB(gene)` | `CUB([gene])` |
| `Her gende bir genomik DNA, bir transkript ve bir protein olmalı (çok izoformlu genler için Gene.isoforms()); hesap durduruldu` | Gen izoformlara ayrılmamış | `CUB(gene.isoforms())` |
| `…: NM_… transkriptinde NP_… proteiniyle örtüşen ORF yok; hesap durduruldu` | Protein bu genetik kodla transkriptten çevrilemiyor | `genetic_code`'u ve kaydı kontrol edin |
| `CUB: 'yok' geni yok` (`KeyError`) | Yanlış gen adı | `cub.names` |
| `trna_counts verilmeli: …` | tAI için tRNA verisi yok | `read_trna_counts(...)` |
| `method='chi2' için reference ve background farklı gen setleri olmalı` | χ² için iki set gerekir | `Synonymous(cub, reference=high, background=low)` |
| `denominator 'synonymous', 'optimal' olmalı, … verildi` | Geçersiz seçenek | Listelenen değerlerden birini kullanın |
| `Farklı genetik kodlu CUB'lar birleştirilemez` | `a + b`'de kodlar farklı | Aynı `genetic_code` |
| `demo: 'p1' proteini zaten var` | Aynı kimlik iki kez eklendi | Farklı kimlik verin |
| `UserWarning: … aynı CDS'i veren izoformlar çıkarıldı (A -> B)` | `isoforms()` kopyaları eledi | Bilgi amaçlı. Kopyalar bilerek tutulmak isteniyorsa `unique_cds=False` |
| `UserWarning: … proteini olmayan transkriptler: NR_…` | Kodlamayan transkript | Bilgi amaçlı |
| `CubWarning: CPB: referansta hiç geçmeyen 119 kodon çifti hesaptan çıkarıldı` | Referans set küçük | Daha geniş referans. Susturmak için `warnings.filterwarnings("ignore", category=CubWarning)` |
| ENC `nan` | Bir dejenerelik sınıfında veri yok (çok kısa gen) | Daha uzun gen ya da havuzlanmış set (`gene=None`) |
| Nötrallik `nan` | İkiden az gen | Çok genli panel |

```python
>>> CUB(ncbi.gene("MEFV"))
TypeError: genes bir liste olmalı, Gene verildi; tek gen için genes=[gene] yazın
>>> CUB([Gene("x", DNA("ATG"))])
ValueError: Her gende bir genomik DNA, bir transkript ve bir protein olmalı (...); hesap durduruldu:
  x: 0 transkript; 0 protein
>>> CUB(isoforms, genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL)    # nükleer gen, yanlış kod
ValueError: MEFV_NM_000243.3: NM_000243.3 transkriptinde NP_000234.1 proteiniyle örtüşen ORF yok; hesap durduruldu
```

---

## 9. Uçtan uca örnekler

Aşağıdaki dört senaryo, çerçevenin parçalarını gerçek bir analiz akışında
birleştirir. Hepsi proje kökünden çalıştırılmıştır.

### 9.1 Santral dogma: DNA'dan proteine

Bir DNA parçasında ORF bulup mRNA'ya, oradan proteine geçmek ve proteini
tanımlamak:

```python
from crn import DNA, ProteinProperties

dna = DNA("CCATGGCTAAGCTGCTGGAGTGGTACTGATTT", name="demo")

orfs = dna.find_orfs(minimum_length_of_orf=9)
orf = max(orfs, key=len)
print(orf, orf.features.to_dict())
# demo_ORF1: ATGGCTAAGCTGCTGGAGTGGTACTGA <27> [DNA] 5'->3'
# {'orf_start': [3], 'orf_end': [29], 'orf_strand': ['+'], 'orf_frame': [2]}

mrna = orf.transcription(template=False)     # ORF kodlayan zincirdir
print(mrna)              # demo_ORF1: AUGGCUAAGCUGCUGGAGUGGUACUGA <27> [RNA] 5'->3' (?)

protein = mrna.translation()
print(protein)           # demo_ORF1: MAKLLEWY* <9> [PRT] N->C

props = ProteinProperties(protein)           # stop (*) hesaba girmez
print(props.length(), round(props.molecular_weight(), 2), props.formula_string())
# 8 1053.28 C51H76N10O12S1
print(round(props.isoelectric_point(), 2), round(props.gravy(), 3), props.extinction_coefficient())
# 5.75 0.212 6990

print(mrna.reverse_transcription())          # cDNA (ilk zincir)
# demo_ORF1: TACCGATTCGACGACCTCACCATGACT <27> [DNA] 3'->5'
```

### 9.2 Yüksek ve düşük ifadeli genlerle panel analizi (NCBI)

Soru: MEFV'nin kodon kullanımı yüksek ifadeli genlere mi, düşük ifadeli
genlere mi benziyor? Kodon kullanımını mutasyon baskısı mı, seçilim mi
belirliyor?

```python
import matplotlib.pyplot as plt
from crn.fetch import Ncbi
from crn.cub import CUB, EffectiveCodons, read_trna_counts

ncbi = Ncbi(verbose=False)

def first_nm(symbol):
    """Genin ilk NM_ izoformu."""
    return next(g for g in ncbi.gene(symbol).isoforms()
                if g.features.first("transcript_id").startswith("NM_"))

# 1. Gen setleri
mefv = CUB([first_nm("MEFV")])
high = CUB([first_nm(s) for s in ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
                                   "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")])
low = CUB([first_nm(s) for s in ("OR1A1", "OR2W1", "OR5A1", "OR7A5", "OR1D2", "OR2J2",
                                  "OR4D1", "OR6A2", "OR10A2", "OR51E2", "OR2T4", "OR1E1")])
panel = mefv + high + low                      # CUB: 25 gen; yeniden sayım yok

# 2. Tüm ölçümler tek tabloda
trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
rows = panel.table(reference=high, background=low,       # Fop, CAI, CBI (+ χ²)
                   host=high + low,                       # RCDI
                   genome=high + low,                     # CPB
                   trna_counts=trna,                      # tAI
                   file_name="panel.tsv")

# 3. Grafikler
eff = EffectiveCodons(panel)
fig, axs = plt.subplots(2, 2, figsize=(11, 10))
eff.plot_enc(ax=axs[0, 0])
eff.plot_enc_ratio(ax=axs[0, 1])
eff.plot_neutrality(ax=axs[1, 0], labels=True)
eff.plot_pr2(ax=axs[1, 1])
fig.savefig("f_grafikleri.png", dpi=300)
print(eff.neutrality()["slope"], eff.neutrality()["r2"])      # 0.1288 0.0709
```

`panel.tsv` 25 satır × 33 sütundur. Bir kesit:

```
gene                  GC3s    ENC    CAI    tAI     pI
MEFV_NM_000243.3     0.695  45.44  0.725  0.339   8.30
RPL3_NM_000967.4     0.711  44.68  0.800  0.370  10.19
RPL4_NM_000968.4     0.460  54.76  0.711  0.357  11.07
...
OR1E1_NM_003553.3    0.667  47.02  0.716  0.340   7.98
```

MEFV satırının tamamı:

```python
{'gene': 'MEFV_NM_000243.3', 'L': 781, 'GC%': 59.027, 'GC1': 0.613, 'GC2': 0.457, 'GC12': 0.535,
 'GC3': 0.7, 'GC3s': 0.695, 'ENC': 45.437, 'ENC_expected': 53.07, 'ENC_ratio': 0.144,
 'SCUO': 0.171, 'CDC': 0.117, 'RCBS': 0.167, 'Fop': 0.163, 'CAI': 0.725, 'CBI': 0.207,
 'Fop_chi2': 0.027, 'CBI_chi2': -0.145, 'RCDI': 1.28, 'CPB': 0.165, 'tAI': 0.339,
 'length': 781, 'MW': 86443.77, 'pI': 8.298, 'net_charge_7': 6.149, 'GRAVY': -0.645,
 'aromaticity': 0.059, 'aliphatic_index': 73.201, 'extinction_280': 62310,
 'IVYWREL': 0.373, 'CvP': 4.481, 'TOP_IDP': 0.175}
```

**Yorum:**

- MEFV'nin CAI'si (0,725) tek başına ayırt edici değil: yüksek ifadeli
  RPL4 (0,711) ile düşük ifadeli OR1E1 (0,716) de aynı düzeyde. χ² optimal
  kodonlarla hesaplanan CBI ise negatif (−0,145): MEFV, yüksek ifadeli
  genlerde düşüklere göre anlamlı fazla kullanılan kodonları tercih etmiyor.
- ENC (45,4), GC3s'ten beklenenin (53,1) altında. Nötrallik eğimi 0,13 ve
  R² = 0,07: GC3 değişiminin yalnızca ~%13'ü mutasyon baskısıyla açıklanıyor,
  panel genelinde seçilim baskın.
- RCDI 1,28: kodon kullanımı `host` olarak verilen setten (high + low) biraz
  uzak (1 = aynı).

### 9.3 Bir genin bütün izoformları (GENCODE) ve .crn ile saklama

Soru: MEFV'nin izoformları arasında kodon kullanımı ve protein özellikleri
nasıl değişiyor? Ayrıca indirme ve izoform ayırma işini bir kez yapıp
saklamak istiyoruz.

```python
import warnings
from pathlib import Path
from crn.fetch import Gencode
from crn.cub import CUB, read_trna_counts

warnings.filterwarnings("ignore", category=UserWarning)      # isoforms() bilgi uyarıları

# 1. GENCODE'dan geni al, izoformlara ayır, her izoformu .crn olarak sakla
gc = Gencode(verbose=False)                      # yerel v50
isoforms = gc.gene("MEFV").isoforms()            # 17 transkript -> 12 benzersiz CDS
out = Path("mefv_iso"); out.mkdir(exist_ok=True)
paths = []
for g in isoforms:
    p = out / f"{g.name}.crn"
    g.to_crn(p)
    paths.append(p)

# 2. Sonraki çalıştırmalarda yalnızca dosyalardan: ağ ve ORF araması yok
cub = CUB(paths)                                 # CUB: 12 gen, L=5828

# 3. Referans: §9.2'deki high seti (NCBI). Kaynaklar yalnızca referans için karışıyor
rows = cub.table(reference=high, trna_counts=read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa"),
                 file_name="mefv_gencode.tsv")
for r in sorted(rows, key=lambda r: -r["CAI"]):
    print(f"{r['gene']:26}{r['L']:>5}{r['GC3s']:>7.3f}{r['ENC']:>7.2f}"
          f"{r['CAI']:>7.3f}{r['tAI']:>7.3f}{r['pI']:>7.2f}")
```

```
gene                          L   GC3s    ENC    CAI    tAI     pI
MEFV_ENST00000539145.5       95  0.849  35.74  0.847  0.365   9.05
MEFV_ENST00000572244.5      101  0.816  35.63  0.831  0.356   9.63
MEFV_ENST00000536980.5      402  0.716  45.25  0.770  0.348   6.74
MEFV_ENST00000541159.5      445  0.685  46.54  0.761  0.342   6.66
MEFV_ENST00000536379.5      570  0.680  46.23  0.755  0.346   7.56
MEFV_ENST00000339854.9      601  0.669  46.62  0.754  0.346   6.80
MEFV_ENST00000538326.5      490  0.740  43.67  0.727  0.333   9.02
MEFV_ENST00000542898.5      644  0.710  45.40  0.727  0.339   6.99
MEFV_ENST00000537682.5      613  0.722  44.84  0.727  0.338   7.85
MEFV_ENST00001141017.1      781  0.695  45.44  0.725  0.339   8.30
MEFV_ENST00000956137.1      780  0.694  45.45  0.725  0.338   8.17
MEFV_ENST00000570511.5      306  0.769  39.95  0.707  0.333   9.18
```

**Yorum:** En yüksek GC3s ve CAI değerleri, en kısa iki izoformda (95 ve 101
kodon) görülüyor. Bu izoformların ENC'si de en düşük olanlar. Kısa dizilerde
bu ölçümler örneklem azlığı nedeniyle uç değerlere kayar. Bu yüzden
izoformlar arası farkları biyolojik olarak yorumlamadan önce uzunluğu hesaba
katın. Tam boy
izoform (ENST00001141017.1, 781 aa) NCBI'daki NM_000243.3 ile aynı
değerleri verir (GC3s 0,695, ENC 45,44, CAI 0,725): iki kaynak aynı CDS'te
uyumludur.

### 9.4 Kendi verinizle: FASTA'dan CUB

İndirme yapmadan, elinizdeki CDS ve protein FASTA dosyalarından analiz.
`CUB`'ın tek istediği, her biri 1 DNA + 1 RNA + 1 PRT içeren `Gene`
nesneleridir.

```python
from crn import Gene, DNA, RNA, PRT, ProteinProperties
from crn.cub import CUB

rnas = {r.name: r for r in RNA.from_fasta_file("cds.fa")}     # >g1 AUG...  >g2 AUG...
prts = {p.name: p for p in PRT.from_fasta_file("prot.fa")}    # >g1 MAAAAKK >g2 MLLLLEE

genes = []
for name, rna in rnas.items():
    g = Gene(name, DNA(rna.seq.replace("U", "T"), name=name))
    g.add_transcript(rna)
    g.add_protein(prts[name])
    genes.append(g)

cub = CUB(genes)                                  # CUB: 2 gen, L=14
for r in cub.table():
    print(r["gene"], r["L"], round(r["GC3s"], 3), r["ENC"], round(r["pI"], 2), round(r["GRAVY"], 3))
# g1 7 0.5   nan 10.0 0.186
# g2 7 0.833 nan 4.24 1.443

props = ProteinProperties(list(prts.values()))
print(props.names, round(props.molecular_weight(), 2), props.chain_count())
# ['g1', 'g2'] 1549.95 2
```

ENC'nin `nan` çıkması beklenen bir sonuçtur: 7 kodonluk genlerde her
dejenerelik sınıfı için veri yoktur (bkz. §8). Gerçek genlerde bu sorun
olmaz.

**Dosyanız DNA alfabesindeyse** `DNA.from_fasta_file(...)` ile okuyup
`dna.transcription(template=False)` ile RNA'ya çevirin.

---

### Ek: Hızlı başvuru

```python
from crn import DNA, RNA, PRT, Gene, GeneticCode, ProteinProperties
from crn.fetch import Ncbi, Gencode
from crn.cub import (CUB, Composition, Synonymous, Adaptation, EffectiveCodons,
                     Divergence, CodonContext, read_trna_counts)

iso = Ncbi().gene("MEFV").isoforms()          # veri -> izoformlar
cub = CUB(iso)                                # sayımlar bir kez
Composition(cub).gc3s()                       # B
Synonymous(cub, reference=high).rscu()        # C
Adaptation(cub, reference=high).cai()         # D
ProteinProperties(cub).summary()              # E
EffectiveCodons(cub).enc()                    # F
Divergence(cub, reference=genome).rcdi()      # G
CodonContext(cub, reference=genome).cpb()     # H
cub.table(reference=high, file_name="t.tsv")  # hepsi bir arada
```
