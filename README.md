# crn Çerçeve Rehberi

Bu rehber `crn` çerçevesini baştan sona anlatır:

1. **crn ne işe yarar?** Çerçevenin amacı, katmanları ve tasarım ilkeleri.
2. **`crn.io`**: diziler, genler ve dosyalar. Her modül ve içindekiler.
3. **`crn.fetch`**: NCBI ve GENCODE'dan veri çekme.
4. **`crn.cub`**: kodon kullanım yanlılığı ve protein özellikleri (A–H).
5. **Ayrıntılı kullanım**: çerçevenin en etkili biçimde nasıl kullanılacağı,
   uçtan uca senaryolarla.
6. **İyi uygulamalar, sık hatalar ve sınırlılıklar.**

Örneklerdeki çıktılar gerçek çalıştırmalardan alınmıştır (insan, NCBI
RefSeq ve GENCODE v50). Her metodun tüm parametreleri ve hata durumları için
başvuru belgeleri:

- [crn_kullanim_kilavuzu.md](docs/crn_kullanim_kilavuzu.md): `io` ve `fetch`
  sınıflarının metot metot başvurusu,
- [crn_fetch_ncbi.md](docs/crn_fetch_ncbi.md): `Ncbi` ve `NcbiSession` için
  ayrıntılı kılavuz,
- [crn_io_features.md](docs/crn_io_features.md), [crn_io_type.md](docs/crn_io_type.md).

---

## İçindekiler

- [1. crn nedir, ne işe yarar?](#1-crn-nedir-ne-işe-yarar)
  - [1.1. Amaç](#11-amaç)
  - [1.2. Üç katman](#12-üç-katman)
  - [1.3. Veri akışı](#13-veri-akışı)
  - [1.4. Tasarım ilkeleri](#14-tasarım-ilkeleri)
  - [1.5. Kurulum ve içe aktarma](#15-kurulum-ve-içe-aktarma)
- [2. crn.io — diziler, genler, dosyalar](#2-crnio--diziler-genler-dosyalar)
  - [2.1. crn_io_type](#21-crn_io_type--molekül-ve-zincir-tipleri)
  - [2.2. crn_io_const](#22-crn_io_const--sabitler-ve-genetik-kodlar)
  - [2.3. crn_io_features](#23-crn_io_features--ek-bilgi-deposu)
  - [2.4. crn_io_seqs](#24-crn_io_seqs--sequence-fasta-ve-crn-biçimi)
  - [2.5. crn_io_dna](#25-crn_io_dna--dna)
  - [2.6. crn_io_rna](#26-crn_io_rna--rna)
  - [2.7. crn_io_prt](#27-crn_io_prt--protein)
  - [2.8. crn_io_gene](#28-crn_io_gene--gene)
- [3. crn.fetch — veri kaynakları](#3-crnfetch--veri-kaynakları)
  - [3.1. crn_fetch_session](#31-crn_fetch_session--ncbisession)
  - [3.2. crn_fetch_ncbi](#32-crn_fetch_ncbi--ncbi)
  - [3.3. crn_fetch_gencode](#33-crn_fetch_gencode--gencode)
  - [3.4. Ncbi mi, Gencode mu?](#34-ncbi-mi-gencode-mu)
- [4. crn.cub — kodon kullanım yanlılığı](#4-crncub--kodon-kullanım-yanlılığı)
  - [4.1. Ortak yapı](#41-ortak-yapı)
  - [4.2. A. Ham sayımlar — CUB](#42-a-ham-sayımlar--cub)
  - [4.3. B. Nükleotit bileşimi — Composition](#43-b-nükleotit-bileşimi--composition)
  - [4.4. C. Sinonim kodon tercihi — Synonymous](#44-c-sinonim-kodon-tercihi--synonymous)
  - [4.5. D. Referans setine dayalı indeksler — Adaptation](#45-d-referans-setine-dayalı-indeksler--adaptation)
  - [4.6. E. Protein özellikleri (E-A … E-G)](#46-e-protein-özellikleri-e-a--e-g)
  - [4.7. F. Etkin kodon sayısı ve grafikler — EffectiveCodons](#47-f-etkin-kodon-sayısı-ve-grafikler--effectivecodons)
  - [4.8. G. Entropi ve sapma — Divergence](#48-g-entropi-ve-sapma--divergence)
  - [4.9. H. tRNA ve kodon bağlamı — CodonContext](#49-h-trna-ve-kodon-bağlamı--codoncontext)
- [5. Ayrıntılı kullanım: en etkili biçim](#5-ayrıntılı-kullanım-en-etkili-biçim)
- [6. İyi uygulamalar ve sık hatalar](#6-iyi-uygulamalar-ve-sık-hatalar)
- [7. Sınırlılıklar](#7-sınırlılıklar)

---

## 1. crn nedir, ne işe yarar?

### 1.1. Amaç

`crn`, moleküler biyolojinin santral dogmasını (DNA → RNA → protein) Python
nesneleriyle modelleyen bir biyoinformatik çerçevesidir. Üç iş için
tasarlanmıştır:

1. **Dizileri doğru modellemek.** DNA, RNA ve protein ayrı sınıflardır. Her
   biri kendi alfabesini denetler, yazılış yönünü (5'→3' / 3'→5') bilir ve
   ek bilgilerini (`features`) taşır. Dönüşümler biyolojik süreçlerin
   adlarıyla yapılır: `transcription`, `reverse_transcription`,
   `translation`, `find_orfs`.
2. **Gerçek veriyi tek komutla almak.** NCBI RefSeq ve GENCODE'dan bir genin
   genomik DNA'sı, tüm transkriptleri ve proteinleri tek bir `Gene`
   nesnesinde toplanır. İndirilen her şey yerel `sys/` klasöründe saklanır;
   ikinci istekte ağa çıkılmaz.
3. **Kodon kullanımını ve proteinleri ölçmek.** `crn.cub`, kodon kullanım
   yanlılığının (codon usage bias, CUB) literatürdeki 34 temel ölçümünü
   (A–H bölümleri) ve 34 protein özelliğini (E-A … E-G) hesaplar. Ölçümler
   bilinen değerlerle doğrulanmıştır (ExPASy ProtParam, UniProt, elle
   hesaplama).

Tipik bir soru: *"MEFV geninin kodon kullanımı ribozomal genlere mi, yoksa
koku reseptörü genlerine mi benziyor; yanlılığı mutasyon baskısıyla mı,
seçilimle mi açıklanır?"* `crn` ile bu soru, verinin indirilmesinden
grafiklerin çizilmesine kadar aynı çerçevede yanıtlanır.

### 1.2. Üç katman

```
crn/
├── io/       Veri modeli: dizi sınıfları, Gene kapsayıcısı, sabitler, dosya biçimleri
├── fetch/    Veri kaynakları: NCBI (E-utilities) ve GENCODE (FTP), yerel sys/ önbelleği
└── cub/      Analiz: kodon kullanım yanlılığı (A–H) ve protein özellikleri (E-A … E-G)
```

Katmanlar tek yönlü bağımlıdır: `fetch` ve `cub`, `io`'yu kullanır; `io`
hiçbirine bağlı değildir. `cub` ise `fetch`'e bağlı değildir. Ona hangi
kaynaktan gelirse gelsin bir `Gene` listesi vermek yeterlidir: NCBI,
GENCODE, FASTA'dan elle kurulmuş bir gen ya da `.crn` dosyası olabilir.

| Katman | Modüller | Başlıca sınıflar |
|---|---|---|
| `crn.io` | `crn_io_type`, `crn_io_const`, `crn_io_features`, `crn_io_seqs`, `crn_io_dna`, `crn_io_rna`, `crn_io_prt`, `crn_io_gene` | `Molecule`, `Strand`, `GenomicStrand`, `Maps`, `GeneticCode`, `Features`, `Sequence`, `DNA`, `RNA`, `PRT`, `Gene` |
| `crn.fetch` | `crn_fetch_session`, `crn_fetch_ncbi`, `crn_fetch_gencode` | `NcbiSession`, `Ncbi`, `Gencode` |
| `crn.cub` | `crn_cub_counts` (A), `crn_cub_composition` (B), `crn_cub_synonymous` (C), `crn_cub_reference` (D), `crn_cub_protein_*` (E-A … E-G), `crn_cub_enc` (F), `crn_cub_divergence` (G), `crn_cub_context` (H) | `CUB`, `Composition`, `Synonymous`, `Adaptation`, `ProteinMass`, `ProteinAtoms`, `ProteinCharge`, `ProteinHydropathy`, `ProteinOptics`, `ProteinStability`, `ProteinDisorder`, `EffectiveCodons`, `Divergence`, `CodonContext` |

### 1.3. Veri akışı

```
   NCBI RefSeq ─┐                         ┌─ FASTA / .crn dosyası
   GENCODE ─────┤  crn.fetch              │  (Sequence.from_fasta_file,
                │  Ncbi.gene("MEFV")      │   Gene.from_crn)
                ▼  Gencode.gene("MEFV")   ▼
          ┌─────────────────────────────────────────┐
          │ Gene                                    │  crn.io
          │  ├── dna          genomik DNA (1)       │
          │  ├── transcripts  {id: RNA} (n)         │
          │  ├── proteins     {id: PRT} (m)         │
          │  └── features     gene_id, aliases ...  │
          └─────────────────────────────────────────┘
                │  gene.match_transcripts_to_proteins()
                ▼
          [Gene, Gene, ...]   her biri: 1 DNA + 1 transkript + 1 protein
                │  CUB(genes=[...])   transkript → CDS (ORF ile)
                ▼
          CUB ──► Composition, Synonymous, Adaptation, EffectiveCodons,
                  Divergence, CodonContext, Protein*          crn.cub
                │
                ▼
          sayılar, sözlükler, matplotlib grafikleri
```

### 1.4. Tasarım ilkeleri

Çerçevenin her yerinde aynı kurallar geçerlidir; bunları bilmek kullanımı
kolaylaştırır.

1. **Ek bilgiler `features` içindedir.** Bir dizinin ya da genin kimliği,
   kaynağı, CDS koordinatı gibi bilgiler sınıf alanı değil, anahtar → değer
   listesi deposudur (`Features`). Yeni bir bilgi taşımak için sınıfı
   değiştirmek gerekmez.
2. **Root alanları istisnadır.** Yalnızca dizinin biyolojik anlamını
   değiştiren alanlar sınıfın kendisindedir: `strand` (yazılış yönü),
   RNA'da `genomic_strand`, tümleme tablosu `maps`.
3. **Diziler büyük harfe çevrilir ve doğrulanır.** `DNA("atgx")` hata verir:
   `ValueError: DNA dizisinde geçersiz karakter(ler): X`.
4. **Kopyalar bağımsızdır.** `seq[:]`, `complementary()`, `transcription()`
   gibi işlemler yeni nesne döndürür ve `features`'ın kopyasını taşır.
5. **Hatalar sessizce yutulmaz.** Belirsiz durumlar ya `ValueError` /
   `KeyError` ile durur, ya da `UserWarning` / `CubWarning` ile bildirilir.
6. **Sabitler tek yerdedir.** Alfabeler, kodon listesi, genetik kodlar,
   kütleler, pKa kümeleri, hidropati ölçekleri, tRNA wobble cezaları
   `crn/io/crn_io_const.py` içindedir. Fonksiyonlara parametre olarak
   girer ve değiştirilebilir.
7. **Analiz sınıfları birbirini kullanır.** `cub` içindeki her bölüm önceki
   bölümlerin metotlarını yeniden kullanır. Örneğin ENC oranı, B'deki GC3s'i
   ve F'deki ENC'yi birleştirir. Aynı formül iki yerde yazılmamıştır.

### 1.5. Kurulum ve içe aktarma

Gereksinimler: Python ≥ 3.12; bağımlılıklar `tqdm`, `requests`,
`matplotlib` (`pyproject.toml` içinde, otomatik kurulur).

```bash
python -m venv .venv
```

```bash
source .venv/bin/activate.fish
```

```bash
pip install -e .
```

Betikleri her zaman sanal ortamın Python'u ile çalıştırın; sistem Python'unda
bağımlılıklar yoktur:

```bash
.venv/bin/python test/cub_test.py
```

İçe aktarma:

```python
# Veri modeli (crn ya da crn.io üzerinden)
from crn import DNA, RNA, PRT, Gene, Features, Strand, GenomicStrand, GeneticCode, Maps

# Veri kaynakları
from crn.fetch import Ncbi, Gencode, NcbiSession

# Analiz
from crn.cub import (CUB, CubWarning, Composition, Synonymous, Adaptation,
                     EffectiveCodons, Divergence, CodonContext, read_trna_counts,
                     ProteinMass, ProteinAtoms, ProteinCharge, ProteinHydropathy,
                     ProteinOptics, ProteinStability, ProteinDisorder)

# Sabitler
from crn.io.crn_io_const import CODONS, AMINO_ACIDS, RNA_BASES, HYDROPATHY_SCALES
```

---

## 2. crn.io — diziler, genler, dosyalar

`crn.io`, çerçevenin veri modelidir. Diğer iki katman yalnızca bu modülün
nesnelerini üretir (`fetch`) ya da tüketir (`cub`).

```
Sequence (crn_io_seqs)          ortak davranış: doğrulama, dilimleme, FASTA, CRN
 ├── DNA  (crn_io_dna)          tümleme, transkripsiyon, ORF arama
 ├── RNA  (crn_io_rna)          ters transkripsiyon, translasyon, ORF eşleme
 └── PRT  (crn_io_prt)          protein
Gene     (crn_io_gene)          1 DNA + n RNA + m PRT + features; izoform ayırma
Features (crn_io_features)      anahtar → değer listesi deposu
crn_io_type / crn_io_const      enum'lar, alfabeler, genetik kodlar, tablolar
```

### 2.1. crn_io_type — molekül ve zincir tipleri

Üç `Enum` tanımlar:

| Enum | Üyeler | Ne işe yarar |
|---|---|---|
| `Molecule` | `DNA`, `RNA`, `PRT`, `SEQ` | Dizinin molekül tipi. CRN dosyasında hangi sınıfın okunacağını belirler. |
| `Strand` | `FivePrime` (5'→3'), `ThreePrime` (3'→5') | Dizinin **yazılış yönü**. Tümleme ve transkripsiyon yönü çevirir; translasyon 3'→5' yazılmış RNA'yı önce ters çevirir. |
| `GenomicStrand` | `Plus` (+), `Minus` (−), `Unstranded` (.), `Unknown` (?) | Yalnızca RNA'da: transkriptin genomun hangi zincirinden geldiği. |

`Strand` ile `GenomicStrand` farklı şeylerdir. MEFV genomun (−) zincirindedir
(`GenomicStrand.Minus`), ama mRNA'sı yine 5'→3' yazılır
(`Strand.FivePrime`).

### 2.2. crn_io_const — sabitler ve genetik kodlar

Çerçevenin bütün sabitleri buradadır. Bir sabit eksikse buraya eklenir ve
fonksiyonlara parametre olarak verilir.

**Alfabeler ve kodonlar**

| Sabit | Değer / anlam |
|---|---|
| `DNA_ALPHABET`, `RNA_ALPHABET`, `PRT_ALPHABET` | IUPAC kodları dahil izin verilen harfler (boşluk `-`, protein için `*` stop) |
| `DNA_BASES`, `RNA_BASES` | `"TCAG"`, `"UCAG"` (NCBI tablo sırası) |
| `CODONS` | 64 RNA kodonu, `UUU, UUC, …, GGG` |
| `AMINO_ACIDS` | 20 standart amino asit, `"ACDEFGHIKLMNPQRSTVWY"` |
| `STOP_SYMBOL` | `"*"` |

**`GeneticCode`**: `STANDARD` (NCBI 1) ve `VERTEBRATE_MITOCHONDRIAL`
(NCBI 2). Her üyenin `table` (kodon → amino asit), `starts` ve `stops`
özellikleri vardır. ORF arama, translasyon ve tüm CUB ölçümleri bu
nesneyle çalışır:

```python
from crn import RNA, GeneticCode
RNA("AUGUGAAAA").translation()                                   # M*K
RNA("AUGUGAAAA").translation(GeneticCode.VERTEBRATE_MITOCHONDRIAL)  # MWK  (UGA = Trp)
```

**`Maps`**: `str.translate` tabloları: DNA/RNA tümleme (IUPAC dahil),
transkripsiyon (kalıp ya da kodlayan zincirden), ters transkripsiyon,
`CODING_TO_RNA_MAP` (T→U), `RNA_TO_CODING_MAP` (U→T).

**Protein tabloları** (E bölümü kullanır): `AMINO_ACID_GROUPS`,
`RESIDUE_AVERAGE_MASS`, `RESIDUE_MONOISOTOPIC_MASS`, `WATER_*`,
`ELEMENTS`, `AMINO_ACID_FORMULAS`, `WATER_FORMULA`, `POSITIVE_RESIDUES`,
`NEGATIVE_RESIDUES`, `IONIZABLE_POSITIVE`, `IONIZABLE_NEGATIVE`,
`PKA_SETS` (EMBOSS, Lehninger, Bjellqvist), `HYDROPATHY_SCALES`
(Kyte–Doolittle, Hopp–Woods, Eisenberg, Wimley–White, Janin),
`BOMAN_SCALE`, `EXTINCTION_COEFFICIENTS_280`, `RESIDUE_SPECIFIC_VOLUMES`,
`RESIDUE_VOLUMES`, `RESIDUE_MAX_ASA`, `IVYWREL_RESIDUES`, `CVP_CHARGED`,
`CVP_POLAR`, `TURN_FORMING_RESIDUES`, `TOP_IDP_SCALE`.

**tRNA ve izoform sabitleri**: `TAI_S_VALUES` (dos Reis wobble cezaları),
`TAI_WOBBLE_PAIRS` (kodonun 3. bazını okuyan antikodonlar),
`TRANSCRIPT_PRIORITY` (`("NM_", "XM_")`: aynı CDS'li izoformlardan hangisi
tutulur).

### 2.3. crn_io_features — ek bilgi deposu

`Features`, her anahtara bir **değer listesi** bağlayan küçük bir depodur.
Her `Sequence` ve her `Gene` bir `features` taşır.

| Metot | Anlam |
|---|---|
| `add(key, value)` | Listeye ekler (çok değerli bilgi, ör. `aliases`) |
| `set(key, value)` | Listeyi tek değerle değiştirir |
| `get(key, default=None)` | Değer listesi |
| `first(key, default=None)` | İlk değer (tek değerli bilgiler için en sık kullanılan) |
| `keys()`, `items()`, `to_dict()` | Dolaşma ve dışa aktarma |
| `copy()`, `remove(key)`, `clear()` | Bağımsız kopya, silme |
| `in`, `len`, `for` | Varlık, sayı, anahtarlar |

```python
f = Features()
f.add("alias", "FMF"); f.add("alias", "MEF"); f.set("chromosome", "16")
f.get("alias")         # ['FMF', 'MEF']
f.first("alias")       # 'FMF'
f.to_dict()            # {'alias': ['FMF', 'MEF'], 'chromosome': ['16']}
```

Çerçevenin kendisinin yazdığı başlıca anahtarlar:

| Anahtar | Nerede | Anlam |
|---|---|---|
| `source`, `gene_id`, `gene_name`, `description`, `aliases`, `chromosome`, `genomic_region` | Gene (fetch) | Gen kimliği ve konumu |
| `transcript_id`, `protein_id`, `predicted`, `title` | RNA / PRT (fetch) | Kayıt kimlikleri |
| `utr5`, `cds`, `utr3` | RNA (fetch) | Transkript üzerinde 1 tabanlı `(başlangıç, bitiş)` |
| `orf_start`, `orf_end`, `orf_strand`, `orf_frame` | `find_orfs`, `match_orf` | ORF'nin konumu |
| `isoform_of`, `transcript_id`, `protein_id` | `match_transcripts_to_proteins` | İzoformun kaynağı |

### 2.4. crn_io_seqs — Sequence, FASTA ve CRN biçimi

`Sequence`, `DNA`, `RNA` ve `PRT`'nin ortak temel sınıfıdır:

- **Özellikler:** `seq`, `name`, `features`.
- **Dizi gibi davranır:** `len(s)`, `s[i]` (harf), `s[a:b]` (aynı tipte
  yeni nesne, features kopyasıyla), `in` (alt dizi), `==` (dizi ya da metin
  ile), `for`, `hash`.
- **FASTA:** `cls.from_fasta_file(yol)` her zaman liste döndürür (`.gz`
  desteklenir, ad başlığın ilk kelimesidir). `s.to_fasta(yol)` diziyi
  yazar; features FASTA'ya yazılmaz.
- **CRN:** `s.to_crn(yol)` / `cls.from_crn(yol)`. CRN, çerçevenin kendi
  ikili biçimidir. Dizi, ad, root alanları (strand …) ve **tüm features**
  (int, float, tuple, enum dahil) yazılır, sonda CRC32 denetimi vardır.
  pickle kullanılmaz.

```python
dna = DNA("atgGCTtgcAAataa", name="ornek")
dna.features.set("kaynak", "elle")
dna.to_fasta("ornek.fa")   # >ornek this sequence contains 1 feature but not listed
dna.to_crn("ornek.crn")
back = DNA.from_crn("ornek.crn")
back.features.first("kaynak"), back.strand       # ('elle', Strand.FivePrime)
DNA.from_fasta_file("ornek.fa")                  # [ornek: ATGGCTTGCAAATAA <15> [DNA] 5'->3']
```

> **Kural:** Analize girecek veriyi FASTA ile değil, CRN ile saklayın. FASTA
> CDS koordinatlarını, kimlikleri ve yönü kaybeder.

### 2.5. crn_io_dna — DNA

`DNA(seq, name="", strand=Strand.FivePrime, maps=Maps.DEFAULT_COMPLEMENT_MAP_FOR_DNA)`

| Metot | Ne yapar |
|---|---|
| `complement()` | **Yerinde** tümler, yönü çevirir, `self` döner |
| `complementary()` | Tümleyiciyi yeni nesne olarak döner |
| `transcription(template=True, genomic_strand=…, maps=None)` | `template=True`: bu DNA kalıp zincirdir, tümleyici RNA üretilir. `template=False`: kodlayan zincirdir, yalnızca T→U yapılır. |
| `find_orfs(minimum_length_of_orf=-1, genetic_code=…, alternative_starts=False)` | İki zincirde, altı çerçevede tüm ORF'leri bulur |
| `from_sequence(seq)` | Başka bir `Sequence`'ten DNA kurar |

```python
dna = DNA("atgGCTtgcAAataa", name="ornek")
dna                         # ornek: ATGGCTTGCAAATAA <15> [DNA] 5'->3'
dna.complementary()         # ornek: TACCGAACGTTTATT <15> [DNA] 3'->5'
dna[0:6]                    # ornek: ATGGCT <6> [DNA] 5'->3'
mrna = dna.transcription(template=False, genomic_strand=GenomicStrand.Plus)
mrna                        # ornek: AUGGCUUGCAAAUAA <15> [RNA] 5'->3' (+)
dna.transcription()         # ornek: UACCGAACGUUUAUU <15> [RNA] 3'->5' (?)

g = DNA("CCATGAAATTTGGGTAACCATGCCCTAGTT", name="g")
for orf in g.find_orfs():
    print(orf.name, orf.seq, orf.features.first("orf_start"), orf.features.first("orf_end"))
# g_ORF1 ATGAAATTTGGGTAA 3 17
# g_ORF2 ATGCCCTAG 20 28
len(g.find_orfs(minimum_length_of_orf=12))   # 1
```

NCBI ve GENCODE'dan gelen gen dizileri **kodlayan zincirdir**; bu yüzden
gen DNA'sından mRNA'ya geçerken `transcription(template=False)` kullanılır.

### 2.6. crn_io_rna — RNA

`RNA(seq, name="", strand=Strand.FivePrime, genomic_strand=GenomicStrand.Unknown, maps=…)`

| Metot | Ne yapar |
|---|---|
| `complement()`, `complementary()` | DNA'daki gibi; `genomic_strand` değişmez |
| `reverse_transcription(maps=…)` | Tümleyici cDNA (ilk zincir) |
| `translation(genetic_code=…)` | İlk nükleotitten üçerli okur; stop `*`, belirsiz kodon `X` |
| `match_orf(protein, genetic_code=…)` | Çevirisi proteinle birebir aynı olan ORF'yi RNA olarak döner (yoksa `None`) |

`match_orf`, CDS koordinatı bilinmeyen ya da güvenilmeyen bir transkriptte
proteini kodlayan bölgeyi bulur. CUB hesaplarının temelidir:

```python
r = RNA("GGCAUGAAAUUUGGGUAAGC", name="t1")
m = r.match_orf("MKFG")
m                                    # t1_ORF1: AUGAAAUUUGGGUAA <15> [RNA] 5'->3' (?)
m.features.first("cds")              # (1, 15)   artık RNA'nın tamamı CDS
m.features.first("orf_start"), m.features.first("orf_end")   # (4, 18)  eski RNA üzerindeki konum
r.match_orf("MKW")                   # None
```

### 2.7. crn_io_prt — protein

`PRT(seq, name="")`: protein dizisi; alfabe 20 amino asit, IUPAC belirsizlik
kodları, U (selenosistein), O (pirolizin), `*` ve `-`'dir. Yönü N→C'dir.
Protein hesapları (kütle, pI, GRAVY …) `crn.cub`'daki E bölümündedir.

### 2.8. crn_io_gene — Gene

`Gene(name, dna=None)` bir gene ait her şeyi bir arada tutar:

```
Gene (ör. MEFV)
 ├── dna          : DNA            genomik dizi (tek)
 ├── transcripts  : {id: RNA}      transkript varyantları
 ├── proteins     : {id: PRT}      proteinler
 └── features     : Features       gene_id, aliases, genomic_region ...
```

| Metot | Ne yapar |
|---|---|
| `add_transcript(rna, id=None)`, `transcript(id)`, `remove_transcript(id)`, `transcripts` | Transkript yönetimi (kimlik verilmezse RNA'nın adı) |
| `add_protein(prt, id=None)`, `protein(id)`, `remove_protein(id)`, `proteins` | Protein yönetimi |
| `match_transcripts_to_proteins(genetic_code=…, unique_cds=True, priority=TRANSCRIPT_PRIORITY)` | Geni izoformlarına ayırır (aşağıda) |
| `to_crn(yol)`, `Gene.from_crn(yol)` | Genin tamamını tek dosyada saklar |

**`match_transcripts_to_proteins`.** Bir gende transkriptlerle proteinler
arasındaki eşleşme kendiliğinden bilinmez. Bu metot her transkript için
`RNA.match_orf` ile proteinleri dener. Eşleşen her (transkript, protein)
çifti, genomik DNA'nın kopyasıyla birlikte yeni bir `Gene` olur. Böylece
her izoform tam olarak **bir DNA, bir transkript ve bir protein** taşır;
`CUB` bu yapıyı ister.

- Yeni genlerin adı `<gen>_<transkript>` olur (ör. `MEFV_NM_000243.3`).
  Features'a `isoform_of`, `transcript_id` ve `protein_id` yazılır.
- Hiçbir proteinle eşleşmeyen transkriptler (kodlamayan RNA'lar) ve eşsiz
  kalan proteinler `UserWarning` ile bildirilir.
- `unique_cds=True` (varsayılan): aynı CDS'i veren izoformlardan yalnızca
  biri tutulur. Önce `NM_`, sonra `XM_` tercih edilir; eşitlikte ilk
  eklenen kalır. Çıkarılanlar `UserWarning` ile bildirilir. Böylece yalnızca
  UTR'leri farklı olan transkriptler çok genli ölçümleri çarpıtmaz.

```python
gene = Gene("ORNEK", DNA("CCATGAAATTTGGGTAACC"))
gene.add_transcript(RNA("CCAUGAAAUUUGGGUAACC", name="T1"))
gene.add_transcript(RNA("AUGAAAUUUGGGUAA"), "T2")     # aynı CDS, farklı UTR
gene.add_protein(PRT("MKFG", name="P1"))
gene                                   # ORNEK: DNA <19>, 2 transkript, 1 protein
iso = gene.match_transcripts_to_proteins()
# UserWarning: ORNEK: aynı CDS'i veren izoformlar çıkarıldı (T2 -> T1)
[(i.name, list(i.transcripts), list(i.proteins)) for i in iso]
# [('ORNEK_T1', ['T1'], ['P1'])]
```

---

## 3. crn.fetch — veri kaynakları

`crn.fetch`, gerçek genleri `crn.io` nesnelerine dönüştürür. İki kaynak
vardır ve ikisi de aynı arayüzü sunar:

| Metot | Döndürdüğü |
|---|---|
| `gene(gen)` | `Gene`: genomik DNA + tüm transkriptler + proteinler |
| `dna(gen)` | Genin genomik DNA'sı (kodlayan zincir, 5'→3') |
| `rna(transkript)` | Transkript; features'ta `utr5`, `cds`, `utr3` |
| `protein(kimlik)` | Protein; transkript kimliği de verilebilir |
| `region(kromozom, başlangıç, bitiş, zincir)` | Herhangi bir genomik bölge |

Her ikisi de veriyi **yerel `sys/` klasörüne** yazar ve sonraki isteklerde
oradan okur. Bir kez indirilen gen için ağa bir daha çıkılmaz.

### 3.1. crn_fetch_session — NcbiSession

NCBI E-utilities için HTTP oturumudur. `Ncbi` tarafından kullanılır; doğrudan
nadiren gerekir.

- **Hız sınırı:** anahtarsız saniyede 3, API anahtarıyla 10 istek.
  İş parçacığı güvenlidir.
- **Yeniden deneme:** 429 ve 5xx hatalarında geri çekilmeli otomatik deneme
  (`max_retries=5`), bağlantı havuzu.
- **Parametreler:** `NcbiSession(api_key=None, e_mail=None, timeout=30, max_retries=5, ...)`.
  Context manager olarak da kullanılabilir.

```python
import os
from crn.fetch import Ncbi, NcbiSession

session = NcbiSession(api_key=os.environ.get("NCBI_API_KEY"), e_mail="sen@ornek.com")
ncbi = Ncbi(organism="human", session=session)
```

> API anahtarını koda yazmayın; ortam değişkeninden (`NCBI_API_KEY`) okuyun.

### 3.2. crn_fetch_ncbi — Ncbi

`Ncbi(organism="human", session=None, sys_dir=None, verbose=True)`.
Organizmalar: `human`, `mouse`, `rat`, `zebrafish`, `fruit_fly`, `worm`.

Önbellek yapısı (`sys/ncbi/`): `symbols/` (sembol → Gene ID),
`genes/<id>.json` (genin özeti), `nuccore/` (transkript FASTA ve özellik
tablosu), `protein/`, `region/`.

```python
ncbi = Ncbi(organism="human", verbose=False)
ncbi.search_gene("MEFV")          # '4210'
mefv = ncbi.gene("MEFV")          # ya da ncbi.gene("4210")
mefv                              # MEFV: DNA <14607>, 2 transkript, 2 protein
mefv.features.to_dict()
# {'source': ['NCBI'], 'gene_id': ['4210'], 'gene_name': ['MEFV'],
#  'description': ['MEFV innate immunity regulator, pyrin'], 'organism': ['Homo sapiens'],
#  'chromosome': ['16'], 'map_location': ['16p13.3'],
#  'aliases': ['FMF', 'MEF', 'PAAND', 'TRIM20'],
#  'genomic_region': [('NC_000016.10', 3242027, 3256633, '-')]}

for tid, rna in mefv.transcripts.items():
    print(tid, len(rna), rna.genomic_strand, rna.features.first("cds"), rna.features.first("protein_id"))
# NM_000243.3    3506 GenomicStrand.Minus (47, 2392) NP_000234.1
# NM_001198536.2 3048 GenomicStrand.Minus (47, 1384) NP_001185465.2

rna = ncbi.rna("NM_000243")                    # sürüm numarası olmadan da olur
start, end = rna.features.first("cds")
cds = rna[start - 1:end]
cds.translation().seq.rstrip("*") == ncbi.protein("NM_000243").seq   # True

ncbi.region("NC_000016.10", 3242028, 3242060, "-", name="bolge")
# bolge: CTGGGCCTCATGAATAAAAAAATTAGAAAATCT <33> [DNA] 5'->3'
ncbi.update("MEFV", force=True)               # önbelleği yenile
```

Ayrıntılı kullanım için: [crn_fetch_ncbi.md](docs/crn_fetch_ncbi.md).

### 3.3. crn_fetch_gencode — Gencode

`Gencode(species="human", release=None, sys_dir=None, verbose=True)`.
Türler: `human`, `mouse`. İlk kullanımda GENCODE sürümünün genomu,
transkriptleri, proteinleri ve GTF'si indirilip dizinlenir. İnsan için bu
yaklaşık 6 GB'tır ve bir kez yapılır. Sonraki aramalar dizin üzerinden
anında yapılır.

| Metot | Ne yapar |
|---|---|
| `update(force=False)`, `is_ready()` | İndirme/hazırlama; hazır mı |
| `installed_releases()`, `latest_release()`, `release` | Sürüm yönetimi |
| `gene`, `dna`, `rna`, `protein`, `region` | Yukarıdaki ortak arayüz |

```python
gc = Gencode(verbose=False)
gc.release, gc.is_ready()          # ('v50', True)
g = gc.gene("MEFV")                # MEFV: DNA <16563>, 17 transkript, 15 protein
g.features.to_dict()
# {'source': ['GENCODE'], 'release': ['v50'], 'gene_id': ['ENSG00000103313.16'],
#  'gene_name': ['MEFV'], 'gene_type': ['protein_coding'], 'chromosome': ['chr16'],
#  'start': [3240392], 'end': [3256954], 'chromosome_strand': ['-']}
gc.rna("MEFV-201")                 # transkript adı, ENST ya da sürümsüz ENST
# ENST00000219596.6: CUACCAGAAG...  <3506> [RNA] 5'->3' (-)
gc.protein("MEFV-201")             # ENSP00000219596.1 ... <781> [PRT] N->C
gc.region("chr16", 3242028, 3242060, "-")
```

### 3.4. Ncbi mi, Gencode mu?

| | `Ncbi` | `Gencode` |
|---|---|---|
| Kapsam | RefSeq (NM_/XM_/NR_/XR_), 6 organizma | GENCODE/Ensembl, insan ve fare |
| İlk kullanım | Gen başına birkaç istek; KB'lar | Tüm sürüm bir kez indirilir (~6 GB) |
| Transkript sayısı | Az, kürasyonlu (MEFV: 2) | Çok, deneysel kanıtlı ve öngörülen (MEFV: 17) |
| En uygun | Az sayıda gen, kürasyonlu izoformlar | Genom ölçeğinde analiz, gen tipine göre tarama |

İki kaynağın koordinatları ve transkript kümeleri farklıdır; bir analizde
tek kaynak kullanın.

---

## 4. crn.cub — kodon kullanım yanlılığı

### 4.1. Ortak yapı

`crn.cub` bir çekirdek ve onun üzerine kurulan ölçüm sınıflarından oluşur:

```
CUB (A, crn_cub_counts)                 genleri hazırlar, kodonlara böler, ham sayımlar
 ├── Composition      (B)               nükleotit bileşimi
 ├── Synonymous       (C)               RSCU, optimal kodon, Fop          [reference, background]
 ├── Adaptation       (D)               w, CAI, CBI                       [reference, background]
 ├── ProteinMass … ProteinDisorder (E)  protein özellikleri
 ├── EffectiveCodons  (F)               ENC, nötrallik, PR2, grafikler
 ├── Divergence       (G)               SCUO, CDC, RCDI, RCBS             [reference]
 └── CodonContext     (H)               tAI, dinükleotit ρ, CPS/CPB       [reference]
```

Bütün sınıflar aynı kurallarla çalışır:

1. **Girdi her zaman bir `CUB` nesnesidir.** `Composition(cub)`,
   `Adaptation(cub, reference=high)` …
2. **`gene` parametresi:** verilmezse tüm genler birlikte, verilirse
   yalnızca o gen hesaplanır. Gen adı `CUB` içindeki addır (ör.
   `"MEFV_NM_000243.3"`).
3. **`codons`, `amino_acids`, ölçekler, pKa kümeleri** parametredir ve
   varsayılanları `crn_io_const`'tan gelir.
4. **Referans setleri** ikinci bir `CUB`'dır. `reference` genellikle yüksek
   ifade edilen genler (ribozomal proteinler), `background` düşük ifade
   edilen genlerdir. Verilmezse `cub`'ın kendisi kullanılır.
5. **Genetik kod** `CUB(genes, genetic_code=…)` ile bir kez seçilir; diğer
   sınıflar ondan alır. Referansın kodu farklıysa `ValueError` verilir.

### 4.2. A. Ham sayımlar — CUB

`CUB(genes: List[Gene | str | Path], genetic_code=GeneticCode.STANDARD)`

**Gen kuralı:** listedeki her gen **tam olarak bir genomik DNA, bir
transkript ve bir protein** içermelidir. Uymayan gen varsa hesap başlamadan
durur:

```python
CUB(genes=[ncbi.gene("MEFV")])
# ValueError: Her gende bir genomik DNA, bir transkript ve bir protein olmalı
#   (çok izoformlu genler için Gene.match_transcripts_to_proteins); hesap durduruldu:
#   MEFV: 2 transkript (NM_000243.3, NM_001198536.2); 2 protein ...
CUB(genes=iso[0])
# TypeError: genes bir liste olmalı, Gene verildi; tek gen için genes=[gene] yazın
```

**CDS'e indirgeme:** her genin transkriptinde proteinle birebir örtüşen ORF
bulunur (`RNA.match_orf`) ve transkriptin yerine konur. Gen **yerinde**
düzenlenir ve `CubWarning` ile bildirilir. Örtüşen ORF yoksa `ValueError`
verilir. Kodlamayan RNA'lar için kodon kullanımı tanımsız olduğundan hesaba
girmez.

```python
iso = ncbi.gene("MEFV").match_transcripts_to_proteins()
cub = CUB(genes=iso)
# CubWarning: MEFV_NM_000243.3: gen yeniden düzenlendi (NM_000243.3, NP_000234.1
#   proteiniyle örtüşen ORF'ye indirgendi: 47-2392).
cub                                  # CUB: 2 gen, L=1226
```

| No | Ölçüm | Metot |
|---|---|---|
| 1 | Kodon sayısı N_kodon (stoplar dahil) | `codon_counts` |
| 2 | Toplam kodon L (stop hariç) | `total_codons` |
| 3 | Amino asit sayısı | `amino_acid_counts` |
| 4 | Kodon frekansı N_kodon / L | `codon_freqs` |
| 5 | Binde frekans | `codon_per_thousand` |

Yardımcılar: `codons(gene)` (okuma sırasıyla kodon listesi), `genes`,
`genetic_code`, `sense_codons`, `codon_families` (amino asit → eş anlamlı
kodonlar), `synonymous_codons`, `fourfold_codons` (standart kodda 32).

```python
name = "MEFV_NM_000243.3"
cub.codon_counts(name)["CUG"]              # 52
cub.total_codons(name)                     # 781
cub.amino_acid_counts(name)["L"]           # 71
cub.codon_per_thousand(name)["CUG"]        # 66.58
```

### 4.3. B. Nükleotit bileşimi — Composition

| No | Ölçüm | Metot |
|---|---|---|
| 6 | GC% | `gc_percent` |
| 7–10 | GC1, GC2, GC12, GC3 | `gc1`, `gc2`, `gc12`, `gc3` |
| 11 | GC3s (eş anlamlı kodonlarda) | `gc3s` |
| 12 | A3s, U3s, G3s, C3s | `base3s` |
| 13 | AT3s | `at3s` |

Temel metotlar `nucleotide_counts(gene, position, codons, bases)` ve
`nucleotide_freqs` diğer bölümlerde de kullanılır. `summary(gene)` hepsini
tek sözlükte verir:

```python
Composition(cub).summary("MEFV_NM_000243.3")
# {'GC%': 59.03, 'GC1': 0.6133, 'GC2': 0.4571, 'GC12': 0.5352, 'GC3': 0.7004,
#  'GC3s': 0.6945, 'A3s': 0.1475, 'U3s': 0.158, 'G3s': 0.3708, 'C3s': 0.3238, 'AT3s': 0.3055}
```

### 4.4. C. Sinonim kodon tercihi — Synonymous

`Synonymous(cub, reference=None, background=None)`

| No | Ölçüm | Metot |
|---|---|---|
| 14 | Kodon oranı N_kodon / N_aa | `codon_ratio` |
| 15 | RSCU | `rscu` |
| 16 | Tercih edilen (RSCU > 1,6) / az kullanılan (< 0,6) kodonlar | `preferred_codons`, `underused_codons` |
| 17 | Fop | `fop`, `optimal_codons` |

Optimal kodonlar iki yöntemle seçilir (`method`):

- `"rscu"`: referansın RSCU'su eşiği (1,6) aşan kodonlar.
- `"chi2"`: Ikemura / CodonW yöntemi. Referans ile `background` arasında 2×2
  χ² testinde anlamlı biçimde fazla kullanılan kodonlar.

Fop paydası da seçilir (`denominator`): `"synonymous"` (tüm eş anlamlı
kodonlar) ya da `"optimal"` (CodonW: yalnızca optimal kodonu olan amino
asitler).

```python
syn = Synonymous(cub, reference=high, background=low)
len(syn.optimal_codons()), len(syn.optimal_codons(method="chi2"))     # (6, 2)
syn.fop(name)                                                        # 0.1632
syn.fop(name, method="chi2", denominator="optimal")                  # 0.2442
```

### 4.5. D. Referans setine dayalı indeksler — Adaptation

`Adaptation(cub, reference=None, background=None)`

| No | Ölçüm | Metot |
|---|---|---|
| 18 | w (göreli uyum) | `relative_adaptiveness(absent_count=0.5)` |
| 19 | CAI (Sharp ve Li 1987) | `cai` |
| 20 | CBI (Bennetzen ve Hall 1982) | `cbi` |

```python
ada = Adaptation(cub, reference=high, background=low)
ada.cai(name)                       # 0.7254
ada.cbi(name, method="chi2")        # -0.1454
```

### 4.6. E. Protein özellikleri (E-A … E-G)

Proteinler de `CUB`'daki genlerden okunur. Böylece aynı izoform için kodon
ve protein ölçümleri yan yana durur. Yedi sınıf vardır ve hepsi `cub`
alır:

| Alt bölüm | Sınıf | Ölçümler |
|---|---|---|
| E-A | `ProteinMass` | `length`, `composition`, `group_composition`, `molecular_weight`, `monoisotopic_mass`, `mean_residue_mass`, `residue_counts`, `disulfide_pairs`, `scale_average` |
| E-B | `ProteinAtoms` | `formula`, `formula_string`, `atom_count`, `atoms_per_residue`, `carbon_oxidation_state` |
| E-C | `ProteinCharge` | `charged_counts`, `simple_net_charge`, `fcr`, `ncpr`, `net_charge(ph)`, `isoelectric_point(pka=…)`, `titration_curve` |
| E-D | `ProteinHydropathy` | `gravy`, `hydropathy(scale=…)`, `aromaticity`, `aliphatic_index`, `boman_index` |
| E-E | `ProteinOptics` | `extinction_coefficient(reduced=…)`, `absorbance_01`, `partial_specific_volume`, `mean_residue_volume`, `mean_residue_surface` |
| E-F | `ProteinStability` | `ivywrel`, `cvp_bias`, `solubility_cv`, `solubility_probability`, `composition_entropy` |
| E-G | `ProteinDisorder` | `uversky`, `top_idp`, `features`, `linear_score(weights=…)` |

```python
ProteinMass(cub).molecular_weight(name)          # 86443.77  (UniProt O15553 ile aynı)
ProteinAtoms(cub).formula_string(name)           # 'C3779H6038N1104O1167S26'
ProteinCharge(cub).isoelectric_point(name)       # 8.3   (Bjellqvist, ExPASy)
ProteinHydropathy(cub).gravy(name)               # -0.645
ProteinOptics(cub).extinction_coefficient(name)  # 62310
ProteinDisorder(cub).uversky(name)
# {'mean_net_charge': 0.0064, 'mean_hydropathy': 0.4283, 'boundary': 0.4156, 'disordered': False}
```

### 4.7. F. Etkin kodon sayısı ve grafikler — EffectiveCodons

`EffectiveCodons(cub)`

| No | Ölçüm | Metot |
|---|---|---|
| 23 | ENC (Wright 1990), F değerleri | `enc`, `homozygosity` |
| 24 | Beklenen ENC: 2 + s + 29 / (s² + (1−s)²) | `expected_enc(gc3s)` |
| 25 | ENC oranı | `enc_ratio`, `enc_points` |
| 26 | Nötrallik: GC12 ~ GC3 regresyonu | `neutrality` |
| 27 | PR2: G3/(G3+C3), A3/(A3+U3) (dört katlı kodonlar) | `pr2`, `pr2_points` |

Literatür biçiminde grafikler (matplotlib):

| Metot | Grafik |
|---|---|
| `plot_enc` | ENC–GC3s: genler ve Wright'ın beklenen eğrisi |
| `plot_enc_ratio` | ENC oranı histogramı |
| `plot_neutrality` | Nötrallik: regresyon doğrusu, denklem, R², y = x |
| `plot_pr2` | PR2: 0–1 eksenler, 0,5 çizgileri |

Her `plot_*` metodu `ax` (var olan bir eksen), `labels` (gen adları) ve
`file_name` (300 dpi kayıt) alır ve `Axes` döndürür.

```python
eff = EffectiveCodons(cub)
eff.enc(name), eff.enc_ratio(name)       # (45.437, 0.1438)
eff.pr2(name)                            # {'x': 0.5303, 'y': 0.5}
eff.neutrality()                         # tek gende slope/intercept/r2 = nan (doğru davranış)
```

### 4.8. G. Entropi ve sapma — Divergence

`Divergence(cub, reference=None)`

| No | Ölçüm | Metot |
|---|---|---|
| 28 | SCUO (Wan 2004) | `scuo` |
| 29 | CDC (kosinüs uzaklığı) | `cdc` |
| 30 | RCDI (Mueller 2006; referans = konak) | `rcdi` |
| 31 | RCBS (Roymondal 2009) | `rcbs` |

Yardımcı: `expected_codon_freqs` (konum baz oranlarından beklenen kodon
oranları).

```python
div = Divergence(cub, reference=CUB(genes=high.genes + low.genes))
div.scuo(name), div.rcdi(name)           # (0.1713, 1.2799)
```

### 4.9. H. tRNA ve kodon bağlamı — CodonContext

`CodonContext(cub, reference=None)`

| No | Ölçüm | Metot |
|---|---|---|
| 32 | tAI (dos Reis 2004) | `trna_weights(trna_counts, s_values, kingdom)`, `tai` |
| 33 | Dinükleotit ρ (Karlin) | `dinucleotide_counts`, `rho(frame=None/"12"/"23"/"31")` |
| 34 | CPS ve CPB (Coleman 2008) | `codon_pair_counts`, `codon_pair_scores`, `cpb` |

tAI için tRNA gen kopya sayıları gerekir. `read_trna_counts` bunları
GtRNAdb FASTA dosyasından okur (insan: `sys/gtrnadb/hg38-tRNAs.fa`, 431
gen, 47 antikodon).

```python
from crn.cub import CodonContext, read_trna_counts
trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
ctx = CodonContext(cub, reference=CUB(genes=high.genes + low.genes))
ctx.tai(name, trna_counts=trna)          # 0.3386
ctx.rho(name)["CG"]                      # 0.4667  (CpG az temsil edilmiş, < 0.78)
ctx.rho(name, frame="31")["CG"]          # 0.4557  (kodon sınırında)
ctx.cpb(name)                            # 0.1649
```

---

## 5. Ayrıntılı kullanım: en etkili biçim

Bu bölüm çerçeveyi gerçek bir çalışma akışı içinde gösterir. Önerilen düzen
şudur:

```
1. Veri al (fetch)  →  2. İzoformlara ayır  →  3. CRN'e kaydet  →
4. Gen setlerini kur (CUB)  →  5. Ölçüm sınıflarını bir kez kur  →
6. Tek gen / çok gen ölçümleri  →  7. Grafikler ve tablolar
```

### 5.1. Santral dogma: bir dizinin yolculuğu

Çerçevenin dizi modeli tek bir örnekte:

```python
from crn import DNA, GenomicStrand, GeneticCode

gene_dna = DNA("CCATGGCTTGCAAATAAGG", name="mini")        # kodlayan zincir
mrna = gene_dna.transcription(template=False, genomic_strand=GenomicStrand.Plus)
orf = gene_dna.find_orfs()[0]                              # ATGGCTTGCAAATAA
protein = orf.transcription(template=False).translation()  # MACK*
cdna = mrna.reverse_transcription()                        # 3'->5' ilk zincir cDNA
same = cdna.complementary()                                # mRNA ile aynı anlamlı DNA
```

Önemli noktalar:

- **Kodlayan zincir:** Gen dizileri kodlayan zincirdir. mRNA'ya geçiş
  `template=False` ile yapılır; `template=True` (varsayılan) kalıp
  zincirden RNA polimeraz gibi tümleyici üretir.
- **`find_orfs`:** Her iki zinciri tarar; ORF'nin konumu ve zinciri
  features'ta yazar.
- **`translation`:** İlk nükleotitten okur. CDS'i önce kesin
  (`rna[start-1:end]`) ya da `match_orf` kullanın.

### 5.2. Bir geni indirmek ve izoformlarına ayırmak

```python
import warnings
from crn.fetch import Ncbi

ncbi = Ncbi(organism="human", verbose=False)
mefv = ncbi.gene("MEFV")                      # ilk çağrıda indirir, sonra sys/'den okur
isoforms = mefv.match_transcripts_to_proteins()
for g in isoforms:
    print(g.name, g.features.first("transcript_id"), g.features.first("protein_id"))
# MEFV_NM_000243.3     NM_000243.3     NP_000234.1
# MEFV_NM_001198536.2  NM_001198536.2  NP_001185465.2
```

GENCODE ile aynı işlem:

```python
from crn.fetch import Gencode
gc = Gencode(verbose=False)
giso = gc.gene("MEFV").match_transcripts_to_proteins()   # 17 transkript → 12 izoform
```

Uyarıları kaydedip incelemek, sessiz bir hesap akışı için iyi bir yoldur:

```python
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    isoforms = ncbi.gene("RPL10").match_transcripts_to_proteins()
for w in caught:
    print(w.category.__name__, w.message)
# UserWarning RPL10: ... aynı CDS'i veren izoformlar çıkarıldı (NM_001303625.1, NM_006013.5 -> NM_001303624.2)
```

### 5.3. Gen panelini kurmak ve CRN ile saklamak

Çok genli analizlerde indirme ve izoform ayırma bir kez yapılmalıdır.
Sonuçlar `.crn` dosyalarına yazılır; `CUB` dosya yollarını da kabul eder.

```python
from pathlib import Path
from crn.cub import CUB

PANEL = Path("data/panel")
PANEL.mkdir(parents=True, exist_ok=True)

def first_nm(symbol):
    """Genin ilk NM_ izoformu (kürasyonlu RefSeq)."""
    return [g for g in ncbi.gene(symbol).match_transcripts_to_proteins()
            if g.features.first("transcript_id").startswith("NM_")][0]

HIGH = ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
        "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")      # yüksek ifadeli
LOW = ("OR1A1", "OR2W1", "OR5A1", "OR7A5", "OR1D2", "OR2J2",
       "OR4D1", "OR6A2", "OR10A2", "OR51E2", "OR2T4", "OR1E1")    # düşük ifadeli

for symbol in HIGH + LOW:
    first_nm(symbol).to_crn(PANEL / f"{symbol}.crn")

# Sonraki oturumlarda ağa ve ORF aramasına gerek kalmadan:
high = CUB(genes=[PANEL / f"{s}.crn" for s in HIGH])
low = CUB(genes=[PANEL / f"{s}.crn" for s in LOW])
```

`CUB` genleri CDS'e indirgedikten sonra kaydederseniz (`g.to_crn` çağrısı
`CUB(...)` sonrası), sonraki okumada "yeniden düzenlendi" uyarısı da
gelmez. Çünkü transkript zaten CDS'tir.

### 5.4. Tek gen için tam CUB profili

Ölçüm sınıflarını **bir kez** kurun ve yeniden kullanın:

```python
from crn.cub import (Composition, Synonymous, Adaptation, EffectiveCodons,
                     Divergence, CodonContext, read_trna_counts)

cub = CUB(genes=isoforms)
name = "MEFV_NM_000243.3"
host = CUB(genes=high.genes + low.genes)

comp = Composition(cub)
syn = Synonymous(cub, reference=high, background=low)
ada = Adaptation(cub, reference=high, background=low)
eff = EffectiveCodons(cub)
div = Divergence(cub, reference=host)
ctx = CodonContext(cub, reference=host)
trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")

profile = {
    "L":     cub.total_codons(name),
    "GC3s":  comp.gc3s(name),
    "Fop":   syn.fop(name, method="chi2", denominator="optimal"),
    "CAI":   ada.cai(name),
    "CBI":   ada.cbi(name, method="chi2"),
    "ENC":   eff.enc(name),
    "ENCr":  eff.enc_ratio(name),
    "SCUO":  div.scuo(name),
    "RCDI":  div.rcdi(name),
    "tAI":   ctx.tai(name, trna_counts=trna),
    "CpG ρ": ctx.rho(name)["CG"],
    "CPB":   ctx.cpb(name),
}
for key, value in profile.items():
    print(f"{key:6} {value:.4f}" if isinstance(value, float) else f"{key:6} {value}")
# L      781
# GC3s   0.6945
# Fop    0.2442
# CAI    0.7254
# CBI    -0.1454
# ENC    45.4370
# ENCr   0.1438
# SCUO   0.1713
# RCDI   1.2799
# tAI    0.3386
# CpG ρ  0.4667
# CPB    0.1649
```

**Nasıl okunur:**

- **ENC:** 45 civarı, orta düzeyde yanlılık (20 aşırı, 61 yanlılık yok).
- **ENC oranı:** 0,14 > 0. Yanlılık GC3s'le beklenenden fazla, yani mutasyon
  baskısı dışında da etken var.
- **CAI:** 0,73. Ribozomal genlerin kodon tercihine orta düzeyde uyum.
- **CBI:** Negatif. χ² ile bulunan optimal kodonlar rastgeleden az
  kullanılıyor.
- **CpG ρ:** 0,47. Omurgalılarda beklenen CpG azlığı.

### 5.5. Referans setinin seçimi

Referansa bağlı ölçümlerde (C, D, G'de RCDI, H'de CPS) sonuç referansa
bağlıdır. Öneriler:

| Ölçüm | Önerilen referans |
|---|---|
| Fop, CBI, CAI (`reference`) | Yüksek ifadeli genler (ribozomal proteinler, EEF1A1) |
| χ² optimal kodonlar (`background`) | Düşük ifadeli genler (ör. koku reseptörleri) |
| RCDI, CPS (`reference`) | Konak genomunun geniş bir ORF seti |

`reference` verilmezse `cub`'ın kendisi kullanılır. Bu durumda tek gende
ölçüm genin kendini kendisiyle karşılaştırmasına dönüşür; anlamı sınırlıdır.
`method="chi2"` için `reference` ile `background` farklı setler olmalıdır.

### 5.6. Çok gen: tablo ve grafikler

Bir panelde her gen için satır üretin:

```python
panel = CUB(genes=[isoforms[0]] + high.genes + low.genes)
comp, eff = Composition(panel), EffectiveCodons(panel)
ada = Adaptation(panel, reference=high)

print(f"{'gen':22} {'GC3s':>6} {'ENC':>6} {'CAI':>6}")
for g in panel.genes:
    print(f"{g.name:22} {comp.gc3s(g.name):6.3f} {eff.enc(g.name):6.2f} {ada.cai(g.name):6.3f}")
# MEFV_NM_000243.3        0.695  45.44  0.725
# RPL3_NM_000967.4        0.711  44.68  0.800
# RPL4_NM_000968.4        0.460  54.76  0.711
# ...
# OR1E1_NM_003553.3       0.667  47.02  0.716
```

Dört F grafiğini tek şekilde çizin:

```python
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, figsize=(11, 10))
eff.plot_enc(ax=axs[0, 0])
eff.plot_enc_ratio(ax=axs[0, 1])
eff.plot_neutrality(ax=axs[1, 0], labels=True)
eff.plot_pr2(ax=axs[1, 1])
fig.tight_layout()
fig.savefig("f_grafikleri.png", dpi=300)
```

Tek grafik için: `eff.plot_enc(file_name="enc.png")`.

**Paneldeki sonuçlar (25 gen):**

- **ENC–GC3s:** Genler Wright eğrisinin altında. Yanlılık yalnızca GC3s ile
  açıklanmıyor.
- **Nötrallik:** Eğim 0,129, R² = 0,07. GC3'teki değişimin ancak yaklaşık
  %13'ü mutasyon baskısıyla açıklanıyor; seçilim baskın.
- **PR2:** Noktalar sol alt çeyrekte. Dört katlı 3. konumlarda C > G ve
  U > A kullanılıyor.

Ham veriye de erişebilirsiniz: `eff.enc_points()`, `eff.neutrality()`,
`eff.pr2_points()`. Bu sözlük listeleri pandas'a (`pd.DataFrame(...)`)
doğrudan aktarılabilir.

### 5.7. Protein profili

```python
from crn.cub import (ProteinMass, ProteinAtoms, ProteinCharge, ProteinHydropathy,
                     ProteinOptics, ProteinStability, ProteinDisorder)

mass, charge = ProteinMass(cub), ProteinCharge(cub)
hyd, opt = ProteinHydropathy(cub), ProteinOptics(cub)

print(mass.length(name), round(mass.molecular_weight(name), 2))   # 781 86443.77
print(ProteinAtoms(cub).formula_string(name))                     # C3779H6038N1104O1167S26
print(round(charge.isoelectric_point(name), 2))                   # 8.3
print(round(charge.isoelectric_point(name, pka="EMBOSS"), 2))     # başka pKa kümesi
print(charge.net_charge(name, ph=7.4))
print(round(hyd.gravy(name), 3))                                  # -0.645
print(hyd.hydropathy(name, scale="Hopp-Woods"))
print(opt.extinction_coefficient(name), opt.extinction_coefficient(name, reduced=True))
print(ProteinStability(cub).ivywrel(name))                        # 0.373
curve = charge.titration_curve(name)                              # [(pH, Q), ...]
```

### 5.8. tRNA havuzu ve kodon bağlamı

```python
trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
w = ctx.trna_weights(trna)                         # {kodon: w}
max(w, key=w.get), min(w, key=w.get)               # ('UGC', 'CUA')
ctx.tai(name, trna_counts=trna)                    # ökaryot (varsayılan)
ctx.tai(name, trna_counts=bakteri_trna, kingdom="prokaryote")   # AUA için lizidin kuralı

rho = ctx.rho(name)                                # 16 dinükleotit
{k: round(v, 2) for k, v in rho.items() if v < 0.78 or v > 1.23}   # az/çok temsil edilenler
ctx.dinucleotide_counts(name, frame="31")          # kodon sınırındaki çiftler
scores = ctx.codon_pair_scores()                   # referanstaki tüm çiftlerin CPS'si
```

Panelde tAI: yüksek ifadeli ortalama 0,355, koku reseptörleri 0,332.
İnsanda tRNA havuzuna uyumun etkisi küçüktür.

### 5.9. Başka genetik kod

Genetik kod bir kez, `CUB` kurulurken seçilir. Bütün aileler, stoplar,
dört katlı kodonlar ve ENC sınıfları bu koda göre kurulur:

```python
from crn import GeneticCode
mt = CUB(genes=mito_genes, genetic_code=GeneticCode.VERTEBRATE_MITOCHONDRIAL)
```

Nükleer bir geni mitokondriyal kodla vermek ORF eşleşmesini bozar ve
`ValueError` ("örtüşen ORF yok") verir. Bu doğru davranıştır.

### 5.10. Kendi verinizi kullanmak

`CUB` yalnızca `Gene` ister; veri `fetch`'ten gelmek zorunda değildir:

```python
from crn import Gene, DNA, RNA, PRT

cds = "AUGGCUGCCGCAGCGUAA"
g = Gene("sentetik", DNA(cds.replace("U", "T")))
g.add_transcript(RNA(cds), "t1")
g.add_protein(PRT("MAAAA"), "p1")
cub = CUB(genes=[g])
```

FASTA'dan: `RNA.from_fasta_file("cds.fa")` ve `PRT.from_fasta_file("prot.fa")`
ile okuyup her çift için bir `Gene` kurun.

---

## 6. İyi uygulamalar ve sık hatalar

**İyi uygulamalar**

1. **Her zaman izoformlara ayırın.** `CUB`'a `ncbi.gene(...)` sonucunu
   doğrudan değil, `match_transcripts_to_proteins()` çıktısını verin.
2. **`unique_cds=True`'yu kapatmayın.** Aynı CDS'li kopyalar referans
   setlerini, regresyonu ve CPS'yi çarpıtır.
3. **İndirmeyi ve ORF aramasını bir kez yapın**, sonucu `.crn` olarak
   saklayın (§5.3).
4. **Ölçüm sınıflarını bir kez kurun**, döngüde yalnızca metot çağırın.
5. **Referansları açıkça verin.** Varsayılan (`cub`'ın kendisi) tek gende
   anlamsız sonuç verir.
6. **Tek gen listesi:** `genes=[gene]`. `genes=gene` `TypeError` verir.
7. **Uyarıları okuyun.** `CubWarning` gen düzenlemesini, eksik referans
   kodonlarını (RCDI) ve çiftlerini (CPB) bildirir. Seçici susturmak için:
   `warnings.filterwarnings("ignore", category=CubWarning)`.
8. **API anahtarını ortam değişkeninde tutun** (`NCBI_API_KEY`).
9. **Tek kaynak kullanın:** NCBI ile GENCODE'u aynı analizde karıştırmayın.

**Sık hatalar**

| Belirti | Neden | Çözüm |
|---|---|---|
| `ValueError: Her gende bir genomik DNA, bir transkript ve bir protein olmalı` | Gen bölünmemiş | `match_transcripts_to_proteins()` |
| `TypeError: genes bir liste olmalı` | `genes=gene` | `genes=[gene]` |
| `ValueError: … örtüşen ORF yok` | Protein transkriptten çevrilemiyor (yanlış genetik kod, kodlamayan RNA) | Genetik kodu ve kaydı kontrol edin |
| `ValueError: DNA dizisinde geçersiz karakter(ler)` | Alfabede olmayan harf | Diziyi temizleyin ya da doğru sınıfı kullanın (RNA'da U) |
| Nötrallik `nan` | İkiden az gen | Çok genli panel kullanın |
| ENC `nan` | Bir dejenerelik sınıfı için hiç veri yok (çok kısa gen) | Daha uzun gen ya da gen seti |
| `KeyError: CUB: '…' geni yok` | Gen adı yanlış | Adlar: `[g.name for g in cub.genes]` |
| `ValueError: trna_counts verilmeli` | tAI için tRNA verisi yok | `read_trna_counts(...)` |

---

## 7. Sınırlılıklar

- **ORF eşleme** yalnızca transkriptin (+) zincirinde ve tam eşleşmeyle
  yapılır. Selenoproteinler (UGA = Sec) ve okunan stop kodonları eşleşmez.
- **Genetik kodlar:** yalnızca standart ve omurgalı mitokondriyal kod
  tanımlıdır. Yenisi `GeneticCode`'a eklenebilir.
- **Kodon sayımları önbelleğe alınmaz.** Her ölçüm sayımları yeniden yapar.
  Yüzlerce gende sorun olmaz; genom ölçeğinde (~20.000 gen) yavaşlık
  beklenir.
- **tAI** gen kopya sayısını tRNA bolluğunun vekili kabul eder. s değerleri
  maya için optimize edilmiştir; tür özgü s (stAI) yoktur.
- **CDC** Zhang ve ark.'nın bootstrap anlamlılık testini içermez. Beklenen
  bileşim doğrudan konum baz oranlarından kurulur.
- **Wilkinson–Harrison** sabitleri Davis ve ark. (1999) yeniden
  kalibrasyonundandır ve E. coli'de aşırı ifade için geçerlidir. Özgün
  makaleyle ayrıca doğrulanmalıdır.
- **Protein-Sol gibi skorlar:** yayımlanmış model ağırlıkları pakette yoktur.
  `ProteinDisorder.linear_score`'a dışarıdan verilir.
- **NCBI ve GENCODE** farklı transkript kümeleri ve koordinatlar verir;
  sonuçlar kaynağa göre değişir.
