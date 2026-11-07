# crn: Nükleik asit ve protein dizilerinin modellenmesi ile kodon kullanım yanlılığının çözümlenmesi için bir Python çerçevesi

`crn`, moleküler biyolojinin santral dogmasında yer alan DNA, RNA ve protein düzeylerinin tipli nesnelerle temsil edildiği, kamuya açık genomik kaynaklardan elde edilen verilerin tek bir gen modeli altında birleştirildiği ve kodon kullanım yanlılığının (codon usage bias, CUB) genden gen grubuna uzanan ölçeklerde nicelendirildiği bir biyoinformatik çerçevedir. Çerçevede dizi nesneleri, veri erişim katmanı ve ölçüm katmanı birbirinden ayrılmıştır; ölçümler, verinin hangi kaynaktan alındığından bağımsız olarak `Gene` nesneleri üzerinden yürütülür. Literatürde yaygın olarak kullanılan CUB ölçümlerine ek olarak iki gen grubunun kodon düzeyinde karşılaştırılması için bir diferansiyel kodon kullanımı yöntemi sunulmuştur. Bu yöntemde Fisher kesin testi, gen düzeyinde permütasyon testiyle birlikte değerlendirilir; çünkü kodon sayımlarının havuzlanması durumunda Fisher testinin rastgele bölünmüş gen gruplarında dahi yüksek oranda anlamlı sonuç ürettiği gözlenmiştir (bkz. [Bölüm 8](#8-diferansiyel-kodon-kullanımı)).

## İçindekiler

1. [Amaç ve kapsam](#1-amaç-ve-kapsam)
2. [Mimari](#2-mimari)
3. [Tasarım ilkeleri](#3-tasarım-ilkeleri)
4. [Kurulum](#4-kurulum)
5. [Veri modeli](#5-veri-modeli)
6. [Veri kaynakları](#6-veri-kaynakları)
7. [Kodon kullanımı ölçümleri](#7-kodon-kullanımı-ölçümleri)
8. [Diferansiyel kodon kullanımı](#8-diferansiyel-kodon-kullanımı)
9. [Protein özellikleri](#9-protein-özellikleri)
10. [Örnek iş akışı](#10-örnek-iş-akışı)
11. [Doğrulama ve yinelenebilirlik](#11-doğrulama-ve-yinelenebilirlik)
12. [Sınırlılıklar](#12-sınırlılıklar)
13. [Belgeler](#13-belgeler)
14. [Kaynakça](#14-kaynakça)

[Ek A. Önceki arayüzden geçiş](#ek-a-önceki-arayüzden-geçiş)

## 1. Amaç ve kapsam

Kodon kullanım yanlılığı çalışmalarında bir sonucun güvenilirliği, ölçüm formülünden önce, hangi transkriptin hangi proteinle eşleştirildiğine ve sayımların hangi birimde yapıldığına bağlıdır. Çerçevenin tasarımında bu nedenle her ölçümün proteinle doğrulanmış bir kodlama dizisi (CDS) üzerinde hesaplanması esas alınmıştır: bir transkriptin ölçüme girebilmesi için çevirisi kaynak veritabanındaki proteinle birebir örtüşen bir açık okuma çerçevesi taşıması gerekir.

Çerçevenin kapsamı üç işlevle sınırlandırılmıştır. Birinci olarak DNA, RNA ve protein dizileri alfabe denetimi yapan ayrı sınıflarla modellenmiş; transkripsiyon, ters transkripsiyon, çeviri ve açık okuma çerçevesi taraması bu sınıfların metotları olarak tanımlanmıştır. İkinci olarak NCBI RefSeq, GENCODE, HGNC ve Human Phenotype Ontology (HPO) kaynaklarına erişim sağlanmış, indirilen veriler yerel bir depoda saklanarak sonraki sorgularda ağ erişimi gereksinimi ortadan kaldırılmıştır. Üçüncü olarak kodon kullanımı, nükleotit bileşimi, tRNA uyumu, kodon bağlamı ve protein özelliklerine ilişkin ölçümler ile iki gen grubunun karşılaştırılmasına yönelik istatistiksel yöntemler uygulanmıştır.

Dizi hizalama, varyant çağırma, ifade düzeyi nicelemesi ve filogenetik çıkarım çerçevenin kapsamı dışında bırakılmıştır.

## 2. Mimari

Paket, veri modeli, veri erişimi ve ölçüm olmak üzere üç katmandan oluşturulmuştur.

```text
crn/
├── constants.py   alfabeler, kodonlar, tümleme tabloları, GeneticCode
├── tables.py      kütleler, pKa kümeleri, hidropati ölçekleri, tAI cezaları
├── features.py    Features: anahtar -> değer listesi deposu
├── seq.py         Sequence, NucleicAcid, DNA, RNA, PRT; Molecule, Strand, GenomicStrand
├── gene.py        Gene: genomik DNA, transkriptler, proteinler; isoforms()
├── fileio.py      FASTA ve CRN dosya biçimleri
├── protein.py     ProteinProperties
├── fetch/         Ncbi, NcbiSession, Gencode, Hgnc, Hpo; yerel sys/ deposu
└── cub/           CUB, ölçüm sınıfları, diferansiyel analiz, uyum çözümlemesi, raporlama
```

Katmanlar arasındaki bağımlılık tek yönlüdür. `fetch` ve `cub` paketleri veri modelini kullanır, veri modeli ise bu paketlerden hiçbirine bağımlı değildir. `cub` paketinin `fetch` paketine bağımlılığı da bulunmamaktadır; bu sayede ölçüm katmanına kullanıcı tarafından kurulmuş sentetik genler ile veritabanlarından alınmış genler aynı biçimde verilebilir. Verinin çerçeve içindeki akışı aşağıda gösterilmiştir.

```text
  HPO ── semptom → gen listesi ──┐
  HGNC ── aile, ortolog ─────────┤
                                 ▼
  NCBI RefSeq ───┐    Ncbi().gene(ad)            ┌── FASTA / .crn dosyası
  GENCODE ───────┤    Gencode().gene(ad)         │   read_fasta, Gene.from_crn
                 ▼                               ▼
          ┌─────────────────────────────────────────┐
          │ Gene                                    │
          │  ├── dna          genomik DNA (1)       │
          │  ├── transcripts  {kimlik: RNA} (n)     │
          │  ├── proteins     {kimlik: PRT} (m)     │
          │  └── features     gene_id, aliases ...  │
          └─────────────────────────────────────────┘
                 │  gene.isoforms()
                 ▼
          [Gene, Gene, ...]     her biri 1 DNA, 1 transkript, 1 protein
                 │  CUB(izoformlar)   transkript → proteinle doğrulanmış CDS
                 ▼
          CUB ──► Composition, Synonymous, Adaptation, EffectiveCodons,
             │    Divergence, CodonContext, CorrespondenceAnalysis,
             │    ProteinProperties
             │
             ├──► cub.table()                 gen başına ölçüm tablosu (TSV)
             └──► DifferentialCodonUsage      iki grubun kodon düzeyinde karşılaştırılması
```

## 3. Tasarım ilkeleri

Çerçevenin kurallarından beş tanesi, kullanıcı tarafından gözlenen davranışı doğrudan belirlediği için burada açıklanmıştır.

Bir diziye ilişkin kimlik, kaynak ve koordinat gibi ek bilgiler sınıf alanı olarak değil, `Features` adı verilen anahtar-değer listesi deposunda tutulur; dolayısıyla yeni bir bilginin taşınması için sınıf tanımının değiştirilmesi gerekmez. Bu kuralın istisnası, dizinin biyolojik yorumunu değiştiren alanlardır: yazılış yönünü belirten `strand` ile RNA'da genomik zinciri belirten `genomic_strand` sınıf düzeyinde tanımlanmıştır.

Diziler nesne oluşturulurken büyük harfe çevrilir ve alfabeye göre doğrulanır. Geçersiz karakter içeren bir dizi `ValueError` ile reddedilir. Belirsiz durumlar sessizce geçilmez; işlem durdurulur ya da `UserWarning` veya `CubWarning` ile bildirilir.

Yeni nesne döndüren işlemlerde girdi nesneleri değiştirilmez. Dilimleme, `complementary()` ve `transcription()` gibi işlemler `features` deposunun kopyasını taşır; `CUB` sınıfı da kendisine verilen genleri kırpmaz, CDS'e indirgenmiş kopyalar üzerinde çalışır.

Sayısal sabitler `constants.py` ve `tables.py` dosyalarında toplanmıştır ve ilgili fonksiyonlara parametre olarak verilebilir. Böylece bir pKa kümesinin ya da hidropati ölçeğinin değiştirilmesi için kaynak kodda düzenleme yapılması gerekmez.

Maliyetli hesaplar tek sefer yapılır. Kodon sayımları `CUB` kurulurken gen başına bir kez üretilir ve ölçümler bu 64 kodonluk sayımdan türetilir; optimal kodonlar, göreli uyum değerleri ve kodon çifti skorları ilk kullanımda hesaplanıp saklanır. Bu sayede ölçümler genom ölçeğindeki gen panellerine yeniden sayım yapılmaksızın uygulanabilir.

## 4. Kurulum

Çerçeve Python 3.12 ve üzeri sürümleri gerektirir. Bağımlılıklar `tqdm`, `requests`, `matplotlib`, `scipy` (1.18.0 ve üzeri) ve `numpy` (2.0 ve üzeri) paketlerinden oluşur. `matplotlib` yalnızca grafik üretilirken, `scipy` ve `numpy` yalnızca diferansiyel kodon testinde yüklenir. Paket, proje kökünde aşağıdaki komutlarla kurulur.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Veriler varsayılan olarak proje kökündeki `sys/` klasörüne yazılır. Bu konum `CRN_SYS_DIR` ortam değişkeniyle ya da veri sınıflarının `sys_dir` parametresiyle değiştirilebilir. NCBI API anahtarı kullanılacaksa anahtarın kaynak koda yazılmaması, `NCBI_API_KEY` ortam değişkeninde tutulması önerilir.

## 5. Veri modeli

Dizi sınıfları ortak `Sequence` tabanından türetilmiştir; `DNA` ve `RNA` sınıfları ayrıca `NucleicAcid` tabanını paylaşır. Sınıfların alfabeleri ve başlıca işlemleri Tablo 1'de verilmiştir.

Tablo 1. Dizi sınıfları.

| Sınıf | Alfabe | Sınıf alanları | Başlıca işlemler |
|---|---|---|---|
| `Sequence` | Denetlenmez | | `seq`, `name`, `features`, dilimleme, FASTA ve CRN okuma-yazma |
| `DNA` | `ACGT`, IUPAC kodları, `-` | `strand` | `complement`, `complementary`, `reverse_complement`, `transcription`, `find_orfs` |
| `RNA` | `ACGU`, IUPAC kodları, `-` | `strand`, `genomic_strand` | `complement`, `complementary`, `reverse_complement`, `reverse_transcription`, `translation`, `match_orf` |
| `PRT` | 20 amino asit, `BZJXUO*-` | | `seq`, `name`, `features` |

`strand` alanı dizinin yazılış yönünü ifade eder. `complement()` ve `complementary()` metotları tümleyici zinciri özgün dizinin altına hizalı olarak 3'→5' yönünde, `reverse_complement()` metodu ise karşı zinciri 5'→3' yönünde üretir.

`Gene` sınıfında bir genomik DNA ile istenen sayıda transkript ve protein bir arada tutulur. Kodon düzeyindeki ölçümler izoform başına tanımlı olduğundan, `gene.isoforms()` metodu ile her transkript, çevirisi bir proteinle birebir örtüşen ORF üzerinden o proteinle eşleştirilir ve bir DNA, bir transkript ve bir protein içeren yeni `Gene` nesneleri üretilir. Eşleştirmede önce kaynağın bildirdiği transkript-protein bağı denenir. Aynı CDS'i veren kopyalardan yalnızca biri tutulur, eşleşmeyen transkriptler `UserWarning` ile raporlanır.

Veriler iki dosya biçiminde saklanır. FASTA biçiminde yalnızca ad ve dizi korunur; tip, yön ve `features` bilgileri yitirilir. CRN biçimi ise pickle kullanılmaksızın tanımlanmış, CRC32 ile korunan ikili bir biçimdir ve dizi, gen ile bunların listelerini kayıpsız olarak saklar. Yorumlanabilir nesne yüklenmediği için `.crn` dosyalarının okunması sırasında kod yürütülmez.

## 6. Veri kaynakları

Bütün kaynaklarda aynı erişim düzeni izlenmiştir: istenen veri önce yerel depoda aranır, bulunamazsa indirilir ve saklanır. Kaynakların özellikleri Tablo 2'de özetlenmiştir.

Tablo 2. Veri kaynakları.

| Sınıf | Kaynak | İçerik | İndirme birimi |
|---|---|---|---|
| `Ncbi` | NCBI E-utilities, RefSeq | Genomik DNA, transkriptler, proteinler; altı organizma | Gen başına; 3 istek/sn, API anahtarıyla 10 istek/sn |
| `Gencode` | GENCODE sürümleri | Genom, transkriptler, proteinler, GTF; insan ve fare | Sürüm başına; 1,3 GB indirme, açılmış hâlde 5 GB |
| `Hgnc` | HGNC, HCOP | Onaylı gen sembolleri, gen aileleri, 10 türün ortologları | Tablo olarak |
| `Hpo` | Human Phenotype Ontology | Ontoloji, gen-fenotip ve hastalık notları | Sürüm başına; yaklaşık 66 MB, SHA-256 ile doğrulanır |

`Ncbi` ve `Gencode` sınıflarında aynı adlı `gene`, `dna`, `rna`, `protein`, `region` ve `update` metotları tanımlanmıştır. İki kaynak aynı gen için farklı sınırlar, transkript kümeleri ve kimlikler verdiğinden, tek bir çözümlemede kaynakların karıştırılmaması gerekir. `Gencode.canonical` metodu geni yalnızca ana transkripti ve proteiniyle döndürür; ana transkript `MANE_Select` etiketinden, bu etiket yoksa `Ensembl_canonical` etiketinden belirlenir. Tek gene yönelik sorgularda `Ncbi`, genom ölçeğindeki çözümlemelerde ise ağ hız sınırına tabi olmayan `Gencode` tercih edilmelidir.

`Hpo` sınıfı, bir fenotip teriminden gen listesi elde edilmesini sağlar. Terimler alt terimleriyle birlikte sorgulanır ve sonuçlar sıklık, ilişki türü, kanıt kodu, başlangıç yaşı ve kalıtım biçimine göre süzülebilir. Bu belgede esas alınan v2026-09-01 sürümü 19894 terim, 5283 gen, 12880 hastalık ve 330074 not içermektedir. Hastalıkta bulunmadığı açıkça belirtilen 826 `NOT` kaydı ile hiç gözlenmediği bildirilen 3083 `0/N` kaydı yükleme sırasında dışlanır. HPO'da gen ile fenotip arasındaki bağ doğrudan değildir, hastalık üzerinden kurulur; elde edilen gen listeleri bu çerçevede yorumlanmalıdır.

## 7. Kodon kullanımı ölçümleri

`CUB` sınıfı, her biri bir DNA, bir transkript ve bir protein içeren `Gene` listesi ya da `.crn` dosya yollarıyla kurulur. Kurulum sırasında proteinle örtüşen CDS belirlenir ve kodonlar sayılır. Metotlarda `gene` parametresi verilmediğinde tüm genler birlikte, verildiğinde yalnızca ilgili gen hesaplanır. Standart genetik kod ile omurgalı mitokondriyal kod `GeneticCode` üzerinden seçilebilir. Ölçümler Tablo 3'te bölümlere ayrılmış olarak listelenmiştir.

Tablo 3. Ölçüm sınıfları.

| Bölüm | Sınıf | Ölçümler | Kaynak |
|---|---|---|---|
| A | `CUB` | Kodon ve amino asit sayıları, frekanslar, binde kullanım | |
| B | `Composition` | GC%, GC1, GC2, GC12, GC3, GC3s, A3s, U3s, G3s, C3s | |
| C | `Synonymous` | Kodon oranı, RSCU, tercih edilen ve optimal kodonlar, Fop | Ikemura, 1981 |
| D | `Adaptation` | Göreli uyum (w), CAI, CBI | Sharp ve Li, 1987; Bennetzen ve Hall, 1982 |
| E | `ProteinProperties` | Kütle, formül, yük, hidropati, optik, stabilite, düzensizlik | |
| F | `EffectiveCodons` | ENC, beklenen ENC, nötrallik analizi, PR2 | Wright, 1990; Sueoka, 1988; Sueoka, 1995 |
| G | `Divergence` | SCUO, CDC, RCDI, RCBS | Wan ve ark., 2004; Zhang ve ark., 2012; Mueller ve ark., 2006; Roymondal ve ark., 2009 |
| H | `CodonContext` | tAI, dinükleotit göreli bolluğu (ρ), CPS, CPB | dos Reis ve ark., 2004; Karlin ve Burge, 1995; Coleman ve ark., 2008 |
| J | `DifferentialCodonUsage` | İki grubun kodon düzeyinde karşılaştırılması | Benjamini ve Hochberg, 1995 |
| | `CorrespondenceAnalysis` | RSCU matrisinin uyum çözümlemesi | |

Her ölçümün formülü ve özgün kaynakla ilişkisi `crn/cub/` modüllerinin başındaki tablolarda ve [doc/cub/](doc/cub/README.md) belgelerinde verilmiştir.

Referansa dayalı ölçümlerde sonuç, referans setinin seçimine doğrudan bağlıdır. Fop, CAI ve CBI için yüksek ifadeli genlerden (örneğin ribozomal proteinler ve EEF1A1) oluşan bir set, ki-kare temelli optimal kodon belirlemesinde buna ek olarak düşük ifadeli bir arka plan seti kullanılmalıdır. RCDI için konak genomunun, CPS ve CPB için ise genomun geniş bir ORF setinin verilmesi gerekir. Referans belirtilmediğinde `cub` nesnesinin kendisi kullanılır; tek genlik bir nesnede bu durum genin kendisiyle karşılaştırılması anlamına gelir. tAI hesabında tRNA gen kopya sayıları, GtRNAdb'den alınan FASTA dosyasından `read_trna_counts` ile okunur.

Bir gen panelinin bütün skaler ölçümleri `cub.table()` metoduyla gen başına bir satır olarak üretilir ve `file_name` parametresiyle TSV biçiminde yazılır. Verilen referanslara göre tabloya CAI, RCDI, CPB ve tAI gibi sütunlar eklenir. `correlation_matrix` ve `plot_correlations` fonksiyonlarıyla ölçümler arasındaki korelasyonlar incelenebilir. ENC-GC3s, ENC oranı, nötrallik, PR2, RSCU, CAI dağılımı, uyum çözümlemesi ve volkan grafikleri için çizim metotları tanımlanmıştır; tamamı `ax` ve `file_name` parametrelerini kabul eder.

## 8. Diferansiyel kodon kullanımı

`DifferentialCodonUsage` sınıfında iki gen grubu kodon başına kurulan 2×2 olumsallık tablolarıyla karşılaştırılır. Tabloda ilgili kodonun sayımı, karşılaştırma havuzunun geri kalanıyla birlikte değerlendirilir. Havuzun tanımı üç düzeyde yapılmıştır: `global` düzeyinde bütün anlamlı kodonlar, `family` düzeyinde aynı amino asidin kodonları, `subfamily` düzeyinde aynı amino asidi kodlayan ve ilk iki bazı ortak olan kodonlar kullanılır. `global` düzeyindeki sonuçlar amino asit bileşimi farkıyla karışabildiği için sinonim tercihin esas ölçütü olarak `family` düzeyi alınmıştır. Stop kodonları hiçbir düzeye dâhil edilmez, tek kodonlu amino asitler aile ve alt aile düzeylerinde test edilmez.

Etki büyüklüğü Haldane-Anscombe düzeltmesi uygulanmış olasılık oranının ikili logaritmasıyla (`log2_or`) ifade edilir. Anlamlılık iki yönlü Fisher kesin testiyle hesaplanır ve çoklu test düzeltmesi Benjamini-Hochberg yöntemiyle (Benjamini ve Hochberg, 1995) yapılır. Genom ölçeğindeki sayımlarda p değerleri çift duyarlıklı kayan nokta aralığının altına inebildiğinden, 1e-200 eşiğinin altındaki değerler log uzayında yeniden hesaplanır.

Fisher testi grup içi sayımları havuzlar ve her kodonu bağımsız gözlem kabul eder. Oysa genler kodon kullanımı bakımından kendi aralarında farklılaştığı için bu varsayım geçerli değildir. Söz konusu durum, kontrol havuzunun 70'er genlik gruplara 12 kez rastgele bölünmesiyle sınanmıştır (`test/dco/dco_6.py`). Gruplar arasında gerçek bir fark bulunmamasına karşın Fisher testinde aile düzeyindeki 59 kodonun 53'üne kadarı anlamlı bulunmuş, aynı bölünmelerde gen düzeyinde permütasyon testiyle hiçbir kodon anlamlı bulunmamıştır. Gerçek fark eklenen pozitif kontrolde permütasyon testiyle 59 kodonun tamamı saptanmış, sıfır hipotezi altında `perm_p < 0,05` oranı %6,05 olarak ölçülmüştür.

Bu bulgular doğrultusunda `permutation` metodunda gen etiketleri karıştırılarak `log2_or` değerinin sıfır dağılımı oluşturulur ve `robust` metodunda yalnızca iki testten de geçen kodonlar raporlanır. `interpret` metodunda üretilen yorum da bu listeye dayanır. Permütasyonda elde edilebilecek en küçük p değeri `1 / (n + 1)` olduğundan, karıştırma sayısı `n`, test sayısının `1 / q` katından az seçilmemelidir; aile düzeyinde 59 test ve q = 0,05 için bu alt sınır 1180'dir. Grup büyüklüğü de ulaşılabilecek en küçük p değerini sınırlar: 4 gene karşı 4 genlik bir karşılaştırmada ayrık bölünme sayısı 35 olduğundan en küçük p değeri 1/35'tir.

Permütasyon testi genler arası çeşitliliği hesaba katar, ancak gruplar arasındaki gen uzunluğu ve bölgesel GC farkını gidermez. GC3s bakımından farklılaşan iki grupta saptanan kodonlar çeviri seçiliminden değil, mutasyon yanlılığından kaynaklanıyor olabilir. Kontrol grubunun bu nedenle uzunluk ve GC içeriğine göre eşleştirilmesi gerekir; eşleştirmenin bir uygulaması `test/dco/dco_ortak.py` dosyasında verilmiştir.

## 9. Protein özellikleri

`ProteinProperties` sınıfı bir `CUB` nesnesi, bir `PRT` nesnesi ya da `PRT` listesiyle kurulabilir ve iki yoldan aynı protein için aynı sonuç elde edilir. Hesaplanan özellikler yedi grupta toplanmıştır (Tablo 4). Birden çok protein birlikte değerlendirildiğinde her biri ayrı zincir sayılır.

Tablo 4. Protein özellikleri.

| Bölüm | Kapsam | Metotlar |
|---|---|---|
| E-A | Sayım ve kütle | `length`, `composition`, `molecular_weight`, `monoisotopic_mass`, `residue_counts`, `disulfide_pairs` |
| E-B | Atom bileşimi | `formula`, `formula_string`, `atom_count`, `carbon_oxidation_state` |
| E-C | Yük | `charged_counts`, `fcr`, `ncpr`, `net_charge`, `isoelectric_point`, `titration_curve` |
| E-D | Hidropati | `gravy`, `hydropathy`, `aromaticity`, `aliphatic_index`, `boman_index` |
| E-E | Optik ve hacim | `extinction_coefficient`, `absorbance_01`, `partial_specific_volume`, `mean_residue_volume` |
| E-F | Stabilite ve çözünürlük | `ivywrel`, `cvp_bias`, `solubility_cv`, `solubility_probability`, `composition_entropy` |
| E-G | Düzensizlik | `uversky`, `top_idp`, `composition_features`, `linear_score` |

Metotların tam listesi ve kullanılan ölçekler [doc/protein.md](doc/protein.md) ile [doc/tables.md](doc/tables.md) belgelerinde verilmiştir.

## 10. Örnek iş akışı

Aşağıdaki kod parçalarında bir genin izoformlara ayrılması, referans setine göre ölçülmesi ve iki gen grubunun karşılaştırılması gösterilmiştir. Yorum satırlarındaki değerler kodun çalıştırılmasıyla elde edilmiştir ve kullanılan veri sürümüne bağlıdır.

```python
from crn.fetch import Ncbi
from crn.cub import CUB, Adaptation, Composition, EffectiveCodons

ncbi = Ncbi(verbose=False)
mefv = ncbi.gene("MEFV")
print(mefv)  # -> MEFV: DNA <14607>, 2 transkript, 2 protein

izoformlar = mefv.isoforms()
cub = CUB(izoformlar)
print(cub.names)  # -> ['MEFV_NM_000243.3', 'MEFV_NM_001198536.2']
print(cub)  # -> CUB: 2 gen, L=1226
```

Referans seti, yüksek ifadeli genlerin birinci RefSeq izoformlarından kurulur.

```python
def ilk_nm(sembol):
    return next(g for g in ncbi.gene(sembol).isoforms()
                if g.features.first("transcript_id").startswith("NM_"))

yuksek = CUB([ilk_nm(s) for s in ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
                                  "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")])

ad = "MEFV_NM_000243.3"
print(round(Adaptation(cub, reference=yuksek).cai(ad), 4))  # -> 0.7254
print(round(EffectiveCodons(cub).enc(ad), 2))  # -> 45.44
print(round(Composition(cub).gc3s(ad), 3))  # -> 0.695
```

İki grubun karşılaştırılmasında GENCODE ana transkriptleri kullanılmıştır. Örnekte ribozomal ve koku reseptörü genleri, grupları tanıtmak amacıyla küçük tutulmuştur; bu gruplar uzunluk ve GC içeriğine göre eşleştirilmemiştir ve sonuç biyolojik bir bulgu olarak değerlendirilmemelidir.

```python
from crn.fetch import Gencode
from crn.cub import DifferentialCodonUsage

gc = Gencode(verbose=False)

def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

ribozomal = kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2", "RPL5", "RPS8"])
koku = kur(["OR1A1", "OR2A2", "OR5A1", "OR10A2", "OR4D1", "OR6A2", "OR1E1", "OR2B2"])

dcu = DifferentialCodonUsage(ribozomal, koku, name_a="ribozom", name_b="koku")
print(dcu)  # -> DifferentialCodonUsage: ribozom (8 gen) / koku (8 gen)
print(len(dcu.significant("family")), len(dcu.robust("family", n=2000)))  # -> 4 2
```

Semptomdan gen grubuna geçiş `Hpo` sınıfıyla yapılır. `pack` metodu ile gruplar tek bir CRN dosyasına yazılır ve `DifferentialCodonUsage.from_crn` ile okunur.

```python
from crn.fetch import Hpo

hpo = Hpo(verbose=False)
print(len(hpo.symbols("Ataxia")), len(hpo.symbols("Ataxia", min_frequency=0.3)))  # -> 1119 669
```

<!-- calistirma -->
```python
rapor = hpo.pack({"nobet": "Seizure", "ataksi": "Ataxia"}, "gruplar.crn")
dcu = DifferentialCodonUsage.from_crn("gruplar.crn", "nobet", "ataksi")
print(dcu.interpret())
```

Modüllerin birlikte kullanıldığı diğer uçtan uca senaryolar [doc/senaryolar.md](doc/senaryolar.md) belgesinde sunulmuştur.

## 11. Doğrulama ve yinelenebilirlik

Çerçevenin kullanımı `test/` klasöründeki çalıştırılabilir betiklerle örneklenmiştir: dizi işlemleri `test/sequence/`, veri kaynakları `test/ncbi/` ve `test/gencode/`, fenotip temelli gen grupları `test/hpo/`, diferansiyel analiz ise `test/cub/` ve `test/dco/` klasörlerindedir. Betikler proje kökünden sanal ortamın yorumlayıcısıyla çalıştırılır ve ürettikleri çıktılar `test/*/cikti/` altına yazılır. Bu klasörler ile `sys/` veri deposu sürüm denetimine dâhil edilmemiştir.

```bash
.venv/bin/python test/dco/dco_6.py
```

`test/dco/` betikleri diferansiyel yöntemin hastalık gen grupları üzerindeki uygulamasını ve Bölüm 8'de aktarılan rastgele bölünme deneyini içerir. Permütasyon sonuçları `seed` parametresine göre yinelenebilirdir. Sonuçların yinelenmesi için kaynak sürümleri (GENCODE v50, HPO v2026-09-01 ve NCBI kayıt sürümleri) birlikte raporlanmalıdır, çünkü RefSeq kayıtları istek anında geçerli sürümle indirilir.

`doc/` belgelerindeki kod blokları çalıştırılarak denetlenmiş, gösterilen çıktılar gerçek çıktılarla karşılaştırılmıştır.

## 12. Sınırlılıklar

ORF eşleme yalnızca transkriptin kendi zincirinde, AUG ile başlayan ve proteinle tam örtüşen çerçevelerde yapılır. Selenosistein ve pirolizin üretilmediği için selenoproteinlerin ana transkriptleri proteinleriyle eşleşmez: `gc.canonical("GPX1").isoforms()` boş liste döndürür. Okunan stop kodonları ve AUG dışı başlangıçlar da aynı nedenle eşleşmez. Genetik kod olarak yalnızca standart (NCBI 1) ve omurgalı mitokondriyal (NCBI 2) kodlar tanımlanmıştır; yeni bir kod `GeneticCode` sınıfına 64 harflik NCBI dizisiyle eklenebilir. Dairesel diziler desteklenmez, dolayısıyla başlangıç noktasının üzerinden geçen ORF'ler bulunamaz.

tAI hesabında gen kopya sayısı tRNA bolluğunun vekili olarak kabul edilir ve wobble cezaları *Saccharomyces cerevisiae* için optimize edilmiş değerlerdir. Bu nedenle tAI mutlak bir değer olarak değil, aynı koşullarda hesaplanmış genler arasında karşılaştırma amacıyla kullanılmalıdır. CDC uygulamasında Zhang ve ark. (2012) tarafından önerilen bootstrap anlamlılık testi yer almamaktadır. CPB değerinin yorumlanabilmesi için referansın genom ölçeğinde olması gerekir, aksi durumda referansta bulunmayan kodon çiftleri hesap dışında kalır.

`plot_volcano` Fisher sonuçlarını gösterir, `robust` işaretini içermez. HPO ilişkilerinde etiket güvenilirliği sınırlıdır ve `min_frequency` süzgeci sıklığı bilinmeyen satırları (v2026-09-01 sürümünde fenotip satırlarının %14,7'si) da dışlar.

Bazı ölçek değerleri ile Wilkinson-Harrison sabitleri özgün kaynaklarıyla doğrulanmayı beklemektedir; ilgili değerler [doc/tables.md](doc/tables.md#8-kaynağıyla-doğrulanması-gereken-değerler) belgesinde listelenmiştir. Protein-Sol gibi skorların model ağırlıkları pakete dâhil edilmemiştir, `linear_score` metoduna dışarıdan verilir. Sınırlılıkların tamamı ve sık karşılaşılan hata iletileri [doc/sinirliliklar.md](doc/sinirliliklar.md) belgesinde toplanmıştır.

## 13. Belgeler

Başvuru belgeleri `doc/` klasöründedir. `crn/` altındaki her modülün bu klasörde aynı yolda bir karşılığı bulunur ve her sınıf, fonksiyon, metot ve sabit ilgili dosyada açıklanmıştır.

| Belge | İçerik |
|---|---|
| [doc/README.md](doc/README.md) | Dizin ve modül haritası |
| [doc/seq.md](doc/seq.md), [doc/gene.md](doc/gene.md), [doc/fileio.md](doc/fileio.md) | Veri modeli ve dosya biçimleri |
| [doc/fetch/README.md](doc/fetch/README.md) | `Ncbi`, `Gencode`, `Hgnc`, `Hpo`, `NcbiSession` |
| [doc/cub/README.md](doc/cub/README.md) | Kodon kullanımı ölçümleri ve grafikler |
| [doc/cub/differential.md](doc/cub/differential.md) | Diferansiyel analizin başvuru belgesi |
| [diferansiyel_kodon_kullanimi.md](diferansiyel_kodon_kullanimi.md) | Diferansiyel analizin yöntemsel gerekçesi |
| [doc/protein.md](doc/protein.md) | `ProteinProperties` |
| [doc/senaryolar.md](doc/senaryolar.md) | Uçtan uca senaryolar |
| [doc/sinirliliklar.md](doc/sinirliliklar.md) | Sınırlılıklar ve hata iletileri |

## 14. Kaynakça

Bennetzen, J. L. ve Hall, B. D. (1982). Codon selection in yeast. *Journal of Biological Chemistry*, 257(6), 3026-3031. PMID: [7037777](https://pubmed.ncbi.nlm.nih.gov/7037777/)

Benjamini, Y. ve Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289-300.

Coleman, J. R., Papamichail, D., Skiena, S., Futcher, B., Wimmer, E. ve Mueller, S. (2008). Virus attenuation by genome-scale changes in codon pair bias. *Science*, 320(5884), 1784-1787. [doi:10.1126/science.1155761](https://doi.org/10.1126/science.1155761)

dos Reis, M., Savva, R. ve Wernisch, L. (2004). Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Research*, 32(17), 5036-5044. [doi:10.1093/nar/gkh834](https://doi.org/10.1093/nar/gkh834)

Ikemura, T. (1981). Correlation between the abundance of *Escherichia coli* transfer RNAs and the occurrence of the respective codons in its protein genes. *Journal of Molecular Biology*, 146(1), 1-21. [doi:10.1016/0022-2836(81)90363-6](https://doi.org/10.1016/0022-2836(81)90363-6)

Karlin, S. ve Burge, C. (1995). Dinucleotide relative abundance extremes: a genomic signature. *Trends in Genetics*, 11(7), 283-290. [doi:10.1016/s0168-9525(00)89076-9](https://doi.org/10.1016/s0168-9525(00)89076-9)

Mueller, S., Papamichail, D., Coleman, J. R., Skiena, S. ve Wimmer, E. (2006). Reduction of the rate of poliovirus protein synthesis through large-scale codon deoptimization causes attenuation of viral virulence by lowering specific infectivity. *Journal of Virology*, 80(19), 9687-9696. [doi:10.1128/JVI.00738-06](https://doi.org/10.1128/JVI.00738-06)

Roymondal, U., Das, S. ve Sahoo, S. (2009). Predicting gene expression level from relative codon usage bias: an application to *Escherichia coli* genome. *DNA Research*, 16(1), 13-30. [doi:10.1093/dnares/dsn029](https://doi.org/10.1093/dnares/dsn029)

Sharp, P. M. ve Li, W. H. (1987). The codon Adaptation Index: a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Research*, 15(3), 1281-1295. [doi:10.1093/nar/15.3.1281](https://doi.org/10.1093/nar/15.3.1281)

Sueoka, N. (1988). Directional mutation pressure and neutral molecular evolution. *Proceedings of the National Academy of Sciences of the United States of America*, 85(8), 2653-2657. [doi:10.1073/pnas.85.8.2653](https://doi.org/10.1073/pnas.85.8.2653)

Sueoka, N. (1995). Intrastrand parity rules of DNA base composition and usage biases of synonymous codons. *Journal of Molecular Evolution*, 40(3), 318-325. [doi:10.1007/BF00163236](https://doi.org/10.1007/BF00163236)

Wan, X. F., Xu, D., Kleinhofs, A. ve Zhou, J. (2004). Quantitative relationship between synonymous codon usage bias and GC composition across unicellular genomes. *BMC Evolutionary Biology*, 4, 19. [doi:10.1186/1471-2148-4-19](https://doi.org/10.1186/1471-2148-4-19)

Wright, F. (1990). The 'effective number of codons' used in a gene. *Gene*, 87(1), 23-29. [doi:10.1016/0378-1119(90)90491-9](https://doi.org/10.1016/0378-1119(90)90491-9)

Zhang, Z., Li, J., Cui, P., Ding, F., Li, A., Townsend, J. P. ve Yu, J. (2012). Codon Deviation Coefficient: a novel measure for estimating codon usage bias and its statistical significance. *BMC Bioinformatics*, 13, 43. [doi:10.1186/1471-2105-13-43](https://doi.org/10.1186/1471-2105-13-43)

Kullanılan veri kaynakları: [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/refseq/), [GENCODE](https://www.gencodegenes.org/), [HGNC](https://www.genenames.org/), [Human Phenotype Ontology](https://hpo.jax.org/), [GtRNAdb](https://gtrnadb.ucsc.edu/).

## Ek A. Önceki arayüzden geçiş

Çerçevenin yeniden yapılandırılması sırasında değiştirilen arayüzler aşağıda verilmiştir.

| Önceki | Güncel |
|---|---|
| `from crn.io import DNA, ...` | `from crn import DNA, ...` |
| `crn.io.crn_io_const` (alfabeler, `GeneticCode`) | `crn.constants` |
| `crn.io.crn_io_const` (kütleler, ölçekler, pKa, tAI) | `crn.tables` |
| `Maps.*` tabloları, `maps=` parametreleri | `crn.constants.DNA_COMPLEMENT`, `RNA_COMPLEMENT`, `DNA_TO_RNA`, `RNA_TO_DNA`; parametre kaldırılmıştır |
| `seq._type` | `seq.molecule` |
| `d.complementary().seq[::-1]` | `d.reverse_complement()` |
| `gene.match_transcripts_to_proteins()` | `gene.isoforms()` |
| `CUB(genes=[...])` (genler yerinde kırpılır) | `CUB([...])` (genler değiştirilmez; `cub.genes` kopyaları verir) |
| `CUB(genes=a.genes + b.genes)` | `a + b` |
| `ProteinMass`, `ProteinAtoms`, `ProteinCharge`, `ProteinHydropathy`, `ProteinOptics`, `ProteinStability`, `ProteinDisorder` | `ProteinProperties` (aynı metot adları) |
| `ProteinDisorder.features()` | `ProteinProperties.composition_features()` |
| Ölçüm metotlarındaki `codons=` parametresi | Kaldırılmıştır (`codon_counts`, `nucleotide_counts`, `pr2` dışında) |
| `RNA.match_orf` sonucunun adı `<ad>_ORF<n>` | RNA'nın kendi adı |
| `crn.fetch.crn_fetch_ncbi` vb. | `crn.fetch.ncbi`, `crn.fetch.gencode`, `crn.fetch.session` |
