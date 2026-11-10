# crn belgeleri

`crn`, moleküler biyolojinin santral dogmasını (DNA → RNA → protein) Python nesneleriyle modelleyen, gerçek veriyi NCBI, GENCODE, HGNC, HPO ve ClinVar'dan alan, kodon kullanımını ve protein özelliklerini ölçen bir kütüphanedir. Bu klasör paketin başvuru belgesidir: `crn/` altındaki her dosyanın burada aynı yolda bir karşılığı vardır, her sınıf, fonksiyon, metot ve sabit bu dosyalarda anlatılır.

## 1. Nereden başlamalı

| Amacınız | Sayfa |
|---|---|
| Bir diziyi yüklemek, çevirmek, ORF bulmak | [seq.md](seq.md) |
| Bir geni transkript ve proteinleriyle tutmak | [gene.md](gene.md) |
| Bir geni NCBI ya da GENCODE'dan indirmek | [fetch/README.md](fetch/README.md) |
| Terminalden çalışmak | [cli.md](cli.md) |
| Semptomdan gen listesi çıkarmak | [fetch/hpo.md](fetch/hpo.md) |
| Klinik varyantlardan gen ve hastalık listesi çıkarmak | [fetch/clinvar.md](fetch/clinvar.md) |
| Kodon kullanımını ölçmek | [cub/README.md](cub/README.md) |
| İki gen grubunu karşılaştırmak | [cub/differential.md](cub/differential.md) |
| Protein özelliklerini hesaplamak | [protein.md](protein.md) |
| Modülleri birlikte kullanmak | [senaryolar.md](senaryolar.md) |
| Bir sonucu raporlamadan önce sınırları görmek | [sinirliliklar.md](sinirliliklar.md) |

## 2. Modül haritası

| Kaynak | Belge | İçerik |
|---|---|---|
| `crn/constants.py` | [constants.md](constants.md) | Alfabeler, kodonlar, tümleme tabloları, `GeneticCode` |
| `crn/tables.py` | [tables.md](tables.md) | Kütleler, pKa, hidropati ölçekleri, tAI cezaları |
| `crn/features.py` | [features.md](features.md) | `Features` |
| `crn/seq.py` | [seq.md](seq.md) | `Sequence`, `NucleicAcid`, `DNA`, `RNA`, `PRT`, `Molecule`, `Strand`, `GenomicStrand` |
| `crn/gene.py` | [gene.md](gene.md) | `Gene` |
| `crn/fileio.py` | [fileio.md](fileio.md) | FASTA ve CRN dosyaları |
| `crn/protein.py` | [protein.md](protein.md) | `ProteinProperties` |
| `crn/cli.py` | [cli.md](cli.md) | `crn` komutu |
| `crn/fetch/__init__.py` | [fetch/README.md](fetch/README.md) | Veri kaynaklarının karşılaştırması |
| `crn/fetch/_store.py` | [fetch/_store.md](fetch/_store.md) | `LocalSource`, ortak dosya yardımcıları |
| `crn/fetch/session.py` | [fetch/session.md](fetch/session.md) | `NcbiSession` |
| `crn/fetch/ncbi.py` | [fetch/ncbi.md](fetch/ncbi.md) | `Ncbi` |
| `crn/fetch/gencode.py` | [fetch/gencode.md](fetch/gencode.md) | `Gencode`, `AmbiguousGeneError` |
| `crn/fetch/hgnc.py` | [fetch/hgnc.md](fetch/hgnc.md) | `Hgnc` |
| `crn/fetch/_geneset.py` | [fetch/_geneset.md](fetch/_geneset.md) | `GeneSetSource`: `Hpo` ve `Clinvar` için ortak gen listesi işlemleri |
| `crn/fetch/hpo.py` | [fetch/hpo.md](fetch/hpo.md) | `Hpo` |
| `crn/fetch/clinvar.py` | [fetch/clinvar.md](fetch/clinvar.md) | `Clinvar` |
| `crn/cub/__init__.py` | [cub/README.md](cub/README.md) | CUB bölümleri ve ortak kurallar |
| `crn/cub/core.py` | [cub/core.md](cub/core.md) | A: `CUB` |
| `crn/cub/composition.py` | [cub/composition.md](cub/composition.md) | B: `Composition` |
| `crn/cub/synonymous.py` | [cub/synonymous.md](cub/synonymous.md) | C: `Synonymous` |
| `crn/cub/adaptation.py` | [cub/adaptation.md](cub/adaptation.md) | D: `Adaptation` |
| `crn/cub/enc.py` | [cub/enc.md](cub/enc.md) | F: `EffectiveCodons` |
| `crn/cub/divergence.py` | [cub/divergence.md](cub/divergence.md) | G: `Divergence` |
| `crn/cub/context.py` | [cub/context.md](cub/context.md) | H: `CodonContext` |
| `crn/cub/ordination.py` | [cub/ordination.md](cub/ordination.md) | `CorrespondenceAnalysis` |
| `crn/cub/differential.py` | [cub/differential.md](cub/differential.md) | J: `DifferentialCodonUsage` |
| `crn/cub/report.py` | [cub/report.md](cub/report.md) | `summary_table`, `correlation_matrix`, `plot_correlations` |
| `crn/cub/_plotting.py` | [cub/_plotting.md](cub/_plotting.md) | Ortak grafik yardımcıları |

## 3. Kurulum

Python 3.12 ya da üstü gerekir. Bağımlılıklar `pyproject.toml` içindedir: tqdm, requests, matplotlib, scipy, numpy.

```text
python -m venv .venv
.venv/bin/pip install -e .
```

Veriler ilk kullanımda proje kökündeki `sys/` klasörüne indirilir. Başka bir disk kullanmak için `CRN_SYS_DIR` ortam değişkenini ayarlayın. Bu klasör git'e girmez.

## 4. Paketin dışa verdikleri

`from crn import ...` ile şunlar alınır (`__all__`): `GeneticCode`, `Features`, `Molecule`, `Strand`, `GenomicStrand`, `Sequence`, `NucleicAcid`, `DNA`, `RNA`, `PRT`, `Gene`, `read_fasta`, `write_fasta`, `read_crn_list`, `write_crn_list`, `ProteinProperties`. Veri kaynakları `crn.fetch`, kodon ölçümleri `crn.cub` altındadır.

```python
from crn import DNA, GeneticCode
from crn.fetch import Gencode
from crn.cub import CUB, EffectiveCodons

gen = Gencode(verbose=False).canonical("HBB").isoforms()[0]
print(gen.name)  # -> HBB_ENST00000335295.4
print(round(EffectiveCodons(CUB([gen])).enc(gen.name), 2))  # -> 34.01

dna = DNA("ATGGCCTAA", name="kisa")
print(dna.transcription(template=False).translation(GeneticCode.STANDARD).seq)  # -> MA*
```

## 5. Sayfaların okunuşu

Her sayfa kaynak dosyanın adını taşır. `###` başlıkları sınıf ya da metot adıdır, altındaki satır parametreleri kaynak koddan alınmış imzayı gösterir. Örneklerdeki `# ->` işareti satırın gerçek çıktısıdır: `sys/` verileriyle çalıştırılarak yazılmıştır. `# ->` taşımayan ağ örnekleri bu otomatik denetimin dışındadır; `update` çağrıları indirme başlattığı için hiç çalıştırılmamıştır. Sayfa sonlarındaki "İç yardımcılar" tabloları doğrudan çağrılması düşünülmeyen ad ve fonksiyonları listeler; bunlar ileride değişebilir.
