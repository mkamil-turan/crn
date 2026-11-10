# `crn.cli`: komut satırı arayüzü

Kaynak: [`crn/cli.py`](../crn/cli.py)

Paket kurulduğunda `crn` komutu `.venv/bin/` altına eklenir (`pyproject.toml`, `[project.scripts]`). Kurulum yapılmadan aynı arayüz `python -m crn.cli` ile çalışır. Arayüz `argparse` ile yazılmıştır ve kütüphanenin üzerindeki ince bir katmandır: hesap yapmaz, girdileri çözer, doğru sınıfı çağırır ve sonucu biçimler. `fetch` ve `cub` paketleri yalnızca kullanan komutta yüklendiği için `crn --help` scipy ve matplotlib beklemez.

<!-- toc -->
- [1. Ortak kurallar](#1-ortak-kurallar)
- [2. Komutlar](#2-komutlar)
- [3. Girdi biçimleri](#3-girdi-biçimleri)
- [4. Örnekler](#4-örnekler)
- [5. Python'dan kullanım](#5-pythondan-kullanım): [`main`](#main) · [`build_parser`](#build_parser)
- [6. Sabitler ve iç yardımcılar](#6-sabitler-ve-iç-yardımcılar)
<!-- /toc -->

## 1. Ortak kurallar

Veri stdout'a ya da `-o` dosyasına, ilerleme, uyarı ve hatalar stderr'e yazılır, böylece komutlar boruyla bağlanabilir. Tablolar sekmeyle ayrılmıştır ve ondalıklar 6 anlamlı basamakla yazılır. Tam değer gerekiyorsa aynı işlem Python'dan yapılmalıdır.

FASTA çıktısı her zaman 5'→3' okunuştadır, çünkü FASTA yön taşımaz ve okunurken 5'→3' kabul edilir. Örneğin `reverse-transcribe` kütüphanede 3'→5' yazılmış cDNA üretir, komut satırında bu dizinin 5'→3' okunuşu görülür. Yön bilgisini korumak için `-o dosya.crn` kullanılmalıdır.

Çıkış kodları: 0 başarı, 1 işlem hatası (kütüphanenin `ValueError`, `KeyError` ya da `OSError` hatası), 2 kullanım hatası, 130 kullanıcı kesmesi. `--debug` tam izleme yığınını gösterir.

İlerleme çubukları yalnızca stderr bir terminale bağlıysa gösterilir, `--quiet` bunları terminalde de kapatır. Uyarılar `--quiet` ile kapanmaz.

`--sys-dir` veri deposunu değiştirir, verilmezse `CRN_SYS_DIR` ya da proje kökündeki `sys/` kullanılır. NCBI API anahtarı `NCBI_API_KEY` ortam değişkeninden okunur.

GENCODE verisi hazır değilse komut 1,3 GB'lık indirmeyi kendiliğinden başlatmaz, hata verir ve `crn update gencode` komutunu önerir. HPO ve HGNC verileri küçük olduğundan ilk kullanımda indirilir.

## 2. Komutlar

| Komut | İşlem |
|---|---|
| `crn seq info` | Ad, tip, uzunluk, GC yüzdesi |
| `crn seq revcomp` | Ters tümleyici |
| `crn seq transcribe` | DNA'dan RNA; `--strand coding` (varsayılan) ya da `template` |
| `crn seq reverse-transcribe` | RNA'dan tümleyici cDNA |
| `crn seq translate` | İlk nükleotitten çeviri; DNA kodlayan zincir sayılır |
| `crn seq orfs` | Altı çerçevede ORF; `--tsv` koordinat tablosu verir |
| `crn protein` | `ProteinProperties.summary` tablosu |
| `crn get dna`, `rna`, `protein`, `region` | Tek kayıt indirme |
| `crn taxonomy AD\|TAXID...` | Organizma adından NCBI Taxonomy ID'sine ya da taxid'den organizmaya |
| `crn gene` | Gen özeti; `--fasta dna\|rna\|cds\|protein`, `--isoforms`, `-o dosya.crn` |
| `crn cub table` | `summary_table`; `--reference`, `--background`, `--host`, `--genome`, `--trna` |
| `crn cub codons` | Kodon sayısı, binde kullanım, RSCU |
| `crn cub plot TÜR` | `enc`, `enc-ratio`, `neutrality`, `pr2`, `rscu`, `rscu-heatmap`, `cai`, `coa`, `correlations` |
| `crn dco A B` | `DifferentialCodonUsage`; `--report robust\|fisher\|significant\|permutation\|levels\|interpret`, `--level`, `--fdr`, `--permutations`, `--plot` |
| `crn hpo search`, `genes`, `symbols`, `phenotypes`, `diseases`, `summary`, `pack` | `Hpo` sorguları ve `--min-frequency`, `--association`, `--evidence` gibi süzgeçler |
| `crn hgnc columns`, `symbols`, `query`, `orthologs` | `Hgnc` sorguları; süzgeçler `sütun=değer`, virgüllü değer öğelerinden biriyle eşleşir |
| `crn clinvar variants`, `genes`, `symbols`, `conditions`, `summary`, `pack` | `Clinvar` sorguları; `--significance`, `--min-stars`, `--effect`, `--condition` gibi süzgeçler |
| `crn update ncbi GEN...`, `gencode`, `hpo`, `hgnc`, `clinvar` | Veri deposunu hazırlar; `--force` yeniden indirir, `clinvar` için `--release` ve `--assembly`, `hgnc` için `--species` (ek HCOP türleri) |

`crn dco`'nun varsayılan raporu `robust`'tur, yani yalnızca Fisher ve permütasyondan birlikte geçen kodonlar yazılır. Gerekçesi için bkz. [differential.md](cub/differential.md#3-neden-permütasyon-gerekir). Her komutun seçenekleri `crn KOMUT --help` ile görülür.

`crn clinvar` süzgeçleri `Clinvar.variants` filtreleriyle aynıdır (bkz. [clinvar.md](fetch/clinvar.md#2-filtreler)). Yinelenebilen seçenekler (`--significance pathogenic --significance likely_pathogenic`) VEYA ile, farklı seçenekler VE ile birleşir. `--max-genes` varsayılan olarak 1'dir ve büyük kopya sayısı değişimlerini sayımlardan çıkarır; `--max-genes 0` sınırı kaldırır. Gen bulunamazsa komut 1 koduyla çıkar.

Organizma `--organism` ya da eşdeğeri `--taxid` ile seçilir ve ad (`mouse`) ya da NCBI Taxonomy ID (`9031`) alır. NCBI kaynağında taxid ile RefSeq'teki herhangi bir organizmanın geni alınır; GENCODE yalnızca insan (9606) ve fareyi (10090) kapsar. Bir türün taxid'i `crn taxonomy` ile bulunur. HCOP ortologlarında tür `taxid=` süzgeciyle seçilir ve hazırlanmamış HCOP türü ilk istekte indirilir.

`crn seq transcribe` varsayılanı kütüphaneden farklıdır: `DNA.transcription` varsayılan olarak girdiyi kalıp zincir sayar, komut ise kodlayan zincir sayar. Terminalde yazılan dizi çoğunlukla gen dizisi olduğundan ve NCBI ile GENCODE gen dizileri kodlayan zincir olduğundan bu seçim yapılmıştır.

## 3. Girdi biçimleri

Dizi girdisi (`GİRDİ`) doğrudan dizi, FASTA dosyası, `.crn` dosyası ya da stdin'dir. Argüman verilmezse stdin okunur. Mevcut dosya adı doğrudan diziye göre önceliklidir. Tip verilmezse yalnızca `ACGTUN-` harflerinden oluşan dizi nükleik asit sayılır. M, K, W gibi IUPAC harfleri aynı zamanda amino asit olduğundan bu tür belirsiz dizilerde `--type` verilmelidir.

Gen girdisi (`GENLER`) şunlardan biridir: gen sembolü, virgülle ayrılmış semboller, `.crn` dosyası, satır başına bir sembol içeren metin dosyası ya da stdin (`-`). Semboller `--source` kaynağından alınır ve izoformlara ayrılır. `--source gencode --canonical` yalnızca ana transkripti kullanır. Bir sembol bulunamazsa ya da proteinle eşleşen izoformu yoksa (örneğin selenoproteinler) o gen uyarıyla atlanır, hiç gen kalmazsa komut hata verir.

## 4. Örnekler

```bash
crn seq revcomp ATGCGC --raw
# GCGCAT

crn seq translate ATGAAATTTTAA --raw
# MKF*

crn seq orfs ATGAAATTTGGGTAACCCATGCCCTAG --tsv
# name        start  end  strand  frame  length
# seq1_ORF1   1      15   +       0      15
# seq1_ORF2   19     27   +       0      9

crn gene MEFV
# gene  kind     id              length
# MEFV  dna      MEFV            14607
# MEFV  rna      NM_000243.3     3506
# MEFV  rna      NM_001198536.2  3048
# MEFV  protein  NP_000234.1     781
# MEFV  protein  NP_001185465.2  445

crn gene MEFV --isoforms -o mefv.crn
crn cub table mefv.crn --reference RPL3,RPL4,RPS3 -o mefv.tsv
crn cub plot enc mefv.crn RPL3,RPL4 -o enc.png
```

Komutlar boruyla bağlanabilir. Aşağıdaki satır bir HGNC gen ailesinin ölçümlerini üretir.

```bash
crn hgnc symbols gene_group=Aquaporins | crn cub table - --no-protein > aquaporin.tsv
```

İki grubun karşılaştırılması doğrudan genlerle ya da `crn hpo pack` dosyasıyla yapılır.

```bash
crn dco RPL3,RPL4,RPS3,RPS6,EEF1A1,EEF2,RPL5,RPS8 \
        OR1A1,OR2A2,OR5A1,OR10A2,OR4D1,OR6A2,OR1E1,OR2B2 \
        --source gencode --canonical --name-a ribozom --name-b koku \
        --permutations 2000 --report interpret --plot volkan.png

crn hpo pack nobet=Seizure ataksi=Ataxia --association MENDELIAN -o gruplar.crn
crn dco nobet ataksi --file gruplar.crn -o robust.tsv
```

Taxid ile çalışma: önce türün taxid'i bulunur, sonra gen bu taxid'le alınır. HCOP ortologunun Gene ID'si doğrudan `crn gene`'e verilebilir.

```bash
crn taxonomy "Gallus gallus" "fission yeast"
# query          taxid  scientific_name            common_name    rank     ...
# Gallus gallus  9031   Gallus gallus              chicken        species
# fission yeast  4896   Schizosaccharomyces pombe  fission yeast  species

crn gene MSTN --taxid 9031
# gene  kind     id              length
# MSTN  dna      MSTN            6687
# MSTN  rna      NM_001001461.2  2322
# MSTN  protein  NP_001001461.1  375

crn hgnc orthologs human_symbol=MSTN taxid=9031,7955 --fields species,ortholog_symbol,ortholog_entrez_gene
crn gene 798441 --taxid 7955 --isoforms -o mstnb.crn        # zebra balığı mstnb, Gene ID ile
```

ClinVar sorguları varyanttan gene ve hastalığa gider. `pack` grupları hastalık adı ya da kimliğiyle tanımlar.

```bash
crn clinvar variants MEFV --significance pathogenic --min-stars 1 | cut -f1,3,8,9
# allele_id  gene_symbol  hgvs_p     effect
# 17577      MEFV         Met694Val  missense
# 17578      MEFV         Met694Ile  missense

crn clinvar genes --condition OMIM:249100 --significance pathogenic
# ncbi_gene_id  gene_symbol  variants  pathogenic  uncertain  benign  conflicting  max_stars
# 4210          MEFV         12        12          0          0       0            2

crn clinvar conditions "mediterranean fever" | cut -f1,3,4
# condition_name                                    variants  genes
# Familial Mediterranean fever                      1187      1
# Familial Mediterranean fever, autosomal dominant  344       3

crn clinvar symbols --condition cardiomyopathy --significance pathogenic \
        --significance likely_pathogenic --min-stars 2 --effect nonsense |
        crn cub table - --source gencode --canonical > kmp.tsv

crn clinvar pack kmp="hypertrophic cardiomyopathy" lynch="Lynch syndrome" \
        --significance pathogenic --min-stars 2 -o clinvar.crn
crn dco kmp lynch --file clinvar.crn --report fisher -o fisher.tsv
```

## 5. Python'dan kullanım

### main

`main(argv=None)`

Komutu çalıştırır ve çıkış kodunu döndürür. `argv` argüman listesidir, `None` ise `sys.argv` kullanılır. `crn` komutunun giriş noktasıdır. Kullanım hatasında `argparse` `SystemExit` yükseltir.

```python
from crn.cli import main

kod = main(["seq", "revcomp", "ATGC", "--raw"])  # stdout: GCAT
print(kod)  # -> 0
```

### build_parser

`build_parser()`

Bütün komutları içeren `argparse.ArgumentParser` nesnesini kurar. Komut ağacını test etmek ya da belge üretmek için kullanılabilir.

```python
from crn.cli import build_parser

args = build_parser().parse_args(["dco", "A", "B", "--fdr", "0.01"])
print(args.report, args.level, args.fdr, args.permutations)  # -> robust family 0.01 5000
```

## 6. Sabitler ve iç yardımcılar

`GENETIC_CODES` komut satırı adlarını `GeneticCode` üyelerine eşler. `SEQUENCE_TYPES` `--type` değerlerini sınıflara eşler. `FASTA_WIDTH` FASTA satır uzunluğudur. `_NUCLEIC_LETTERS` otomatik tip tanımada kullanılan harf kümesi, `_FILE_SUFFIXES` dosya gibi görünen ama bulunamayan girdiyi diziden ayıran uzantılardır. `PLOT_KINDS` ve `DCO_REPORTS` grafik ve rapor seçenekleridir.

Çıktı: `_log` stderr'e yazar, `_emit` metni stdout'a ya da dosyaya yazar, `_cell` tablo hücresini biçimler, `_write_table` TSV, `_write_lines` satır listesi yazar, `_forward` 5'→3' okunuşu verir, `_fasta_text` FASTA metni kurar, `_write_sequences` çıktı biçimini seçer.

Dizi girdisi: `_load_crn` türüne bakmadan CRN okur, `_read_stdin` terminalde bekleme yerine hata verir, `_parse_text` stdin metnini çözer, `_read_sequences` tek girdiyi okur, `_coerce` diziyi sınıfa çevirir, `_sequences` bir komutun bütün girdilerini toplar.

Veri kaynakları: `_verbose` ilerleme çubuğu kararını verir, `_ncbi` ve `_gencode` kaynak nesnelerini kurar, `_source` seçili kaynağı döndürür, `_fetch_gene` geni alır, `_symbols` gen girdisini sembollere çözer, `_load_genes` `CUB`'a hazır liste kurar, `_cub` `CUB` nesnesini, `_reference_sets` `summary_table` referanslarını üretir. `_hpo`, `_hpo_filters`, `_clinvar`, `_clinvar_filters`, `_hgnc`, `_hgnc_filters`, `_fields` HPO ve HGNC argümanlarını çevirir, `_groups` `AD=DEĞER` argümanlarını grup sözlüğüne çevirir, `_pack_report` `pack` raporunu stderr'e yazar, `_gene_sequences` genin istenen dizilerini seçer.

Komutlar: `_cmd_seq_info`, `_cmd_seq_revcomp`, `_cmd_seq_transcribe`, `_cmd_seq_reverse_transcribe`, `_cmd_seq_translate`, `_cmd_seq_orfs`, `_cmd_protein`, `_cmd_get`, `_cmd_taxonomy`, `_cmd_region`, `_cmd_gene`, `_cmd_update`, `_cmd_cub_table`, `_cmd_cub_codons`, `_cmd_cub_plot`, `_cmd_dco`, `_cmd_hpo_search`, `_cmd_hpo_genes`, `_cmd_hpo_symbols`, `_cmd_hpo_phenotypes`, `_cmd_hpo_diseases`, `_cmd_hpo_summary`, `_cmd_hpo_pack`, `_cmd_clinvar_variants`, `_cmd_clinvar_genes`, `_cmd_clinvar_symbols`, `_cmd_clinvar_conditions`, `_cmd_clinvar_summary`, `_cmd_clinvar_pack`, `_cmd_hgnc_columns`, `_cmd_hgnc_symbols`, `_cmd_hgnc_query`, `_cmd_hgnc_orthologs` her biri `args` alır ve sonucu yazar.

Ayrıştırıcı: `_common_parser`, `_output_parser`, `_raw_parser`, `_code_parser`, `_input_parser`, `_source_parser`, `_hpo_parser`, `_clinvar_parser`, `_reference_parser` alt komutların paylaştığı seçenekleri tanımlar. `_LazyChoices` ClinVar seçenek listelerini (`SIGNIFICANCES`, `EFFECTS`, `ASSEMBLIES`) ilk denetimde yükler; böylece `crn --help` ve ClinVar dışındaki komutlar `crn.fetch`'i (requests, tqdm) yüklemez. `_LazyChoices.__init__` sabitin adını alır, `_LazyChoices._values` sabiti modülden okur, `_LazyChoices.__contains__` ve `_LazyChoices.__iter__` `argparse`'ın seçenek denetimi ile yardım metnine listeyi verir. `_show_warning` Python uyarılarını `crn: uyarı:` satırına çevirir.
