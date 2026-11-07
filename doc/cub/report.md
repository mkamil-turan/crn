# `crn.cub.report`: bütün ölçümler tek tabloda

Kaynak: [`crn/cub/report.py`](../../crn/cub/report.py)

Her gen için bütün skaler ölçümleri tek satıra toplar. Ölçüm sınıflarını tek tek kurmak yerine çoğu analizde başlangıç noktası budur; `CUB.table` aynı fonksiyonu çağırır. Çıktı TSV olarak yazılabilir ya da korelasyon analizine verilebilir.

<!-- toc -->
- [1. Fonksiyonlar](#1-fonksiyonlar): [`summary_table`](#summary_table) · [`correlation_matrix`](#correlation_matrix) · [`plot_correlations`](#plot_correlations)
- [2. İç yardımcılar](#2-iç-yardımcılar)
<!-- /toc -->

## 1. Fonksiyonlar

### summary_table

`summary_table(cub, reference=None, background=None, host=None, genome=None, trna_counts=None, kingdom='eukaryote', protein=True, file_name=None)`

Her gen için bir satır verir: `{"gene": ad, ölçüm: değer, ...}`. Her zaman bulunan sütunlar: L, GC%, GC1, GC2, GC12, GC3, GC3s, ENC, ENC_expected, ENC_ratio, SCUO, CDC, RCBS. İsteğe bağlı sütunlar ilgili veri verilince eklenir.

`cub` ölçülecek settir. `reference` (yüksek ifadeli genler) Fop, CAI, CBI ekler; `background` da verilirse Fop_chi2 ve CBI_chi2 gelir. `host` RCDI, `genome` CPB, `trna_counts` tAI ekler, `kingdom` tAI alemidir. `protein=True` `ProteinProperties.summary` sütunlarını ekler. `file_name` verilirse tablo TSV olarak da yazılır. Satır satır `CubWarning` uyarıları bastırılır.

```python
from crn.fetch import Gencode
from crn.fetch._store import DEFAULT_SYS_DIR
from crn.cub import CUB, summary_table, correlation_matrix, plot_correlations, read_trna_counts

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "HBB", "ACTB", "MEFV", "GAPDH"]))
yuksek = CUB(kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2"]))
dusuk = CUB(kur(["OR1A1", "OR2A2", "OR5A1", "OR10A2"]))
trna = read_trna_counts(DEFAULT_SYS_DIR / "gtrnadb" / "hg38-tRNAs.fa")

satirlar = summary_table(cub, reference=yuksek, background=dusuk, host=yuksek, genome=yuksek,
                         trna_counts=trna, kingdom="eukaryote", protein=False, file_name="ozet.tsv")
print(list(satirlar[0])[:4], "tAI" in satirlar[0], "CBI_chi2" in satirlar[0])  # -> ['gene', 'L', 'GC%', 'GC1'] True True
```

### correlation_matrix

`correlation_matrix(rows, columns=None)`

`summary_table` satırlarından (`rows`) seçili sütunların Pearson korelasyon matrisini verir. `columns` verilmezse `DEFAULT_CORRELATION_COLUMNS` içinden satırlarda bulunan, sayısal ve değişkenliği olan sütunlar alınır. `nan` içeren değerler o çiftte atlanır. `(columns, matrix)` döner.

```python
sutunlar, m = correlation_matrix(satirlar, columns=["GC3s", "ENC"])
print(sutunlar, round(m[0][1], 3))  # -> ['GC3s', 'ENC'] -0.638
```

### plot_correlations

`plot_correlations(cub, columns=None, ax=None, annotate=True, file_name=None, **summary_kwargs)`

CUB indeksleri arasındaki Pearson korelasyonlarının ısı haritasını çizer. `cub` ölçülecek settir. `columns` sütunları elle seçer. `annotate=True` hücrelere değerleri yazar. `**summary_kwargs` `summary_table` parametreleridir: hangi sütunların hesaplanacağını bunlar belirler, ör. `reference=yuksek` CAI'yi ekler. `ax` ve `file_name` [ortak kurala](README.md#1-ortak-kurallar) uyar.

```python
ax = plot_correlations(cub, columns=["GC3s", "ENC", "SCUO"], annotate=True, reference=yuksek, file_name="kor.png")
```

## 2. İç yardımcılar

| Ad | Görevi |
|---|---|
| `_pearson` | İki dizinin Pearson korelasyonu; değişkenlik yoksa `nan` |
