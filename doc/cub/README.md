# `crn.cub`: kodon kullanım yanlılığı

Kaynak: [`crn/cub/__init__.py`](../../crn/cub/__init__.py)

Kodon kullanım yanlılığı (CUB), aynı amino asidi kodlayan eş anlamlı kodonların eşit sıklıkta kullanılmamasıdır. Bu paket literatürdeki temel ölçümleri hesaplar. Bütün ölçüm sınıfları bir `CUB` nesnesi alır: `CUB` genleri hazırlar ve kodonları gen başına bir kez sayar, ölçümler bu sayımlardan hesaplanır. Bu yüzden 20 000 genlik bir analizde de yavaşlama olmaz.

| Bölüm | Sınıf | Ölçümler | Belge |
|---|---|---|---|
| A | `CUB` | Kodon ve amino asit sayıları, frekanslar | [core.md](core.md) |
| B | `Composition` | GC%, GC1-3, GC12, GC3s, A3s/U3s/G3s/C3s | [composition.md](composition.md) |
| C | `Synonymous` | Kodon oranı, RSCU, tercih edilen kodonlar, Fop | [synonymous.md](synonymous.md) |
| D | `Adaptation` | w, CAI, CBI | [adaptation.md](adaptation.md) |
| E | `ProteinProperties` | Kütle, formül, yük, hidropati, optik, kararlılık | [../protein.md](../protein.md) |
| F | `EffectiveCodons` | ENC, beklenen ENC, nötrallik, PR2 | [enc.md](enc.md) |
| G | `Divergence` | SCUO, CDC, RCDI, RCBS | [divergence.md](divergence.md) |
| H | `CodonContext` | tAI, dinükleotit ρ, CPS, CPB | [context.md](context.md) |
| J | `DifferentialCodonUsage` | İki gen grubunun kodon kodon karşılaştırılması | [differential.md](differential.md) |
| | `CorrespondenceAnalysis` | RSCU tablosunun uyum çözümlemesi | [ordination.md](ordination.md) |
| | `summary_table`, `correlation_matrix`, `plot_correlations` | Bütün ölçümler tek tabloda | [report.md](report.md) |
| | Grafik yardımcıları | Ortak biçimler | [_plotting.md](_plotting.md) |

## 1. Ortak kurallar

Ölçüm metotlarında `gene` parametresi aynı anlama gelir: verilmezse bütün genler birlikte (havuzlanmış sayımlarla), verilirse yalnızca o gen hesaplanır. Değer gen adıdır, yani `cub.names` içindeki bir öğe; `isoforms()`'tan gelen genlerde `TP53_ENST00000269305.9` biçimindedir. Bilinmeyen ad `KeyError` verir.

Grafik metotlarının hepsi aynı arayüzü izler: `ax` verilirse onun üzerine, verilmezse yeni şekle çizer, `file_name` verilirse 300 dpi kaydeder ve `Axes` döndürür. Gen başına nokta çizenlerde `labels=True` gen sembolünü yazar, `groups` genleri renklendirir (ayrıntı: [_plotting.md](_plotting.md)). matplotlib yalnızca çizim istendiğinde yüklenir.

Referans gen seti isteyen ölçümlerde (Fop, CAI, CBI, RCDI, CPB) varsayılan `cub`'ın kendisidir. Bu çoğu zaman anlamsızdır: tek gen kendisine göre ölçülmüş olur. Referansı açıkça verin.

## 2. En kısa kullanım

```python
from crn.fetch import Gencode
from crn.cub import CUB, Composition, Synonymous, Adaptation, EffectiveCodons, CubWarning

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

cub = CUB(kur(["TP53", "MEFV", "HBB", "ACTB"]))
yuksek = CUB(kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2"]))

ad = cub.names[0]
print(ad, cub.total_codons(ad))  # -> TP53_ENST00000269305.9 393
print(round(Composition(cub).gc3s(ad), 3))  # -> 0.61
print(round(EffectiveCodons(cub).enc(ad), 2))  # -> 51.23
print(round(Adaptation(cub, reference=yuksek).cai(ad), 3))  # -> 0.677
print(sorted(cub.table(reference=yuksek)[0])[:5])  # -> ['CAI', 'CBI', 'CDC', 'CvP', 'ENC']
```
