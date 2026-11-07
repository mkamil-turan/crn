# Uçtan uca senaryolar

Her senaryo kendi başına çalışır ve modüllerin birlikte nasıl kullanıldığını gösterir. Daha uzun örnekler `test/` altındadır: `test/sequence/`, `test/ncbi/`, `test/gencode/`, `test/hpo/`, `test/dco/`.

<!-- toc -->
- [1. Tek geni incelemek: DNA'dan proteine](#1-tek-geni-incelemek-dnadan-proteine)
- [2. Bir gen grubunun kodon kullanımı](#2-bir-gen-grubunun-kodon-kullanımı)
- [3. Semptomdan gen grubuna](#3-semptomdan-gen-grubuna)
- [4. İki hastalık grubunu karşılaştırmak](#4-iki-hastalık-grubunu-karşılaştırmak)
- [5. Bir gen ailesini ortologlarıyla listelemek](#5-bir-gen-ailesini-ortologlarıyla-listelemek)
- [6. Kayan pencereyle GC profili](#6-kayan-pencereyle-gc-profili)
- [7. Sonucu saklamak](#7-sonucu-saklamak)
<!-- /toc -->

## 1. Tek geni incelemek: DNA'dan proteine

MEFV'nin izoformlarını ayırır, ana izoformun CDS'ini çevirir ve proteinin özelliklerini hesaplar.

```python
import warnings
from crn import ProteinProperties
from crn.fetch import Ncbi

warnings.simplefilter("ignore")
ncbi = Ncbi(verbose=False)
izoformlar = ncbi.gene("MEFV").isoforms()
print([g.name for g in izoformlar][:2])  # -> ['MEFV_NM_000243.3', 'MEFV_NM_001198536.2']

mefv = izoformlar[0]
rna = next(iter(mefv.transcripts.values()))
prt = next(iter(mefv.proteins.values()))
print(rna.features.first("orf_start"), len(prt))  # -> 47 781
print(round(ProteinProperties(prt).molecular_weight(), 1))  # -> 86443.8
```

## 2. Bir gen grubunun kodon kullanımı

Ribozomal genleri, referansa göre ölçer ve sonucu TSV'ye yazar.

```python
from crn.fetch import Gencode
from crn.cub import CUB

gc = Gencode(verbose=False)
def kur(semboller):
    return [gc.canonical(s).isoforms()[0] for s in semboller]

referans = CUB(kur(["RPL3", "RPL4", "RPS3", "RPS6", "EEF1A1", "EEF2"]))
cub = CUB(kur(["TP53", "HBB", "ACTB", "GAPDH"]))
satirlar = cub.table(reference=referans, file_name="ozet.tsv")
en_yuksek = max(satirlar, key=lambda s: s["CAI"])
print(en_yuksek["gene"].split("_")[0], round(en_yuksek["CAI"], 3))  # -> ACTB 0.843
```

## 3. Semptomdan gen grubuna

HPO terimini arar ve gen listesini çıkarır. `cerebellar ataxia` bir eş anlamlıdır, arama onu `Ataxia` terimine götürür.

```python
from crn.fetch import Hpo

hpo = Hpo(verbose=False)
terim = hpo.search("cerebellar ataxia", limit=1)[0]
print(terim["id"], terim["name"])  # -> HP:0001251 Ataxia
semboller = hpo.symbols(terim["id"], min_frequency=0.8)
print(len(semboller) > 10, "ATXN1" in hpo.symbols(terim["id"]))  # -> True True
```

## 4. İki hastalık grubunu karşılaştırmak

İki grubu tek CRN dosyasına yazar, kodon farkını Fisher ve permütasyonla sınar. Örnek genleri elle seçilmiştir ve kontrol grubu eşleştirilmemiştir, bu yüzden sonucun biyolojik bir anlamı yoktur; gerçek analiz için `test/dco/` örneklerine bakın.

```python
from crn.cub import DifferentialCodonUsage

gruplar = {"sca": {"gene": ["ATXN1", "ATXN3", "CACNA1A", "ATXN7", "TBP", "PPP2R2B", "SPTBN2", "ITPR1"]},
           "cmt": {"gene": ["PMP22", "MPZ", "GJB1", "MFN2", "GDAP1", "NEFL", "HSPB1", "SH3TC2"]}}
rapor = hpo.pack(gruplar, "sca_cmt.crn")
print(rapor["groups"], rapor["skipped"])  # -> {'sca': 8, 'cmt': 8} {}

dcu = DifferentialCodonUsage.from_crn("sca_cmt.crn", "sca", "cmt")
print(len(dcu.significant("family")), len(dcu.robust("family", n=2000)))  # -> 0 0
```

Sekizer genlik bu gruplarda iki test de anlamlı kodon bulmaz. Grup büyüdükçe Fisher kolayca anlamlı sonuç üretmeye başlar, [permütasyon bu yüzden gerekir](cub/differential.md#3-neden-permütasyon-gerekir).

## 5. Bir gen ailesini ortologlarıyla listelemek

HGNC'den aileyi alır, GENCODE'dan protein uzunluklarını okur, HCOP'tan fare karşılıklarını bulur.

```python
from crn.fetch import Hgnc

hgnc = Hgnc(verbose=False)
aile = hgnc.symbols(gene_group="Aquaporins")[:3]
for sembol in aile:
    uzunluk = len(next(iter(gc.canonical(sembol).proteins.values())))
    fare = [s["ortholog_symbol"] for s in hgnc.orthologs(human_symbol=sembol, species="mouse")]
    print(sembol, uzunluk, fare[:1])
```

## 6. Kayan pencereyle GC profili

`Sequence.sliding` uzun bir dizide bölgesel değişimi gösterir.

```python
dna = gc.dna("HBB")
gc_profili = [sum(b in "GC" for b in p) / len(p) for p in dna.sliding(200, step_size=200)]
print(len(gc_profili), round(max(gc_profili), 3))  # -> 19 0.55
```

## 7. Sonucu saklamak

İndirme ve izoform arama bir kez yapılır, sonraki çalıştırmalar CRN'den okur.

```python
from crn import write_crn_list, read_crn_list

write_crn_list(kur(["TP53", "HBB"]), "genler.crn")
genler = read_crn_list("genler.crn")
print(CUB(genler))  # -> CUB: 2 gen, L=540
```
