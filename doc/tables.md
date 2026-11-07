# `crn.tables`: sayısal tablolar

Kaynak: [`crn/tables.py`](../crn/tables.py)

Protein ve tRNA ölçümlerinin sabitleri burada toplanır. Bu tablolar fonksiyon içine gömülü değildir: `ProteinProperties` ve `CodonContext` metotları hepsini parametre olarak alır, yani başka bir değer kümesiyle hesap yapmak için kaynak kodu değiştirmek gerekmez.

<!-- toc -->
- [1. Amino asit grupları](#1-amino-asit-grupları)
- [2. Kütle ve formül](#2-kütle-ve-formül)
- [3. Yük ve pKa](#3-yük-ve-pka)
- [4. Hidropati ve enerji ölçekleri](#4-hidropati-ve-enerji-ölçekleri)
- [5. Optik ve hacim](#5-optik-ve-hacim)
- [6. Kararlılık ve düzensizlik](#6-kararlılık-ve-düzensizlik)
- [7. tAI sabitleri](#7-tai-sabitleri)
- [8. Kaynağıyla doğrulanması gereken değerler](#8-kaynağıyla-doğrulanması-gereken-değerler)
<!-- /toc -->

## 1. Amino asit grupları

`AMINO_ACID_GROUPS` grup adından harflere giden sözlüktür. Alifatik, aromatik, küçük, bazik ve asidik gruplar yaygın tanımlardır; hidrofobik ve polar gruplar EMBOSS pepstats'ın "Non-polar" ve "Polar" sınıflarıdır. `ProteinProperties.group_composition` varsayılan olarak bu sözlüğü kullanır.

## 2. Kütle ve formül

| Ad | İçerik | Kaynak |
|---|---|---|
| `RESIDUE_AVERAGE_MASS` | Kalıntı (amino asit − H2O) ortalama kütleleri, Da | ExPASy Compute pI/Mw |
| `RESIDUE_MONOISOTOPIC_MASS` | Kalıntı monoizotopik kütleleri, Da | ExPASy PeptideMass |
| `WATER_AVERAGE_MASS` | Zincir uçlarındaki H ve OH, 18,01524 Da | |
| `WATER_MONOISOTOPIC_MASS` | Aynısı, 18,01056 Da | |
| `ELEMENTS` | Formül demetlerinin sırası: `C H N O S` | |
| `AMINO_ACID_FORMULAS` | Serbest amino asit formülleri, `ELEMENTS` sırasıyla | |
| `WATER_FORMULA` | Peptit bağında ayrılan su, `(0, 2, 0, 1, 0)` | |

## 3. Yük ve pKa

`POSITIVE_RESIDUES` (`KR`) ve `NEGATIVE_RESIDUES` (`DE`) pH 7'de yüklü sayılan kalıntılardır; His ihmal edilir ve kaba yük ölçümleri (FCR, NCPR) bunları kullanır. `IONIZABLE_POSITIVE` (`HKR`) ve `IONIZABLE_NEGATIVE` (`CDEY`) ise Henderson-Hasselbalch hesabında iyonlaşan yan zincirlerdir; uçlar ayrıca `N_term` ve `C_term` anahtarlarıyla eklenir.

`PKA_SETS` üç küme içerir: `EMBOSS`, `Lehninger` ve `Bjellqvist`. Her kümede `N_term`, `C_term`, `K R H D E C Y` değerleri bulunur. Bjellqvist kümesinde `N_term_residue` ve `C_term_residue` sözlükleri de doludur: pKa uçtaki kalıntıya göre değişir, ExPASy Compute pI/Mw böyle hesaplar. Bu yüzden Bjellqvist diğer iki kümeden farklı pI verir ve `crn` onu varsayılan alır.

## 4. Hidropati ve enerji ölçekleri

`HYDROPATHY_SCALES` ExPASy ProtScale değerlerini tutar. Anahtarlar Kyte-Doolittle, Hopp-Woods, Eisenberg, Janin ve Wimley-White ölçekleridir; tam adları `ProteinProperties.hydropathy` çağrısındaki hata mesajında ve aşağıdaki örnekte görülür. Kyte-Doolittle, Eisenberg ve Janin'de büyük değer hidrofobiktir, Hopp-Woods'ta büyük değer hidrofiliktir. Wimley-White bir arayüz ölçeğidir (su → arayüz ΔG, kcal/mol) ve negatif değer arayüzü tercih eder; His nötr değeriyle (0,17) alınmıştır, protonlanmış His için değer 0,96'dır.

`BOMAN_SCALE` yan zincirlerin sikloheksandan suya transfer serbest enerjileridir (kcal/mol; Boman 2003, Radzicka-Wolfenden). Hidrofilik kalıntılar pozitiftir. Pro için kaynakta değer yoktur ve 0 alınmıştır.

```python
from crn.tables import HYDROPATHY_SCALES, PKA_SETS

print(list(HYDROPATHY_SCALES))  # -> ['Kyte-Doolittle', 'Hopp-Woods', 'Eisenberg', 'Wimley-White', 'Janin']
print(list(PKA_SETS))  # -> ['EMBOSS', 'Lehninger', 'Bjellqvist']
```

## 5. Optik ve hacim

| Ad | İçerik | Kaynak |
|---|---|---|
| `EXTINCTION_COEFFICIENTS_280` | 280 nm molar ekstinksiyon (W 5500, Y 1490, disülfit başına 125) | Pace ve ark. 1995 |
| `RESIDUE_SPECIFIC_VOLUMES` | Kısmi spesifik hacimler, cm³/g | Cohn ve Edsall 1943 |
| `RESIDUE_VOLUMES` | Kalıntı hacimleri, Å³ | Zamyatnin 1972 |
| `RESIDUE_MAX_ASA` | Gly-X-Gly'de en büyük erişilebilir yüzey, Å² | Tien ve ark. 2013 |

## 6. Kararlılık ve düzensizlik

| Ad | İçerik | Kaynak |
|---|---|---|
| `IVYWREL_RESIDUES` | Termofilik proteinlerde zenginleşen kalıntılar | Zeldovich ve ark. 2007 |
| `CVP_CHARGED`, `CVP_POLAR` | CvP yanlılığının iki kalıntı kümesi (`DEKR`, `NQST`) | |
| `TURN_FORMING_RESIDUES` | Wilkinson-Harrison modelinde dönüş oluşturan kalıntılar (`NGPS`) | |
| `TOP_IDP_SCALE` | Düzensizlik eğilimi; pozitif değer düzensizliği destekler | Campen ve ark. 2008 |

## 7. tAI sabitleri

`TAI_S_VALUES` wobble cezalarıdır (dos Reis ve ark. 2004, *S. cerevisiae* için optimize). Anahtar `antikodon 34. baz:kodon 3. baz` biçimindedir; `I` A34 deaminasyonundan gelen inozin, `L` prokaryotlarda tRNA-Ile CAU'nun AUA okuyan lizidinli biçimidir. `TAI_WOBBLE_PAIRS` kodonun 3. bazını okuyabilen antikodon bazlarını ve bu eşleşmenin ceza anahtarını verir. Genomdaki A34 olgun tRNA'da I'ya dönüşür, bu yüzden tabloda genomda yazıldığı gibi `A` olarak durur.

## 8. Kaynağıyla doğrulanması gereken değerler

Tabloların tamamı aynı güvende değildir. Kyte-Doolittle değerleri GRAVY ile, ekstinksiyon katsayıları ExPASy ile doğrulanmıştır. Hopp-Woods, Eisenberg, Wimley-White, Janin, Boman, Cohn-Edsall, Zamyatnin, Tien ASA ve TOP-IDP ölçekleri ise bilinen bir sonuçla karşılaştırılamamıştır; yayın öncesinde özgün kaynaklarıyla eşlenmelidir.

Wilkinson-Harrison sabitleri (`crn.protein` içindeki `WH_*`, Davis ve ark. 1999) özgün makaleden doğrulanmamıştır. Ubikitin *E. coli*'de iyi çözündüğü halde model onu çözünmez çıkarır (CV = 1,96). Bunun nedeni modelin kabalığı da sabitlerden birinin yanlış olması da olabilir. Protein-Sol ağırlıkları ise koda hiç konmamıştır: `ProteinProperties.composition_features` özellikleri hesaplar, katsayılar Hebditch ve ark. 2017'nin ek dosyalarından alınıp `ProteinProperties.linear_score`'a verilmelidir.
