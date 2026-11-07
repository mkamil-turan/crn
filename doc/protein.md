# `crn.protein`: protein özellikleri

Kaynak: [`crn/protein.py`](../crn/protein.py)

`ProteinProperties` bir proteinin dizisinden fizikokimyasal özelliklerini hesaplar. Kaynak bir `CUB` olabilir: o zaman proteinler genlerin CDS'lerinden okunan amino asitlerdir. Bir ya da birden çok `PRT` de olabilir: o zaman dizideki 20 standart amino asit sayılır, `*`, `-` ve belirsiz kodlar atlanır.

Bütün ölçüm metotlarında `gene` aynı anlama gelir: `None` ise bütün proteinler birlikte, verilirse yalnızca o protein hesaplanır (`CUB`'da gen adı, `PRT`'de protein adı). Birlikte hesapta her protein ayrı bir zincirdir: kütleye her biri için bir H2O eklenir, uç gruplar ve disülfit eşleri zincir başına sayılır.

Bölüm numaraları (E-A ... E-G) ve ölçüm numaraları `crn.cub` ile ortak listeyi izler. Tablolar [`crn.tables`](tables.md) içindedir ve parametreyle değiştirilebilir.

<!-- toc -->
- [1. Sınıf](#1-sınıf): [`ProteinProperties`](#proteinproperties) · [`names`](#proteinpropertiesnames) · [`chain_count`](#proteinpropertieschain_count) · [`summary`](#proteinpropertiessummary)
- [2. E-A: sayım ve kütle](#2-e-a-sayım-ve-kütle): [`length`](#proteinpropertieslength) · [`composition`](#proteinpropertiescomposition) · [`group_composition`](#proteinpropertiesgroup_composition) · [`scale_average`](#proteinpropertiesscale_average) · [`molecular_weight`](#proteinpropertiesmolecular_weight) · [`monoisotopic_mass`](#proteinpropertiesmonoisotopic_mass) · [`mean_residue_mass`](#proteinpropertiesmean_residue_mass) · [`residue_counts`](#proteinpropertiesresidue_counts) · [`disulfide_pairs`](#proteinpropertiesdisulfide_pairs)
- [3. E-B: atom bileşimi](#3-e-b-atom-bileşimi): [`formula`](#proteinpropertiesformula) · [`formula_string`](#proteinpropertiesformula_string) · [`atom_count`](#proteinpropertiesatom_count) · [`atoms_per_residue`](#proteinpropertiesatoms_per_residue) · [`carbon_oxidation_state`](#proteinpropertiescarbon_oxidation_state)
- [4. E-C: yük](#4-e-c-yük): [`termini`](#proteinpropertiestermini) · [`charged_counts`](#proteinpropertiescharged_counts) · [`simple_net_charge`](#proteinpropertiessimple_net_charge) · [`fcr`](#proteinpropertiesfcr) · [`ncpr`](#proteinpropertiesncpr) · [`net_charge`](#proteinpropertiesnet_charge) · [`isoelectric_point`](#proteinpropertiesisoelectric_point) · [`titration_curve`](#proteinpropertiestitration_curve)
- [5. E-D: hidropati](#5-e-d-hidropati): [`gravy`](#proteinpropertiesgravy) · [`hydropathy`](#proteinpropertieshydropathy) · [`aromaticity`](#proteinpropertiesaromaticity) · [`aliphatic_index`](#proteinpropertiesaliphatic_index) · [`boman_index`](#proteinpropertiesboman_index)
- [6. E-E: optik ve hacim](#6-e-e-optik-ve-hacim): [`extinction_coefficient`](#proteinpropertiesextinction_coefficient) · [`absorbance_01`](#proteinpropertiesabsorbance_01) · [`partial_specific_volume`](#proteinpropertiespartial_specific_volume) · [`mean_residue_volume`](#proteinpropertiesmean_residue_volume) · [`mean_residue_surface`](#proteinpropertiesmean_residue_surface)
- [7. E-F: kararlılık ve çözünürlük](#7-e-f-kararlılık-ve-çözünürlük): [`ivywrel`](#proteinpropertiesivywrel) · [`cvp_bias`](#proteinpropertiescvp_bias) · [`solubility_cv`](#proteinpropertiessolubility_cv) · [`solubility_probability`](#proteinpropertiessolubility_probability) · [`composition_entropy`](#proteinpropertiescomposition_entropy)
- [8. E-G: düzensizlik](#8-e-g-düzensizlik): [`uversky`](#proteinpropertiesuversky) · [`top_idp`](#proteinpropertiestop_idp) · [`composition_features`](#proteinpropertiescomposition_features) · [`linear_score`](#proteinpropertieslinear_score)
- [9. İç yardımcılar ve sabitler](#9-iç-yardımcılar-ve-sabitler)
<!-- /toc -->

## 1. Sınıf

Aşağıdaki örnekler ubikitin üzerinden yürür (76 kalıntı).

### ProteinProperties

`ProteinProperties(source)`

`source` bir `CUB`, bir `PRT` ya da `PRT` listesidir. Başka tip `TypeError`, boş liste ya da yinelenen ad `ValueError` verir. Kaynak metotlarına `ProteinProperties.source` üzerinden erişilir.

```python
from crn import PRT, ProteinProperties

ubq = PRT("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG", name="UBQ")
pp = ProteinProperties(ubq)
print(pp)  # -> ProteinProperties(1 protein)
```

### ProteinProperties.names

`ProteinProperties.names  (özellik)`

Hesaba giren protein adları.

```python
print(pp.names)  # -> ['UBQ']
```

### ProteinProperties.chain_count

`ProteinProperties.chain_count(gene=None)`

Hesaba giren zincir sayısı: `gene` `None` ise protein sayısı, değilse 1.

```python
print(ProteinProperties([ubq, PRT("MK", name="k")]).chain_count())  # -> 2
```

### ProteinProperties.summary

`ProteinProperties.summary(gene=None)`

Başlıca skaler özellikleri ProtParam benzeri tek sözlükte verir. `crn.cub.summary_table` protein sütunlarını buradan alır.

```python
ozet = pp.summary()
print(list(ozet)[:6])  # -> ['length', 'MW', 'pI', 'net_charge_7', 'GRAVY', 'aromaticity']
```

## 2. E-A: sayım ve kütle

### ProteinProperties.length

`ProteinProperties.length(gene=None, amino_acids=AMINO_ACIDS)`

Madde 1. Uzunluk: amino asit sayılarının toplamı, stop hariç. `amino_acids` sayılacak kalıntıları sınırlar.

```python
print(pp.length(), pp.length(amino_acids="KR"))  # -> 76 11
```

### ProteinProperties.composition

`ProteinProperties.composition(gene=None, amino_acids=AMINO_ACIDS)`

Madde 2. Her amino asidin yüzdesi, n_i / N × 100. `amino_acids` listelenecek kalıntılardır.

```python
print(round(pp.composition()["L"], 2), pp.composition(amino_acids="W"))  # -> 11.84 {'W': 0.0}
```

### ProteinProperties.group_composition

`ProteinProperties.group_composition(gene=None, groups=AMINO_ACID_GROUPS)`

Madde 3. `groups` içindeki her grubun oranı (0-1). Varsayılan `AMINO_ACID_GROUPS`'tur.

```python
print(round(pp.group_composition()["aromatic"], 3))  # -> 0.039
print(pp.group_composition(groups={"sulfur": "CM"}))  # -> {'sulfur': 0.013157894736842105}
```

### ProteinProperties.scale_average

`ProteinProperties.scale_average(scale, gene=None)`

Bir kalıntı ölçeğinin bileşimle ağırlıklı ortalaması: Σ(n_i × ölçek_i) / N. `scale` sözlüğünde bulunan amino asitler sayılır. GRAVY, Boman, kalıntı hacmi ve TOP-IDP bu ortak formülü kullanır; kendi ölçeğinizi de buradan uygularsınız.

```python
print(round(pp.scale_average({"K": 1.0, "R": 1.0}), 3))  # -> 1.0
```

### ProteinProperties.molecular_weight

`ProteinProperties.molecular_weight(gene=None, masses=RESIDUE_AVERAGE_MASS, water=WATER_AVERAGE_MASS)`

Madde 4. Ortalama moleküler ağırlık (Da): Σ(n_i × M_i) + zincir başına H2O. `masses` kalıntı kütleleri, `water` su kütlesidir.

```python
print(round(pp.molecular_weight(), 2))  # -> 8564.84
```

### ProteinProperties.monoisotopic_mass

`ProteinProperties.monoisotopic_mass(gene=None, masses=RESIDUE_MONOISOTOPIC_MASS, water=WATER_MONOISOTOPIC_MASS)`

Madde 5. Monoizotopik kütle (Da). `masses` ve `water` varsayılan olarak monoizotopik tablolardır.

```python
print(round(pp.monoisotopic_mass(), 2))  # -> 8559.62
```

### ProteinProperties.mean_residue_mass

`ProteinProperties.mean_residue_mass(gene=None, masses=RESIDUE_AVERAGE_MASS, water=WATER_AVERAGE_MASS)`

Madde 6. Ortalama kalıntı kütlesi: MW / N. `masses` ve `water` `molecular_weight` ile aynıdır.

```python
print(round(pp.mean_residue_mass(), 2))  # -> 112.7
```

### ProteinProperties.residue_counts

`ProteinProperties.residue_counts(gene=None, residues='CWY')`

Madde 7. `residues` içindeki kalıntıların sayısı; varsayılan Cys, Trp, Tyr.

```python
print(pp.residue_counts(), pp.residue_counts(residues="K"))  # -> {'C': 0, 'W': 0, 'Y': 1} {'K': 7}
```

### ProteinProperties.disulfide_pairs

`ProteinProperties.disulfide_pairs(gene=None)`

Madde 7. Olası disülfit bağı sayısı: zincir başına ⌊nCys / 2⌋. Gerçek bağ sayısını değil, kuramsal üst sınırı verir.

```python
print(ProteinProperties(PRT("CCAC")).disulfide_pairs())  # -> 1
```

## 3. E-B: atom bileşimi

### ProteinProperties.formula

`ProteinProperties.formula(gene=None, formulas=AMINO_ACID_FORMULAS, water=WATER_FORMULA, elements=ELEMENTS)`

Madde 8. Atom formülü: Σ(n_i × amino asit formülü) − (N − zincir sayısı) × H2O. `formulas`, `water` ve `elements` tabloları değiştirir.

```python
print(pp.formula())  # -> {'C': 378, 'H': 629, 'N': 105, 'O': 118, 'S': 1}
```

### ProteinProperties.formula_string

`ProteinProperties.formula_string(gene=None, **kwargs)`

Madde 8. Formül metni; 0 olan elementler yazılmaz. `**kwargs` `formula`'ya gider.

```python
print(pp.formula_string())  # -> C378H629N105O118S1
```

### ProteinProperties.atom_count

`ProteinProperties.atom_count(gene=None, **kwargs)`

Madde 9. Toplam atom sayısı. `**kwargs` `formula`'ya gider.

```python
print(pp.atom_count())  # -> 1231
```

### ProteinProperties.atoms_per_residue

`ProteinProperties.atoms_per_residue(gene=None, **kwargs)`

Madde 10. Kalıntı başına atom sayısı.  `**kwargs` `formula`'ya gider.

```python
print(round(pp.atoms_per_residue(), 3))  # -> 16.197
```

### ProteinProperties.carbon_oxidation_state

`ProteinProperties.carbon_oxidation_state(gene=None, charge=0.0, **kwargs)`

Madde 11. Ortalama karbon oksidasyon durumu (Dick): (Z − nH + 3nN + 2nO + 2nS) / nC. `charge` net yük Z'dir, nötr protein için 0. `**kwargs` `formula`'ya gider.

```python
print(round(pp.carbon_oxidation_state(), 4))  # -> -0.2011
```

## 4. E-C: yük

### ProteinProperties.termini

`ProteinProperties.termini(gene=None)`

Her zincirin (ilk, son) kalıntısı. Bjellqvist pKa kümesi uç kalıntıya baktığı için `net_charge` bunu kullanır.

```python
print(pp.termini())  # -> [('M', 'G')]
```

### ProteinProperties.charged_counts

`ProteinProperties.charged_counts(gene=None, positive=POSITIVE_RESIDUES, negative=NEGATIVE_RESIDUES)`

Madde 12. Yüklü kalıntı sayıları. `positive` ve `negative` varsayılan olarak `KR` ve `DE`'dir.

```python
print(pp.charged_counts())  # -> {'positive': 11, 'negative': 11}
```

### ProteinProperties.simple_net_charge

`ProteinProperties.simple_net_charge(gene=None, positive=POSITIVE_RESIDUES, negative=NEGATIVE_RESIDUES)`

Madde 13. Kaba net yük (pH 7): (R + K) − (D + E). His ve uçlar ihmal edilir. `positive` ve `negative` kümeleri değiştirir.

```python
print(pp.simple_net_charge())  # -> 0
```

### ProteinProperties.fcr

`ProteinProperties.fcr(gene=None, positive=POSITIVE_RESIDUES, negative=NEGATIVE_RESIDUES)`

Madde 14. Yüklü kalıntı oranı: (R + K + D + E) / N. `positive` ve `negative` kümeleri değiştirir.

```python
print(round(pp.fcr(), 3))  # -> 0.289
```

### ProteinProperties.ncpr

`ProteinProperties.ncpr(gene=None, positive=POSITIVE_RESIDUES, negative=NEGATIVE_RESIDUES)`

Madde 15. Kalıntı başına net yük: [(R + K) − (D + E)] / N. `positive` ve `negative` kümeleri değiştirir.

```python
print(round(pp.ncpr(), 3))  # -> 0.0
```

### ProteinProperties.net_charge

`ProteinProperties.net_charge(gene=None, ph=7.0, pka=DEFAULT_PKA_SET, positive=IONIZABLE_POSITIVE, negative=IONIZABLE_NEGATIVE)`

Madde 16. pH'a bağlı net yük, Henderson-Hasselbalch ile. Uç gruplar serbest, bütün Cys'ler indirgenmiş sayılır. `ph` hesap pH'ıdır. `pka` `PKA_SETS` içindeki kümenin adı ya da aynı anahtarlarla sözlüktür; geçersiz adda `ValueError` verilir. `positive` ve `negative` iyonlaşan yan zincirlerdir.

```python
print(round(pp.net_charge(ph=7.0), 3), round(pp.net_charge(ph=7.0, pka="EMBOSS"), 3))  # -> -0.398 0.226
```

### ProteinProperties.isoelectric_point

`ProteinProperties.isoelectric_point(gene=None, pka=DEFAULT_PKA_SET, tolerance=PI_TOLERANCE, ph_range=PH_RANGE)`

Madde 17. Teorik pI: Q(pH) = 0 olan pH, `ph_range` aralığında `tolerance` duyarlılıkla ikiye bölerek bulunur. Q iki uçta aynı işaretliyse yakın uç döner. `pka` `net_charge` ile aynıdır; varsayılanlar `DEFAULT_PKA_SET` (Bjellqvist), `PH_RANGE` (0-14) ve `PI_TOLERANCE` (0,0001) sabitleridir.

```python
print(round(pp.isoelectric_point(), 2), round(pp.isoelectric_point(pka="Lehninger"), 2))  # -> 6.56 7.67
```

### ProteinProperties.titration_curve

`ProteinProperties.titration_curve(gene=None, pka=DEFAULT_PKA_SET, step=0.1, ph_range=PH_RANGE)`

Madde 18. `ph_range` içinde `step` aralıklarla (pH, Q) çiftleri. `pka` `net_charge` ile aynıdır. `step` 0 ya da negatifse `ValueError` verilir.

```python
egri = pp.titration_curve(step=7.0)
print([(ph, round(q, 2)) for ph, q in egri])  # -> [(0.0, 13.0), (7.0, -0.4), (14.0, -12.96)]
```

## 5. E-D: hidropati

### ProteinProperties.gravy

`ProteinProperties.gravy(gene=None)`

Madde 19. Kyte-Doolittle değerlerinin ortalaması. Pozitif değer hidrofobik proteini gösterir.

```python
print(round(pp.gravy(), 3))  # -> -0.489
```

### ProteinProperties.hydropathy

`ProteinProperties.hydropathy(gene=None, scale='Kyte-Doolittle')`

Madde 20. Seçilen ölçeğin ortalaması. `scale` `HYDROPATHY_SCALES` içindeki ad ya da sözlüktür, geçersiz adda `ValueError` verilir.

```python
print(round(pp.hydropathy(scale="Hopp-Woods"), 3))  # -> 0.274
```

### ProteinProperties.aromaticity

`ProteinProperties.aromaticity(gene=None, aromatic=AMINO_ACID_GROUPS['aromatic'])`

Madde 21. (F + Y + W) / N. `aromatic` kalıntı kümesidir.

```python
print(round(pp.aromaticity(), 3))  # -> 0.039
```

### ProteinProperties.aliphatic_index

`ProteinProperties.aliphatic_index(gene=None, v=ALIPHATIC_V, il=ALIPHATIC_IL)`

Madde 22. Ikai indeksi: X(A) + a·X(V) + b·[X(I) + X(L)], X mol yüzdesidir. `v` (`ALIPHATIC_V`, 2,9) ve `il` (`ALIPHATIC_IL`, 3,9) katsayılardır. Isıl kararlılıkla ilişkilendirilir.

```python
print(round(pp.aliphatic_index(), 2))  # -> 100.0
```

### ProteinProperties.boman_index

`ProteinProperties.boman_index(gene=None, scale=BOMAN_SCALE)`

Madde 23. Boman indeksi (kcal/mol): Σ(n_i × ΔG_i) / N. Yaklaşık 2,48 üstü yüksek bağlanma potansiyeli sayılır. `scale` varsayılan olarak `BOMAN_SCALE`'dir.

```python
print(round(pp.boman_index(), 3))  # -> 1.974
```

## 6. E-E: optik ve hacim

### ProteinProperties.extinction_coefficient

`ProteinProperties.extinction_coefficient(gene=None, reduced=False, coefficients=EXTINCTION_COEFFICIENTS_280)`

Madde 24. 280 nm molar ekstinksiyon: 5500·nTrp + 1490·nTyr + 125·nCystine. `reduced=True` bütün Cys'leri serbest sayar, `False` hepsini eşleşmiş varsayar; ExPASy ProtParam iki değeri de verir. `coefficients` tabloyu değiştirir.

```python
ins = ProteinProperties(PRT("CCWYY", name="p"))
print(ins.extinction_coefficient(), ins.extinction_coefficient(reduced=True))  # -> 8605 8480
```

### ProteinProperties.absorbance_01

`ProteinProperties.absorbance_01(gene=None, reduced=False, coefficients=EXTINCTION_COEFFICIENTS_280)`

Madde 25. %0,1 (1 mg/mL) çözeltinin 280 nm absorbansı: ε / MW. `reduced` ve `coefficients` önceki metotla aynıdır.

```python
print(round(pp.absorbance_01(), 3))  # -> 0.174
```

### ProteinProperties.partial_specific_volume

`ProteinProperties.partial_specific_volume(gene=None, volumes=RESIDUE_SPECIFIC_VOLUMES, masses=RESIDUE_AVERAGE_MASS)`

Madde 26. Kısmi spesifik hacim (cm³/g), Cohn-Edsall: Σ(n_i·M_i·v_i) / Σ(n_i·M_i). `volumes` ve `masses` tabloları değiştirir.

```python
print(round(pp.partial_specific_volume(), 4))  # -> 0.7483
```

### ProteinProperties.mean_residue_volume

`ProteinProperties.mean_residue_volume(gene=None, volumes=RESIDUE_VOLUMES)`

Madde 27. Ortalama kalıntı hacmi (Å³, Zamyatnin). `volumes` tabloyu değiştirir.

```python
print(round(pp.mean_residue_volume(), 2))  # -> 137.14
```

### ProteinProperties.mean_residue_surface

`ProteinProperties.mean_residue_surface(gene=None, surfaces=RESIDUE_MAX_ASA)`

Madde 27. Ortalama kalıntı yüzey alanı (Å², en büyük ASA, Tien ve ark. 2013). `surfaces` tabloyu değiştirir.

```python
print(round(pp.mean_residue_surface(), 2))  # -> 195.97
```

## 7. E-F: kararlılık ve çözünürlük

### ProteinProperties.ivywrel

`ProteinProperties.ivywrel(gene=None, residues=IVYWREL_RESIDUES)`

Madde 28. IVYWREL oranı: (I + V + Y + W + R + E + L) / N. Termofilik proteinlerde yüksektir. `residues` kümeyi değiştirir.

```python
print(round(pp.ivywrel(), 3))  # -> 0.408
```

### ProteinProperties.cvp_bias

`ProteinProperties.cvp_bias(gene=None, charged=CVP_CHARGED, polar=CVP_POLAR)`

Madde 29. CvP yanlılığı (%): [(D + E + K + R) − (N + Q + S + T)] / N × 100. `charged` ve `polar` kümeleri değiştirir.

```python
print(round(pp.cvp_bias(), 2))  # -> 5.26
```

### ProteinProperties.solubility_cv

`ProteinProperties.solubility_cv(gene=None, turn=TURN_FORMING_RESIDUES, lambda_turn=WH_LAMBDA_TURN, lambda_charge=WH_LAMBDA_CHARGE, offset=WH_CHARGE_OFFSET)`

Madde 30. Wilkinson-Harrison kanonik değişkeni: λ1·(N + G + P + S)/N + λ2·|NCPR − 0,03|. `turn`, `lambda_turn`, `lambda_charge` ve `offset` sırasıyla `TURN_FORMING_RESIDUES`, `WH_LAMBDA_TURN`, `WH_LAMBDA_CHARGE` ve `WH_CHARGE_OFFSET` varsayılanlarını değiştirir. Sabitlerin doğrulanma durumu için [tables.md](tables.md#8-kaynağıyla-doğrulanması-gereken-değerler).

```python
print(round(pp.solubility_cv(), 3))  # -> 1.956
```

### ProteinProperties.solubility_probability

`ProteinProperties.solubility_probability(gene=None, discriminant=WH_DISCRIMINANT)`

Madde 30. Wilkinson-Harrison tahmini: `cv`, `soluble` (CV − CV' <= 0) ve `probability` (tahmin edilen sınıfın olasılığı, `WH_PROBABILITY` katsayılarıyla). `discriminant` CV' eşiğidir, varsayılanı `WH_DISCRIMINANT`. Ubikitinin çözünür olduğu bilindiği halde model "çözünmez" der: model kabadır, sonucu tek başına kanıt saymayın.

```python
print(pp.solubility_probability()["soluble"])  # -> False
```

### ProteinProperties.composition_entropy

`ProteinProperties.composition_entropy(gene=None)`

Madde 31. Bileşimin Shannon entropisi (bit): −Σ p_i log2 p_i. En çok log2 20 ≈ 4,32'dir; düşük değer düşük karmaşıklık gösterir.

```python
print(round(pp.composition_entropy(), 3))  # -> 3.917
```

## 8. E-G: düzensizlik

### ProteinProperties.uversky

`ProteinProperties.uversky(gene=None, hydropathy=NORMALIZED_KYTE_DOOLITTLE)`

Madde 32. Uversky grafiğindeki konum. `mean_net_charge` (|NCPR|), `mean_hydropathy`, `boundary` ve `disordered` döner. Sınır ⟨H⟩ = (⟨R⟩ + `UVERSKY_OFFSET`) / `UVERSKY_SLOPE` doğrusudur. `hydropathy` 0-1 aralığına normalize Kyte-Doolittle ölçeğidir (`NORMALIZED_KYTE_DOOLITTLE`).

```python
print(pp.uversky()["disordered"])  # -> False
```

### ProteinProperties.top_idp

`ProteinProperties.top_idp(gene=None, scale=TOP_IDP_SCALE)`

Madde 33. TOP-IDP ortalaması; pozitif değer düzensizlik eğilimidir. `scale` ölçeği değiştirir.

```python
print(round(pp.top_idp(), 3))  # -> 0.116
```

### ProteinProperties.composition_features

`ProteinProperties.composition_features(gene=None)`

Madde 34. Bileşim tabanlı çözünürlük ve toplanma modellerinin kullandığı özellikler: amino asit oranları, `K-R`, `D-E`, `K+R`, `D+E`, `K+R-D-E`, `K+R+D+E`, `F+W+Y`, `pI`, `abs_charge`, `gravy`, `entropy`, `top_idp`, `fcr`, `ncpr`.

```python
print(round(pp.composition_features()["K+R"], 3))  # -> 0.145
```

### ProteinProperties.linear_score

`ProteinProperties.linear_score(gene=None, weights=None, intercept=0.0)`

Madde 34. Doğrusal skor: `intercept` + Σ ağırlık × özellik. `weights` özellik adından katsayıya giden sözlüktür ve zorunludur: yayımlanmış modellerin (ör. Protein-Sol) katsayıları pakette yoktur, kaynağından alınıp verilmelidir. `weights` verilmezse `ValueError`, bilinmeyen özellik adında `KeyError` verir.

```python
print(round(pp.linear_score(weights={"gravy": 2.0, "fcr": -1.0}, intercept=0.5), 3))  # -> -0.768
```

## 9. İç yardımcılar ve sabitler

| Ad | Görevi |
|---|---|
| `_KD` | Normalizasyon için Kyte-Doolittle tablosunun kısa adı |
| `_Proteins` | `PRT` listesini `CUB`'ın sayım arayüzüne uyarlar, böylece ölçüm kodu iki kaynakta aynı kalır |
| `_Proteins.__init__` | Tipi ve yinelenen adı denetler, 20 standart amino asidi sayar, uçları saklar |
| `_Proteins.names` | Protein adları |
| `_Proteins._selected` | `gene` geçerliyse onu, `None` ise hepsini verir; bilinmeyen adda `KeyError` |
| `_Proteins.amino_acid_counts` | `CUB.amino_acid_counts` karşılığı |
| `_Proteins.termini` | `CUB.termini` karşılığı |
| `_Proteins.chain_count` | `CUB.chain_count` karşılığı |
| `ProteinProperties._counts` | Kaynağın `amino_acid_counts` çağrısının kısayolu |
| `ProteinProperties._pka_set` | `pka` adını `PKA_SETS` kümesine çevirir, sözlüğü olduğu gibi bırakır |
| `ProteinProperties.__repr__` | `ProteinProperties(1 protein)` |
