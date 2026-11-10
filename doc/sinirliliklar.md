# Sınırlılıklar ve sık hatalar

Paketin bilinçli olarak yapmadıkları ve sık karşılaşılan hatalar bu sayfada toplanmıştır. Bir sonucu raporlamadan önce ilgili bölüme bakın.

<!-- toc -->
- [1. Diziler ve dönüşümler](#1-diziler-ve-dönüşümler)
- [2. CUB ölçümleri](#2-cub-ölçümleri)
- [3. Diferansiyel kodon kullanımı](#3-diferansiyel-kodon-kullanımı)
- [4. HPO](#4-hpo)
- [5. ClinVar](#5-clinvar)
- [6. Sabitler](#6-sabitler)
- [7. Sık hatalar](#7-sık-hatalar)
<!-- /toc -->

## 1. Diziler ve dönüşümler

FASTA kayıplıdır: `to_fasta` yalnızca adı ve diziyi yazar, tip, `strand`, `genomic_strand` ve features dosyaya geçmez. FASTA'dan okunan her dizi 5'→3' kabul edilir, bu yüzden 3'→5' yazılmış bir dizi FASTA'ya yazılıp okunursa yönü yanlış yorumlanır. Kayıpsız saklama için CRN kullanın. FASTA adı tek kelime olmalıdır, boşluklu ad `to_fasta`'da hata verir.

Eşitlik yalnızca tip ve diziye bakar: ad, features ve strand karşılaştırılmaz, aynı harf dizisine sahip 5'→3' ve 3'→5' iki DNA eşit sayılır. Diziler hash'lenemez, küme ya da sözlük anahtarı olamaz. Ters dilim yönü değiştirmez: `s[::-1]` harfleri ters çevirir ama `strand` aynı kalır.

RNA'da T kabul edilmez; T içeren dosyaları genel `Sequence` ile okuyup çevirin. Değiştirilmiş nükleotitler (inozin, psödouridin) desteklenmez. Diziler büyük harfe çevrildiği için tekrar bölgelerini gösteren küçük harf maskesi korunmaz.

Belirsiz kodonlar her zaman `X` olur. `UAR` (her durumda stop) ya da `GCN` (her durumda alanin) gibi aslında belirli kodonlar da `X` olur ve `find_orfs` belirsiz bazlı stopları stop olarak tanımaz. Yalnızca iki genetik kod vardır: standart (NCBI 1) ve omurgalı mitokondriyal (NCBI 2).

`translation` ilk nükleotitten okur, başlangıç kodonu ya da stop aramaz, alternatif başlangıç kodonunu (CUG) M'ye çevirmez. Selenosistein ve pirolizin hiçbir zaman üretilmez; bu yüzden selenoproteinlerin ana transkripti proteinle eşleşmez: `gc.canonical("GPX1").isoforms()` boş döner. `Hpo.to_genes` böyle genlerde başka bir transkripte geçer ve `canonical_tag` değerini `other_transcript` yapar.

`find_orfs` dairesel dizileri desteklemez: plazmid ya da mitokondriyal DNA'da başlangıç noktasının üzerinden geçen ORF bulunmaz. İç içe ORF'lerin hepsi döner, kısmi ORF alınmaz. Uzun genomik dizilerde çok sayıda kısa ORF çıkar, `minimum_length_of_orf` kullanın.

## 2. CUB ölçümleri

Hesaplar verilen genetik koda bağlıdır. Mitokondriyal genler standart kodla çevrilemez; `Hpo.to_genes` bunları atlar ve bildirir, kendi kurduğunuz sette `genetic_code` vermeniz gerekir.

CDC'de bootstrap yoktur ve arka plan konum baz oranlarından kurulur (bkz. [divergence.md](cub/divergence.md)). tAI'nin wobble cezaları *S. cerevisiae* için optimize edilmiştir, insan için ayrıca optimize edilmemiştir. Bu nedenle tAI'yi mutlak değer olarak değil, aynı koşulda hesaplanmış genler arasında karşılaştırın.

Küçük gen setlerinde ENC `nan` dönebilir, CPB ise referansta bulunmayan çiftleri atar. CPB için referansın bütün genom olması gerekir.

## 3. Diferansiyel kodon kullanımı

Fisher sonucu tek başına kanıt değildir, [nedeni burada](cub/differential.md#3-neden-permütasyon-gerekir). Rapor edilecek liste `robust`'tur. Grup küçükse permütasyon anlamlılığa matematiksel olarak ulaşamaz: 4'e karşı 4 gende ayrık bölünme sayısı 35'tir, en küçük p 1/35 olur. `interpret` bu durumu söyler, `DifferentialCodonUsage._smallest_p` sınırı hesaplar.

Permütasyon uzunluk ve bölgesel GC farkını gidermez. GC3s farkı olan iki grupta çıkan kodonlar çeviri seçilimi değil mutasyon yanlılığı olabilir; kontrol grubunu uzunluk ve GC'ye göre eşleştirin. `plot_volcano` Fisher sonucunu çizer, `robust` işaretini göstermez.

## 4. HPO

Gen-fenotip ilişkisi hastalık üzerindendir ve ilişki etiketleri güvenilir değildir; ayrıntı [hpo.md](fetch/hpo.md#1-yorumlarken-bilinmesi-gerekenler). `min_frequency` sıklığı bilinmeyen satırları da eler (v2026-09-01 fenotip satırlarının %14,7'si). Eşik gen sayısını hızla düşürür: `Ataxia` 1119 gen verirken `min_frequency=0.3` ile 669 gen kalır.

`Gencode` kapsamı dışındaki semboller, protein kodlamayan genler ve ORF'u proteine eşleşmeyen genler `to_genes` tarafından atlanır. `pack` raporundaki `skipped` sözlüğü nedenleri listeler; bir grubun ne kadarının düştüğünü her zaman kontrol edin.

## 5. ClinVar

`significance="pathogenic"` yalnızca gönderenlerin uzlaştığı varyantları verir. Çelişkili sınıflandırmalar (`conflicting`, 2026-11'de 183.872 varyant) patojen sayılmaz; bir genin patojen varyant sayısı bu yüzden gerçek yükün alt sınırıdır. Yıldızsız (0) sınıflandırmalar da sonuca girer, çalışma kümesi için `min_stars` verin. Ayrıntı [clinvar.md](fetch/clinvar.md#1-yorumlarken-bilinmesi-gerekenler).

`effect` varyant adındaki HGVS protein yazımından kurallarla çıkarılır, etki tahmini yapılmaz. Ekleme bölgesini bozan eşanlamlı bir varyant `synonymous` görünür; 2026-11'de patojen sınıfta 474 eşanlamlı varyant vardır ve bunların patojenliği amino asit değişiminden gelmez. Protein yazımı olmayan varyantlarda `effect` boştur ve `effect` süzgeci bunları dışarıda bırakır.

`max_genes` varsayılan olarak 1'dir. Bu, büyük delesyonların gen sayımlarını şişirmesini önler, ama iki komşu geni birlikte listeleyen küçük varyantları da çıkarır. Sayımlar beklenenden düşükse `max_genes=None` ile karşılaştırın.

Hastalık adları gönderildiği biçimdedir, aynı hastalık birden çok adla geçer. Ad parçasıyla arama hem eksik hem fazla sonuç verebilir; yayımlanacak analizde OMIM ya da MONDO kimliği kullanın.

NCBI aylık arşivler için MD5 yayımlamaz. İndirme boyut ve gzip CRC32 ile denetlenir, SHA-256 manifest'e yazılır; aynı sürümün başka bir makinede aynı veriyi verdiği bu özetle karşılaştırılarak doğrulanır.

## 6. Sabitler

Bazı ölçeklerin ve Wilkinson-Harrison sabitlerinin özgün kaynakla doğrulanması gerekir: [tables.md](tables.md#8-kaynağıyla-doğrulanması-gereken-değerler).

## 7. Sık hatalar

| Belirti | Neden | Çözüm |
|---|---|---|
| `Her gende bir genomik DNA, bir transkript ve bir protein olmalı` | Çok izoformlu gen `CUB`'a verilmiş | `gene.isoforms()` ya da `gc.canonical(...).isoforms()` |
| `genes bir liste olmalı` | Tek gen liste olmadan verilmiş | `CUB([gene])` |
| `Aynı adlı genler var` | İki kaynak ya da iki izoform aynı adla gelmiş | İzoform adlarını kullanın ya da birini çıkarın |
| `'PDE8B' adında birden çok gen var` | `AmbiguousGeneError` | Mesajdaki `gene_id`'yi kullanın |
| `CRN dosyasına yazılamayan değer tipi` | Features'ta numpy tam sayısı, `dict` ya da `set` var | `int()` ile çevirin ya da listeye dönüştürün |
| `bu bir liste dosyası; read_crn_list ile açın` | Yanlış okuyucu | Dosyayı yazan fonksiyonun çifti ile okuyun |
| `method='chi2' için reference ve background farklı gen setleri olmalı` | `background` verilmemiş | İkinci bir `CUB` verin |
| `HPO sürüm bilgisi alınamadı (HTTP 403)` | GitHub API saatlik sınırı | Bekleyin ya da yerel sürümle çalışın |
| `ENC` `nan` | Gen çok kısa | Beklenen davranış; bkz. [enc.md](cub/enc.md) |
