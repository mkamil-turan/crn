# `crn.fetch._geneset`: gen listesi veren kaynakların ortak tabanı

Kaynak: [`crn/fetch/_geneset.py`](../../crn/fetch/_geneset.py)

`Hpo` ve `Clinvar` farklı sorulara yanıt verir. Biri fenotipten, öteki klinik varyanttan yola çıkar, ama ikisinin de sonucu bir gen listesidir. Bu listeyi `CUB`'a verilebilecek `Gene` nesnelerine çevirmek, dizilere ayırmak ve gruplar hâlinde tek CRN dosyasına yazmak iki kaynakta aynı iştir. `GeneSetSource` bu ortak işi taşır. Kullanıcı bu sınıfı doğrudan kullanmaz; metotlarına `Hpo` ve `Clinvar` üzerinden ulaşır.

<!-- toc -->
- [1. Alt sınıfın sağladıkları](#1-alt-sınıfın-sağladıkları)
- [2. Sınıf](#2-sınıf): [`GeneSetSource`](#genesetsource) · [`pack`](#genesetsourcepack) · [`write_tsv`](#genesetsourcewrite_tsv)
- [3. İç yardımcılar](#3-iç-yardımcılar)
<!-- /toc -->

## 1. Alt sınıfın sağladıkları

| Ad | Görevi | `Hpo` | `Clinvar` |
|---|---|---|---|
| `genes` | Filtrelere uyan gen satırları; her satırda `ncbi_gene_id` ve `gene_symbol` bulunur | terimden | varyanttan |
| `_describe` | Kurulan `Gene`'in features'ına kaynağa özgü bilgiyi yazar | HPO terimleri ve hastalıklar | ClinVar özeti |
| `release` | Kullanılan veri sürümü | `v2026-09-01` | `2026-11` |
| `FEATURE_KEYS` | `to_sequences`'ın genden diziye taşıdığı anahtarlar | `hpo_id`, `hpo_name`, `disease_id`, `disease_name` | `clinvar_condition`, `clinvar_variants`, `clinvar_pathogenic` |
| `RELEASE_KEY` | Sürümün yazıldığı features anahtarı | `hpo_release` | `clinvar_release` |
| `GROUP_FILTER` | `pack`'te grup değeri sözlük değilse aktarıldığı filtre | `term` | `condition` |

`SEQUENCE_KINDS` (`dna`, `rna`, `protein`) `to_sequences`'ın, `OVERLAP_POLICIES` (`exclude`, `first`, `error`) `pack`'in kabul ettiği değerlerdir.

## 2. Sınıf

### GeneSetSource

`GeneSetSource(sys_dir=None, verbose=True)`

`LocalSource`'tan türer. `sys_dir` veri klasörü, `verbose` ilerleme çubuklarıdır. Kendisi bir veri kaynağı değildir: `genes`, `_describe` ve `release` alt sınıfta tanımlanmadan çağrılırsa `NotImplementedError` verir.

### GeneSetSource.pack

`GeneSetSource.pack(groups, file_name, source=None, key='group', overlap='exclude', **filters)`

Gen gruplarını tek CRN listesine yazar; dosya `DifferentialCodonUsage.from_crn` ile doğrudan açılır. `groups` grup adından sorgu değerine ya da gruba özel filtre sözlüğüne giden eşlemedir. Değer sözlük değilse `GROUP_FILTER` filtresine aktarılır, yani `Hpo`'da terim, `Clinvar`'da koşul sayılır. `file_name` yazılacak dosya, `source` genlerin kurulacağı kaynaktır (`Gencode`, `Ncbi` ya da `None` ile `sys/` içindeki GENCODE). `key` grup anahtarıdır. `overlap` birden çok gruba giren genler içindir: `"exclude"` hepsinden çıkarır (varsayılan), `"first"` ilk grupta bırakır, `"error"` hata verir. `**filters` bütün gruplara uygulanır.

Dönen rapor `file`, `release`, `groups` (grup → gen sayısı), `overlap` (çakışan semboller) ve `skipped` (atlanan gen → neden) alanlarını içerir. `overlap` geçersizse, gruplar çakışıyorsa (`"error"`) ya da bir grup boş kalırsa `ValueError` verir. Örnekler [hpo.md](hpo.md#hpopack) ve [clinvar.md](clinvar.md#clinvarpack) belgelerindedir.

### GeneSetSource.write_tsv

`statik: GeneSetSource.write_tsv(rows, file_name)`

Sorgu sonucunu (`rows`) `file_name` dosyasına sekmeyle ayrılmış yazar. Satır yoksa `ValueError` verir.

## 3. İç yardımcılar

| Ad | Görevi |
|---|---|
| `GeneSetSource.release` | Alt sınıfın sürüm özelliği |
| `GeneSetSource.genes` | Alt sınıfın gen sorgusu |
| `GeneSetSource._describe` | Alt sınıfın features kancası |
| `GeneSetSource._check_choices` | Değerin izin verilen listede olduğunu denetler, değilse seçenekleri sayan `ValueError` verir |
| `GeneSetSource._row_key` | Gen satırının kimliği: NCBI Gene ID, yoksa sembol. ClinVar'da çok genli varyantlarda görülen sembollerin kimliği olmayabilir |
| `GeneSetSource._gencode` | Varsayılan `Gencode` nesnesini bir kez açar |
| `GeneSetSource._canonical` | Ana transkripti alır, aynı adlı kayıtları ayırır: tek kodlayan kayıt (CAST) ya da HGNC üzerinden Ensembl kimliği (SHOX, PDE8B) |
| `GeneSetSource._isoform` | Genin tek DNA, tek transkript, tek proteinli halini kurar; ana transkriptin ORF'u uymazsa (PANK2) başka transkripte geçer. `Ncbi` kaynağında Gene ID yoksa sembol kullanılır |
| `GeneSetSource._build` | Gen satırlarını `Gene` nesnelerine çevirir, `_describe`'ı çağırır, atlananları nedenleriyle toplar |
| `GeneSetSource._to_genes` | `to_genes`'in ortak gövdesi; atlanan genleri `UserWarning` ile bildirir |
| `GeneSetSource._to_sequences` | `to_sequences`'ın ortak gövdesi; `FEATURE_KEYS` ve `RELEASE_KEY` değerlerini diziye taşır |
