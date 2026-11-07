"""Örnek 5: Küçük RNA'lar: miRNA, snRNA, snoRNA ve fazlası.

Senaryo: Küçük RNA dizileme verinizi hizalamak için tiplere göre referans dizileri
ve miRNA ailelerini çıkarmak istiyorsunuz.

Bu örnekte:
  1. Küçük RNA tiplerinin sayısı, uzunluğu ve GC içeriği hesaplanır.
  2. Gen adlarından miRNA aileleri (MIR ve let-7) gruplanır.
  3. let-7 öncüllerindeki ortak 21-merler sliding ile bulunur.
  4. let-7 öncülleri ve Mt-tRNA'lar tek CRN listesine, miRNA'lar FASTA'ya yazılır.

Çıktılar: test/gencode/cikti/gencode_orn_5/
Not: GENCODE miRNA kayıtları olgun miRNA değil, saç tokası öncülleridir.
"""
import re
from collections import Counter
from pathlib import Path
from statistics import median

from crn import RNA, write_crn_list, write_fasta
from crn.fetch import Gencode

CIKTI = Path(__file__).parent / "cikti" / "gencode_orn_5"
CIKTI.mkdir(parents=True, exist_ok=True)

TIPLER = ["miRNA", "snRNA", "snoRNA", "scaRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "vault_RNA"]

gc = Gencode(verbose=False)


def rnalar(tip: str) -> list:
    sonuc = []
    for t in gc.transcripts(transcript_type=tip):
        rna = gc.rna(t["transcript_id"])
        rna.name = t["transcript_name"].rsplit("-", 1)[0]
        sonuc.append(rna)
    return sonuc


def gc_yuzde(rna: RNA) -> float:
    return 100 * (rna.seq.count("G") + rna.seq.count("C")) / len(rna)


kumeler = {tip: rnalar(tip) for tip in TIPLER}
print(f"{'Tip':<12}{'Kayıt':>7}{'Ortanca':>9}{'En kısa':>9}{'En uzun':>9}{'Ort.GC%':>10}")
for tip, liste in kumeler.items():
    boylar = [len(r) for r in liste]
    print(f"{tip:<12}{len(liste):>7,}{median(boylar):>9.0f}{min(boylar):>9}{max(boylar):>9,}"
          f"{median(gc_yuzde(r) for r in liste):>10.1f}")

aileler = Counter()
for rna in kumeler["miRNA"]:
    eslesme = re.match(r"(MIRLET7|MIR\d+)", rna.name)
    if eslesme:
        aileler[eslesme.group(1)] += 1
print("\nEn kalabalık miRNA aileleri:", ", ".join(f"{a} ({n})" for a, n in aileler.most_common(6)))

let7 = [r for r in kumeler["miRNA"] if r.name.startswith("MIRLET7")]
print(f"\nlet-7 öncülü: {len(let7)}")
kmerler = Counter()
for rna in let7:
    kmerler.update({p.seq for p in rna.sliding(21)})
print("En yaygın 21-mer:")
for kmer, sayi in kmerler.most_common(3):
    print(f"  {kmer}  {sayi}/{len(let7)} öncül")

write_crn_list(let7 + kumeler["Mt_tRNA"], CIKTI / "let7_ve_mt_trna.crn")
adlar = Counter()
miRNA = []
for rna in kumeler["miRNA"]:      # FASTA'da ad boşluksuz ve mümkünse tekil olsun
    adlar[rna.name] += 1
    rna.name = f"{rna.name}_{adlar[rna.name]}" if adlar[rna.name] > 1 else rna.name
    miRNA.append(rna)
write_fasta(miRNA, CIKTI / "mirna_oncu.fasta")
print("\nÇıktılar:", CIKTI)
