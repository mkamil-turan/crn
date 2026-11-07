"""Diferansiyel kodon örneklerinin ortak araçları (dco_1 ... dco_6 bunu kullanır).

Bu dosya tek başına çalıştırılmaz. Her örnek aynı yedi adımı izlediği için adımlar
burada toplandı. Örnekleri okumadan önce bu açıklamayı okumak işi kolaylaştırır.

SORU
    Bir hastalığın genleri, eş anlamlı kodonları diğer genlerden farklı mı seçiyor?
    Protein aynı kalsa bile kodon seçimi çeviri hızını, katlanmayı ve mRNA ömrünü
    etkileyebilir. Bu nedenle kodon kullanımı, genetik değişimin ötesinde bir
    düzenleme katmanı olarak incelenir.

NEDEN DOĞRUDAN KARŞILAŞTIRILMAZ
    Fark başka sebeplerden kolayca doğar. Üç ana tuzak ve burada alınan önlem:

    1. Gen uzunluğu ve bölgesel GC. Otizm genleri gibi birçok sinir sistemi geni çok
       uzundur ve GC'si düşük bölgelerdedir. Üçüncü kodon konumu bu bölgesel GC'yi
       izler. Sonuç seçilim gibi görünür ama yalnızca genomdaki konumu yansıtır.
       Önlem: her hastalık genine CDS uzunluğu ve genomik GC'si en yakın kontrol
       geni seçilir (``eslestir``).

    2. Kodon havuzlama. Fisher testi gruptaki tüm kodonları toplar, her kodonu
       bağımsız sayar. Tek bir uzun gen sonucu sürükleyebilir.
       Önlem: gen etiketleri binlerce kez karıştırılır (``dcu.permutation``). Fark ancak
       gerçek grup, rastgele gruplardan belirgin biçimde ayrılıyorsa "sağlam" sayılır.
       Bu önlem isteğe bağlı değildir: dco_6'daki kalibrasyon, kontrol havuzunun rastgele
       ikiye bölünmesinde bile Fisher'ın onlarca kodonu anlamlı bulduğunu gösterir.

    3. Kontrol seçimi. Kontrol genleri de Mendel hastalık genleridir, böylece
       "iyi çalışılmış gen" etkisi iki grupta da bulunur. Sinir sistemi ve göz notu
       taşımazlar. Ama ifade düzeyi eşleştirilmemiştir, çünkü veri setimizde yoktur.

YEDİ ADIM
    1. HPO'dan hastalık gen grubu çıkarılır.
    2. Kontrol havuzundan eşleştirilmiş kontrol seçilir, denge raporlanır.
    3. İki grup tek CRN dosyasına yazılır (grup bilgisi features'ta).
    4. ``DifferentialCodonUsage`` ile üç düzeyde kodon testi yapılır.
    5. Aile düzeyi sonuçlar gen düzeyinde permütasyonla sınanır (``dcu.robust``).
    6. CUB parametreleri gen başına karşılaştırılır (Mann-Whitney, BH).
    7. Klinik özet ve yorum yazılır.

CUB PARAMETRELERİ NE ANLATIR
    L          Kodon sayısı (protein uzunluğu).
    GC%, GC3s  Genel GC ve üçüncü konumdaki sinonim GC. Kodon yanlılığının en güçlü
               belirleyicisi bölgesel GC'dir.
    ENC        Etkin kodon sayısı. 20 tek kodon seçimi, 61 tamamen eşit kullanım.
               Düşük ENC güçlü yanlılık demektir.
    ENC_ratio  GC3s'nin beklettiği ENC'ye göre oran. Bileşimin açıklamadığı yanlılığı
               gösterir.
    SCUO, CDC  Yanlılığın bilgi teorisine dayalı genel ölçüleri.
    RCBS       Baz bileşimine göre düzeltilmiş kodon yanlılığı.
    Fop        Optimal kodon oranı. Referans: ribozomal protein genleri.
    CAI        Yüksek ifadeli genlere uyum.
    CBI        Yanlılığın rastgele kullanıma göre konumu.
    tAI        Kodonların tRNA gen kopya sayısına uyumu. Düşükse çeviri, en azından
               tRNA havuzu açısından daha yavaş olabilir.

SONUÇ NASIL OKUNUR
    Cliff delta (-1 ile 1): hastalık geninin kontrolden yüksek çıkma eğilimi. 0.15'in
    altı önemsiz, 0.33 orta, 0.47 üstü büyük kabul edilir (Romano ölçütleri). q değeri
    anlamlılığı, delta büyüklüğü verir. Bulgular ilişkidir, hastalık nedeni değildir.
    Bütün çıktılar yeni hipotezdir, klinik karar için kullanılmaz.
"""
from __future__ import annotations

import csv
import math
import random
import re
import urllib.request
import warnings
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import false_discovery_control, mannwhitneyu

from crn import Gene, write_crn_list
from crn.cub import CUB, DifferentialCodonUsage, read_trna_counts
from crn.fetch import Gencode, Hpo

CIKTI = Path(__file__).parent / "cikti"
TRNA_URL = "https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa"
PARAMETRELER = ("L", "GC%", "GC3s", "ENC", "ENC_ratio", "SCUO", "CDC", "RCBS", "Fop", "CAI", "CBI", "tAI")
# Kontrol havuzu: fenotip notu olan, ama sinir sistemi ve göz notu taşımayan genler.
KONTROL_DISI = ["Abnormality of the nervous system", "Abnormality of the eye"]
# HPO'nun ilişki etiketi tekrar genişlemesi ve risk genlerinde güvenilir değildir: SCA2
# "POLYGENIC" görünür, LRRK2 ve GBA1 de öyledir. Yalnızca MENDELIAN seçilirse temel
# genler düşer. Bu yüzden türü bilinmeyen (UNKNOWN) Orphanet kayıtları dışlanır.
ILISKI = ["MENDELIAN", "POLYGENIC"]
TOHUM = 2026
PERMUTASYON = 5000


@dataclass
class Kaynaklar:
    """Her örnekte gereken veri kaynakları."""
    hpo: Hpo
    gencode: Gencode
    referans: CUB                   # CAI/Fop için yüksek ifadeli referans
    trna: Dict[str, int]            # tAI için antikodon başına tRNA gen sayısı
    havuz: List[Gene]               # eşleştirilecek kontrol adayları


def hazirla() -> Kaynaklar:
    """Kaynakları yükler. Dış veri yalnızca ilk çalıştırmada indirilir."""
    print("Veri kaynakları hazırlanıyor...")
    hpo = Hpo(verbose=False)
    gencode = Gencode(verbose=False)

    # Ribozomal protein genleri insanda en yüksek ifadeli gen ailelerinden biridir;
    # CAI ve Fop için alışılmış referanstır.
    referans = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for satir in gencode.genes(gene_type="protein_coding", name=r"^RP[LS]\d+[AXYL]?\d*$"):
            izoformlar = gencode.canonical(satir["gene_id"]).isoforms()
            if izoformlar:
                referans.append(izoformlar[0])
        havuz = hpo.to_genes("Phenotypic abnormality", exclude=KONTROL_DISI)

    yol = hpo.sys_dir / "gtrnadb" / "hg38-tRNAs.fa"
    if not yol.is_file():
        yol.parent.mkdir(parents=True, exist_ok=True)
        istek = urllib.request.Request(TRNA_URL, headers={"User-Agent": "crn"})
        with urllib.request.urlopen(istek, timeout=60) as yanit:
            yol.write_bytes(yanit.read())
    trna = read_trna_counts(yol)
    print(f"  Referans {len(referans)} ribozomal gen, kontrol havuzu {len(havuz)} gen, "
          f"{sum(trna.values())} tRNA geni\n")
    return Kaynaklar(hpo, gencode, CUB(referans), trna, havuz)


def gen_al(k: Kaynaklar, *args, **filtreler) -> List[Gene]:
    """HPO sorgusundan Gene listesi; atlanan genler (kodlamayan vb.) sessizce düşer."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return k.hpo.to_genes(*args, **filtreler)


def gen_adi(gene: Gene) -> str:
    return gene.features.first("gene_name") or gene.name


def _ozellik(gene: Gene) -> Tuple[float, float]:
    """Eşleştirmede kullanılan özellikler: log10(kodon sayısı), genomik GC."""
    dizi = gene.dna.seq_5to3().upper()
    gc = (dizi.count("G") + dizi.count("C")) / len(dizi)
    orf = next(iter(gene.transcripts.values()))
    uzunluk = (orf.features.first("orf_end") - orf.features.first("orf_start") + 1) / 3
    return math.log10(uzunluk), gc


def eslestir(hasta: Sequence[Gene], havuz: Sequence[Gene], tohum: int = TOHUM
             ) -> Tuple[List[Gene], Dict[str, float]]:
    """Her hastalık genine uzunluk ve genomik GC'si en yakın kontrolü seçer.

    Yöntem: açgözlü en yakın komşu, tekrarsız. Özellikler z-skoruna çevrilir, böylece
    iki özellik eşit ağırlık alır. Sıra sabit tohumla karıştırıldığı için sonuç tekrar
    üretilebilir.

    Returns:
        (kontrol genleri, denge). Denge, standart ortalama farkıdır (SMD): 0.1 altı
        iyi kabul edilir. Küçükse iki grup uzunluk ve GC açısından karşılaştırılabilir.
    """
    hasta_adlari = {gen_adi(g) for g in hasta}
    adaylar = [g for g in havuz if gen_adi(g) not in hasta_adlari]
    h, a = np.array([_ozellik(g) for g in hasta]), np.array([_ozellik(g) for g in adaylar])
    ort, ss = a.mean(axis=0), a.std(axis=0)
    hz, az = (h - ort) / ss, (a - ort) / ss

    bos = np.ones(len(adaylar), dtype=bool)
    secilen = []
    sira = list(range(len(hasta)))
    random.Random(tohum).shuffle(sira)
    for i in sira:
        uzaklik = ((az - hz[i]) ** 2).sum(axis=1)
        uzaklik[~bos] = np.inf
        j = int(uzaklik.argmin())
        bos[j] = False
        secilen.append(j)
    k = a[secilen]
    smd = np.abs(h.mean(axis=0) - k.mean(axis=0)) / np.sqrt((h.var(axis=0) + k.var(axis=0)) / 2)
    return [adaylar[j] for j in secilen], {"uzunluk": float(smd[0]), "GC": float(smd[1])}


def paketle(gruplar: Dict[str, Sequence[Gene]], dosya: Path) -> None:
    """Grupları tek CRN gen listesine yazar; grup adı her genin features'ına girer."""
    dosya.parent.mkdir(parents=True, exist_ok=True)
    liste = []
    for ad, genler in gruplar.items():
        for gen in genler:
            gen.features.set("group", ad)
            liste.append(gen)
    write_crn_list(liste, dosya)


def parametreler(cub: CUB, k: Kaynaklar) -> List[Dict[str, object]]:
    """Gen başına CUB tablosu: GC, ENC, Fop, CAI, CBI, tAI ve diğerleri."""
    return cub.table(reference=k.referans, trna_counts=k.trna, protein=False)


def karsilastir(a: List[Dict[str, object]], b: List[Dict[str, object]]) -> List[Dict[str, float]]:
    """CUB parametrelerini gen düzeyinde karşılaştırır.

    Mann-Whitney U, dağılımın normal olmasını istemez ve her gene eşit oy verir.
    Cliff delta = 2U / (n1 n2) - 1: a'daki bir genin b'deki bir genden büyük olma
    olasılığı ile küçük olma olasılığı arasındaki farktır.
    """
    sonuc = []
    for ad in PARAMETRELER:
        x = [float(s[ad]) for s in a if s.get(ad) is not None and not math.isnan(float(s[ad]))]
        y = [float(s[ad]) for s in b if s.get(ad) is not None and not math.isnan(float(s[ad]))]
        if not x or not y:
            continue
        test = mannwhitneyu(x, y, alternative="two-sided")
        sonuc.append({"parametre": ad, "medyan_a": float(np.median(x)), "medyan_b": float(np.median(y)),
                      "delta": 2 * float(test.statistic) / (len(x) * len(y)) - 1, "p": float(test.pvalue)})
    for satir, q in zip(sonuc, false_discovery_control([s["p"] for s in sonuc], method="bh")):
        satir["q"] = float(q)
    return sonuc


def _buyukluk(delta: float) -> str:
    d = abs(delta)
    return "önemsiz" if d < .147 else "küçük" if d < .33 else "orta" if d < .474 else "büyük"


def klinik_ozet(k: Kaynaklar, genler: Sequence[Gene]) -> Dict[str, object]:
    """Grubun tıbbi profili: hastalık sayısı, kalıtım, başlangıç, sık fenotipler."""
    semboller = [gen_adi(g) for g in genler]
    satirlar = k.hpo.annotations("Phenotypic abnormality", gene=semboller, association=ILISKI)
    kalitim, baslangic = Counter(), Counter()
    hastaliklar = {s["disease_id"] for s in satirlar}
    for hastalik in k.hpo.diseases(gene=semboller, association=ILISKI):
        for ad in filter(None, hastalik["inheritance"].split("|")):
            kalitim[ad] += 1
        for ad in filter(None, hastalik["onset"].split("|")):
            baslangic[ad] += 1
    fenotip: Dict[str, set] = {}
    for s in satirlar:
        fenotip.setdefault(s["hpo_name"], set()).add(s["ncbi_gene_id"])
    sik = sorted(fenotip.items(), key=lambda x: -len(x[1]))[:6]
    return {"hastalik": len(hastaliklar), "kalitim": kalitim.most_common(4),
            "baslangic": baslangic.most_common(4),
            "fenotip": [(ad, len(g), round(100 * len(g) / len(semboller))) for ad, g in sik]}


@dataclass
class Sonuc:
    """Bir hastalık-kontrol karşılaştırmasının bütün çıktıları."""
    ad: str
    kontrol: str
    dcu: DifferentialCodonUsage
    hasta: List[Gene]
    karsi: List[Gene]
    denge: Dict[str, float]
    saglam: List[Dict[str, object]]     # Fisher ve permütasyonda birlikte anlamlı kodonlar
    fisher: Dict[str, int]              # düzey -> anlamlı kodon
    fisher_aile: int
    param: List[Dict[str, float]]
    tablo_a: List[Dict[str, object]]
    tablo_b: List[Dict[str, object]]
    klinik: Dict[str, object]


def analiz(k: Kaynaklar, ad: str, hasta: List[Gene], klasor: str,
           karsi: Optional[List[Gene]] = None, kontrol_adi: str = "kontrol") -> Sonuc:
    """Yedi adımın tümünü çalıştırır, tabloları yazar ve sonuçları ekrana döker.

    Args:
        k: ``hazirla`` sonucu.
        ad: Hastalık grubunun adı (dosya ve tablolarda).
        hasta: Hastalık genleri.
        klasor: ``cikti`` altındaki klasör.
        karsi: Verilirse eşleştirme yapılmaz, bu grupla doğrudan karşılaştırılır
            (ör. iki hastalık arasında). Bu durumda uzunluk ve GC eşlenmemiştir.
        kontrol_adi: Karşı grubun adı.
    """
    cikti = CIKTI / klasor
    cikti.mkdir(parents=True, exist_ok=True)
    print(f"=== {ad} ({len(hasta)} gen) / {kontrol_adi} ===")

    # Adım 2: kontrol seçimi
    if karsi is None:
        karsi, denge = eslestir(hasta, k.havuz)
        print(f"Eşleştirme dengesi (SMD, 0.1 altı iyi): uzunluk {denge['uzunluk']:.3f}, GC {denge['GC']:.3f}")
    else:
        denge = {}
        print("Doğrudan karşılaştırma: uzunluk ve GC eşleştirilmedi.")

    # Adım 3: tek CRN dosyası
    paket = cikti / f"{ad}_{kontrol_adi}.crn"
    paketle({ad: hasta, kontrol_adi: karsi}, paket)

    # Adım 4: kodon testi, üç düzey
    dcu = DifferentialCodonUsage.from_crn(paket, ad, kontrol_adi)
    fisher = {}
    for duzey in ("global", "family", "subfamily"):
        fisher[duzey] = sum(s["significant"] for s in dcu.table(duzey, file_name=cikti / f"{duzey}.tsv"))
    print(f"Fisher, anlamlı kodon (q<0.05): küresel {fisher['global']}/61, "
          f"aile {fisher['family']}/59, alt aile {fisher['subfamily']}/59")

    # Adım 5: permütasyon. Aile düzeyinde Fisher'ın bulduklarından hangileri kalıyor?
    saglam = dcu.robust("family", n=PERMUTASYON, seed=TOHUM)
    print(f"Gen düzeyinde permütasyondan sonra kalan aile kodonu: {len(saglam)}/{fisher['family']}")
    en_iyi = sorted(dcu.permutation("family", n=PERMUTASYON, seed=TOHUM), key=lambda s: s["perm_p"])[:3]
    print("  Permütasyonda en güçlü 3 kodon:",
          ", ".join(f"{s['codon']} (p {s['perm_p']:.3f}, q {s['perm_q']:.3f})" for s in en_iyi))
    for s in saglam[:8]:
        print(f"  {s['codon']} ({s['amino_acid']}) pay {s['share_a']:.2f} / {s['share_b']:.2f}, "
              f"OR {s['odds_ratio']:.2f}, Fisher q {s['q']:.1e}, perm q {s['perm_q']:.3f}, lehine {s['favored']}")
    if saglam:
        Hpo.write_tsv(saglam, cikti / "saglam_kodonlar.tsv")

    # Adım 6: CUB parametreleri
    tablo_a, tablo_b = parametreler(dcu.a, k), parametreler(dcu.b, k)
    param = karsilastir(tablo_a, tablo_b)
    Hpo.write_tsv(tablo_a + tablo_b, cikti / "gen_parametreleri.tsv")
    Hpo.write_tsv(param, cikti / "parametre_karsilastirma.tsv")
    print(f"\n{'Parametre':<11}{ad[:9]:>10}{kontrol_adi[:9]:>10}{'delta':>8}{'q':>10}  Etki")
    for s in param:
        print(f"{s['parametre']:<11}{s['medyan_a']:>10.3f}{s['medyan_b']:>10.3f}{s['delta']:>8.2f}"
              f"{s['q']:>10.1e}  {_buyukluk(s['delta'])}{' *' if s['q'] < .05 else ''}")

    # Adım 7: klinik özet
    klinik = klinik_ozet(k, hasta)
    print(f"\nKlinik profil: {klinik['hastalik']} hastalık kaydı")
    print("  Kalıtım:", ", ".join(f"{a} ({n})" for a, n in klinik["kalitim"]))
    print("  Başlangıç:", ", ".join(f"{a} ({n})" for a, n in klinik["baslangic"]))
    print("  Sık fenotipler:", ", ".join(f"{a} %{y}" for a, _, y in klinik["fenotip"][:5]))

    sonuc = Sonuc(ad, kontrol_adi, dcu, hasta, karsi, denge, saglam, fisher, fisher["family"], param,
                  tablo_a, tablo_b, klinik)
    yorum = yorumla(sonuc)
    (cikti / "yorum.txt").write_text(yorum + "\n", encoding="utf-8")
    dcu.plot_volcano("family", file_name=cikti / "yanardag_aile.png")
    print(f"\nYORUM\n{yorum}\n\nÇıktılar: {cikti}\n")
    return sonuc


def yorumla(s: Sonuc) -> str:
    """Sonuçtan kurallı yorum üretir. Yazılan her cümle hesaplanan değere dayanır."""
    satirlar = []
    anlamli = [p for p in s.param if p["q"] < .05]
    if anlamli:
        parca = ", ".join(f"{p['parametre']} ({p['medyan_a']:.3f} / {p['medyan_b']:.3f}, delta {p['delta']:+.2f})"
                          for p in sorted(anlamli, key=lambda x: -abs(x["delta"])))
        satirlar.append(f"Gen düzeyinde {s.ad} ile {s.kontrol} arasında farklı olan parametreler: {parca}.")
    else:
        satirlar.append(f"Gen düzeyinde {s.ad} ile {s.kontrol} arasında anlamlı parametre farkı yok.")
    karisik = [p["parametre"] for p in anlamli if p["parametre"] in ("L", "GC%", "GC3s")]
    if not s.denge and karisik:
        satirlar.append(f"Dikkat: bu karşılaştırma eşleştirilmedi ve {', '.join(karisik)} da farklı. "
                        f"Öteki parametrelerdeki fark büyük olasılıkla uzunluk ve GC'den geliyor, "
                        f"hastalığa özgü kodon seçimi olarak okunmamalıdır.")
    if s.fisher_aile and not s.saglam:
        satirlar.append(f"Fisher {s.fisher_aile} kodon buldu, ama hiçbiri gen düzeyinde permütasyondan "
                        f"geçmedi. Bu fark rastgele gen gruplarında da görülen doğal gen çeşitliliğiyle "
                        f"uyumludur, hastalık grubunun özelliği sayılmamalıdır.")
    elif s.saglam:
        artan = [x["codon"] for x in s.saglam if x["favored"] == s.ad]
        azalan = [x["codon"] for x in s.saglam if x["favored"] != s.ad]
        ucuncu = Counter(c[2] for c in artan)
        satirlar.append(f"Permütasyondan geçen {len(s.saglam)} kodon var. {s.ad} lehine: "
                        f"{', '.join(artan) or 'yok'}; {s.kontrol} lehine: {', '.join(azalan) or 'yok'}.")
        if sum(ucuncu.values()) >= 3:
            gc = ucuncu["G"] + ucuncu["C"]
            satirlar.append(f"{s.ad} lehine kodonların {gc}/{sum(ucuncu.values())} tanesi G/C ile bitiyor; "
                            f"GC3s farkıyla birlikte okunmalı.")
    else:
        satirlar.append("Aile düzeyinde Fisher'ın da permütasyonun da bulduğu kodon yok.")
    if s.denge and max(s.denge.values()) > .1:
        satirlar.append(f"Uyarı: eşleştirme dengesi zayıf (SMD {max(s.denge.values()):.2f}). "
                        f"Kalan fark uzunluk ya da GC'den gelebilir.")
    if anlamli or s.saglam:
        satirlar.append("Bu bulgu ilişkidir, neden değildir. İfade düzeyi eşleştirilmedi; "
                        "klinik karar için kullanılmaz.")
    else:
        satirlar.append(f"Sonuç negatiftir: {len(s.hasta)} genlik bu grupta ortak bir kodon imzası "
                        f"gösterilemedi. Küçük etkiler ya da tek tek genlere özgü farklar kaçmış olabilir.")
    return "\n".join(satirlar)
