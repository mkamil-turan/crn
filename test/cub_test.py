from crn.fetch import Ncbi
from crn.cub import CUB

ncbi = Ncbi(organism="human", verbose=False)

# Her gende bir genomik DNA, bir transkript ve bir protein olmalı:
# çok izoformlu gen, izoform genlerine ayrılır.
mefv = ncbi.gene("MEFV")
isoforms = mefv.isoforms()
for isoform in isoforms:
    print(isoform, isoform.features.first("transcript_id"), isoform.features.first("protein_id"))

cub = CUB(isoforms)                # hesap CDS üzerinden; isoforms değişmez
print(cub)

for isoform in cub.genes:
    for tid, rna in isoform.transcripts.items():
        print(tid, len(rna), rna.features.first("orf_start"), rna.features.first("orf_end"))

# A. Ham sayımlar (tek izoform)
name = "MEFV_NM_000243.3"
counts = cub.codon_counts(gene=name)                           # 1. N_kodon
print("1. N(CUG)    :", counts["CUG"], " stop:", counts["UAA"], counts["UAG"], counts["UGA"])
print("2. L         :", cub.total_codons(gene=name))
print("3. Amino asit:", cub.amino_acid_counts(gene=name))
print("4. f(CUG)    :", round(cub.codon_freqs(gene=name)["CUG"], 4))
print("5. ‰(CUG)    :", round(cub.codon_per_thousand(gene=name)["CUG"], 2))

# B. Nükleotit bileşimi
from crn.cub import Composition

comp = Composition(cub)
for key, value in comp.summary(gene=name).items():
    print(f"{key:5}: {value:.4f}")

# C. Sinonim kodon tercihi
from crn.cub import Synonymous

syn = Synonymous(cub)              # referans verilmezse cub'ın kendisi
rscu = syn.rscu(gene=name)
print("14. Kodon oranı (CUG):", round(syn.codon_ratio(gene=name)["CUG"], 4))
print("15. RSCU (Leu)       :", {c: round(rscu[c], 3) for c in cub.codon_families()["L"]})
print("16. Tercih edilen    :", syn.preferred_codons(gene=name))
print("    Az kullanılan    :", syn.underused_codons(gene=name))
print("17. Fop              :", round(syn.fop(gene=name), 4), syn.optimal_codons())

# D. Referans setine dayalı indeksler
from crn.cub import Adaptation


def first_nm(symbol):
    """Genin ilk NM_ izoformu."""
    return [g for g in ncbi.gene(symbol).isoforms()
            if g.features.first("transcript_id").startswith("NM_")][0]


high = CUB([first_nm(s) for s in ("RPL3", "RPL4", "RPS3", "RPS6", "RPL7", "RPS2",
                                         "EEF1A1", "RPL5", "RPL10", "RPS4X", "RPL13", "RPS8")])
low = CUB([first_nm(s) for s in ("OR1A1", "OR2W1", "OR5A1", "OR7A5", "OR1D2", "OR2J2",
                                        "OR4D1", "OR6A2", "OR10A2", "OR51E2", "OR2T4", "OR1E1")])
ada = Adaptation(cub, reference=high, background=low)
w = ada.relative_adaptiveness()
print("18. w (Leu)   :", {c: round(w[c], 3) for c in cub.codon_families()["L"]})
print("19. CAI       :", round(ada.cai(gene=name), 4))
print("20. CBI (rscu):", round(ada.cbi(gene=name), 4))
print("    CBI (chi2):", round(ada.cbi(gene=name, method="chi2"), 4))

# E. Protein özellikleri (E-A … E-G tek sınıfta; PRT ile de çalışır)
from crn.cub import ProteinProperties

pp = ProteinProperties(cub)
pm = pa = pc = ph = po = ps = pd = pp

# E-A. Protein: sayım ve kütle
print("E1. N            :", pm.length(gene=name))
print("E2. Bileşim (%)  :", {aa: round(v, 2) for aa, v in pm.composition(gene=name).items()})
print("E3. Gruplar      :", {g: round(v, 4) for g, v in pm.group_composition(gene=name).items()})
print("E4. MW (Da)      :", round(pm.molecular_weight(gene=name), 2))
print("E5. Monoizotopik :", round(pm.monoisotopic_mass(gene=name), 2))
print("E6. Ort. kalıntı :", round(pm.mean_residue_mass(gene=name), 3))
print("E7. C/W/Y        :", pm.residue_counts(gene=name), " S-S eşi:", pm.disulfide_pairs(gene=name))

# E-B. Protein: atom bileşimi

print("E8.  Formül         :", pa.formula_string(gene=name), pa.formula(gene=name))
print("E9.  Toplam atom    :", pa.atom_count(gene=name))
print("E10. Atom / kalıntı :", round(pa.atoms_per_residue(gene=name), 3))
print("E11. Z_C            :", round(pa.carbon_oxidation_state(gene=name), 4))

# E-C. Protein: yük ve iyonlaşma

print("E12. Yüklü kalıntı :", pc.charged_counts(gene=name))
print("E13. Kaba net yük  :", pc.simple_net_charge(gene=name))
print("E14. FCR           :", round(pc.fcr(gene=name), 4))
print("E15. NCPR          :", round(pc.ncpr(gene=name), 4))
for pka in ("Bjellqvist", "EMBOSS", "Lehninger"):
    print(f"E16/17. {pka:10}: Q(7) = {pc.net_charge(gene=name, ph=7.0, pka=pka):6.2f}  "
          f"pI = {pc.isoelectric_point(gene=name, pka=pka):.2f}")
print("E18. Titrasyon     :", [(ph, round(q, 1)) for ph, q in pc.titration_curve(gene=name, step=2.0)])

# E-D. Protein: hidrofobiklik ve polarite

print("E19. GRAVY          :", round(ph.gravy(gene=name), 4))
for scale in ("Hopp-Woods", "Eisenberg", "Wimley-White", "Janin"):
    print(f"E20. {scale:14}:", round(ph.hydropathy(gene=name, scale=scale), 4))
print("E21. Aromo          :", round(ph.aromaticity(gene=name), 4))
print("E22. Alifatik indeks:", round(ph.aliphatic_index(gene=name), 2))
print("E23. Boman (kcal/mol):", round(ph.boman_index(gene=name), 3))

# E-E. Protein: optik ve hidrodinamik özellikler

print("E24. ε280 (S-S / indirgenmiş):", po.extinction_coefficient(gene=name),
      po.extinction_coefficient(gene=name, reduced=True))
print("E25. Abs %0,1                :", round(po.absorbance_01(gene=name), 3),
      round(po.absorbance_01(gene=name, reduced=True), 3))
print("E26. v̄ (cm³/g)               :", round(po.partial_specific_volume(gene=name), 4))
print("E27. Hacim (Å³) / ASA (Å²)   :", round(po.mean_residue_volume(gene=name), 2),
      round(po.mean_residue_surface(gene=name), 2))

# E-F. Protein: stabilite ve çözünürlük

print("E28. IVYWREL :", round(ps.ivywrel(gene=name), 4))
print("E29. CvP (%) :", round(ps.cvp_bias(gene=name), 2))
print("E30. W–H     :", ps.solubility_probability(gene=name))
print("E31. Entropi :", round(ps.composition_entropy(gene=name), 4))

# E-G. Protein: düzensizlik

print("E32. Uversky :", pd.uversky(gene=name))
print("E33. TOP-IDP :", round(pd.top_idp(gene=name), 4))
print("E34. Özellikler:", {k: round(v, 4) for k, v in pd.composition_features(gene=name).items() if len(k) > 1})

# F. Etkin kodon sayısı ve grafikler (her gen bir nokta: MEFV + yüksek + düşük ifadeli set)
from crn.cub import EffectiveCodons

eff = EffectiveCodons(cub)
print("F23. ENC          :", round(eff.enc(gene=name), 3))
gc3s = Composition(cub).gc3s(gene=name)
print("F24. Beklenen ENC :", round(eff.expected_enc(gc3s), 3), f"(GC3s = {gc3s:.4f})")
print("F25. ENC oranı    :", round(eff.enc_ratio(gene=name), 4))
print("F27. PR2          :", {k: round(v, 4) for k, v in eff.pr2(gene=name).items()})

panel = EffectiveCodons(CUB([first_nm("MEFV")]) + high + low)
neutral = panel.neutrality()
print("F26. Nötrallik    : eğim", round(neutral["slope"], 4), "kesişim", round(neutral["intercept"], 4),
      "r²", round(neutral["r2"], 4), f"({len(neutral['points'])} gen)")

# G. Entropi ve sapma tabanlı ölçümler (RCDI için konak: yüksek + düşük ifadeli set)
from crn.cub import Divergence

div = Divergence(cub, reference=high + low)
print("G28. SCUO :", round(div.scuo(gene=name), 4))
print("G29. CDC  :", round(div.cdc(gene=name), 4))
print("G30. RCDI :", round(div.rcdi(gene=name), 4))
print("G31. RCBS :", round(div.rcbs(gene=name), 4))

# H. tRNA ve kodon bağlamı (CPS için referans: yüksek + düşük ifadeli set)
from crn.cub import CodonContext

ctx = CodonContext(cub, reference=high + low)
rho = ctx.rho(gene=name)
print("H33. ρ(CpG)       :", round(rho["CG"], 4), " ρ(UpA):", round(rho["UA"], 4),
      " ρ(CpG, 3-1):", round(ctx.rho(gene=name, frame="31")["CG"], 4))
print("H34. CPB          :", round(ctx.cpb(gene=name), 4))
# H32. tAI: tRNA gen kopya sayıları GtRNAdb'den (gtrnadb.ucsc.edu, Hsapi38)
from crn.cub import read_trna_counts

trna = read_trna_counts("sys/gtrnadb/hg38-tRNAs.fa")
print("H32. tAI          :", round(ctx.tai(gene=name, trna_counts=trna), 4),
      f"({sum(trna.values())} tRNA geni, {len(trna)} antikodon)")

# Tüm genler, tüm skaler ölçümler tek tabloda (TSV olarak da yazılabilir)
rows = cub.table(reference=high, background=low, host=high + low, genome=high + low,
                 trna_counts=trna)          # file_name="tablo.tsv" ile dosyaya da yazılır
print("Tablo  :", len(rows), "satır,", len(rows[0]), "sütun")
