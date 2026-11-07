"""Örnek 11: Gerçek referans üzerinde klinik varyantları doğrula.

Senaryo: Ailesel Akdeniz Ateşi panelinizdeki MEFV varyantlarının NCBI referans
transkriptinde beklenen baz ve amino asit değişimini verdiğini kontrol etmek ve
mutant proteinleri saklamak istiyorsunuz.

Bu örnekte:
  1. NM_000243 transkriptinden CDS çıkarılır.
  2. Her varyantta referans baz doğrulanır.
  3. Mutant CDS proteine çevrilir, değişim bulunup beklenenle karşılaştırılır.
  4. Referans ve mutant proteinler (PRT) tek CRN listesine yazılır.

Çıktılar: test/ncbi/cikti/ncbi_orn_11/
Not: Sonuç yalnızca protein değişimini doğrular, patojenite bilgisi vermez.
"""
from pathlib import Path

from crn import write_crn_list, write_fasta
from crn.fetch import Ncbi

CIKTI = Path(__file__).parent / "cikti" / "ncbi_orn_11"
CIKTI.mkdir(parents=True, exist_ok=True)

UC_HARF = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
           "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
           "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}

# (cDNA konumu, referans, alternatif, beklenen protein değişimi)
VARYANTLAR = [(2080, "A", "G", "p.Met694Val"), (2082, "G", "A", "p.Met694Ile"),
              (2040, "G", "C", "p.Met680Ile"), (2177, "T", "C", "p.Val726Ala"),
              (442, "G", "C", "p.Glu148Gln")]

ncbi = Ncbi(verbose=False)
rna = ncbi.rna("NM_000243")
bas, bit = rna.features.first("cds")
cds = rna[bas - 1:bit]
referans = cds.translation()
referans.name = "MEFV_referans"
print(f"{rna.name}: CDS {len(cds)} nt, protein {len(referans) - 1} aa\n")

proteinler = [referans]
print(f"{'cDNA':<12}{'Kodon':<12}{'Protein':<14}{'Beklenen':<14}Sonuç")
for konum, ref, alt, beklenen in VARYANTLAR:
    i = konum - 1
    ref_rna, alt_rna = ref.replace("T", "U"), alt.replace("T", "U")
    if cds[i] != ref_rna:
        raise ValueError(f"c.{konum}: {ref_rna} beklendi, {cds[i]} bulundu; transkript sürümünü denetleyin")

    mutant = cds[:]
    mutant.seq = cds.seq[:i] + alt_rna + cds.seq[i + 1:]
    protein = mutant.translation()

    fark = [j for j, (a, b) in enumerate(zip(referans.seq, protein.seq)) if a != b]
    j = fark[0]
    gosterim = f"p.{UC_HARF[referans[j]]}{j + 1}{UC_HARF[protein[j]]}"
    protein.name = f"MEFV_{gosterim}"
    proteinler.append(protein)

    kodon_i = i // 3 * 3
    kodon = f"{cds.seq[kodon_i:kodon_i + 3]}>{mutant.seq[kodon_i:kodon_i + 3]}"
    print(f"{f'c.{konum}{ref}>{alt}':<12}{kodon:<12}{gosterim:<14}{beklenen:<14}"
          f"{'uyumlu' if gosterim == beklenen else 'UYUMSUZ'}")

write_crn_list(proteinler, CIKTI / "mefv_varyant_proteinleri.crn")
write_fasta(proteinler, CIKTI / "mefv_varyant_proteinleri.fasta")
print("\nÇıktılar:", CIKTI)
