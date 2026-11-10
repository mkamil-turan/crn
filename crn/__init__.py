"""crn: DNA, RNA ve protein dizileri, NCBI/GENCODE/HPO erişimi, kodon kullanım yanlılığı.

::

    crn
    ├── constants  alfabeler, kodonlar, tümleme tabloları, GeneticCode
    ├── tables     kütleler, pKa, hidropati ölçekleri, tAI cezaları ...
    ├── features   Features (anahtar -> değer listesi)
    ├── seq        Sequence, DNA, RNA, PRT, Strand, GenomicStrand
    ├── gene       Gene (DNA + transkriptler + proteinler)
    ├── fileio     FASTA ve CRN dosyaları
    ├── protein    ProteinProperties (E-A … E-G)
    ├── fetch      Ncbi, Gencode, NcbiSession, Hgnc, Hpo
    ├── cub        CUB ve A–H bölümleri
    └── cli        komut satırı arayüzü (crn komutu)
"""
from crn.constants import GeneticCode
from crn.features import Features
from crn.seq import Molecule, Strand, GenomicStrand, Sequence, NucleicAcid, DNA, RNA, PRT
from crn.gene import Gene
from crn.fileio import read_fasta, write_fasta, read_crn_list, write_crn_list
from crn.protein import ProteinProperties

__all__ = ["GeneticCode", "Features", "Molecule", "Strand", "GenomicStrand", "Sequence",
           "NucleicAcid", "DNA", "RNA", "PRT", "Gene", "read_fasta", "write_fasta",
           "read_crn_list", "write_crn_list", "ProteinProperties"]
