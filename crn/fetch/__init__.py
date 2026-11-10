"""Veri kaynakları: NCBI (E-utilities), GENCODE (FTP), HGNC, HPO ve ClinVar; veriler yerel sys/ klasöründe saklanır."""
from crn.fetch.session import NcbiSession
from crn.fetch.ncbi import Ncbi
from crn.fetch.gencode import Gencode
from crn.fetch.hgnc import Hgnc
from crn.fetch.hpo import Hpo
from crn.fetch.clinvar import Clinvar

__all__ = ["NcbiSession", "Ncbi", "Gencode", "Hgnc", "Hpo", "Clinvar"]
