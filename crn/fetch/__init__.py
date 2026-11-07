"""Veri kaynakları: NCBI (E-utilities), GENCODE (FTP), HGNC ve HPO; veriler yerel sys/ klasöründe saklanır."""
from crn.fetch.session import NcbiSession
from crn.fetch.ncbi import Ncbi
from crn.fetch.gencode import Gencode
from crn.fetch.hgnc import Hgnc
from crn.fetch.hpo import Hpo

__all__ = ["NcbiSession", "Ncbi", "Gencode", "Hgnc", "Hpo"]
