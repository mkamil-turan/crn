"""Veri kaynakları: NCBI (E-utilities) ve GENCODE (FTP); veriler yerel sys/ klasöründe saklanır."""
from crn.fetch.session import NcbiSession
from crn.fetch.ncbi import Ncbi
from crn.fetch.gencode import Gencode

__all__ = ["NcbiSession", "Ncbi", "Gencode"]
