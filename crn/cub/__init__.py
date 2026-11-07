"""Kodon kullanım yanlılığı (A–H, J) ve protein özellikleri (E).

==========  =====================  ==================================================
 Bölüm      Sınıf                  Ölçümler
==========  =====================  ==================================================
 A          ``CUB``                kodon/amino asit sayıları, frekanslar
 B          ``Composition``        GC%, GC1-3, GC12, GC3s, A3s/U3s/G3s/C3s
 C          ``Synonymous``         kodon oranı, RSCU (+ çubuk/ısı haritası), Fop
 D          ``Adaptation``         w, CAI (+ dağılım), CBI
 E          ``ProteinProperties``  kütle, formül, yük/pI, hidropati, optik, stabilite
 F          ``EffectiveCodons``    ENC, beklenen ENC, nötrallik, PR2, grafikler
 G          ``Divergence``         SCUO, CDC, RCDI, RCBS
 H          ``CodonContext``       tAI, dinükleotit ρ, CPS/CPB
 J          ``DifferentialCodonUsage``  iki gen grubunun kodon kodon karşılaştırması,
                                   gen düzeyinde permütasyon (``robust``)
 —          ``CorrespondenceAnalysis``  RSCU'nun uyum çözümlemesi (COA) + saçılım
==========  =====================  ==================================================

``CUB.table()`` tüm skaler ölçümleri gen başına tek tabloda verir;
``plot_correlations`` bu indekslerin korelasyon ısı haritasını çizer.

Grafik metotları: ``Synonymous.plot_rscu``, ``plot_rscu_heatmap``;
``Adaptation.plot_cai``; ``EffectiveCodons.plot_enc``, ``plot_enc_ratio``,
``plot_neutrality``, ``plot_pr2``; ``CorrespondenceAnalysis.plot_coa``;
``plot_correlations``; ``DifferentialCodonUsage.plot_volcano``. Hepsi ``ax=``/``file_name=`` alır ve ``Axes`` döner.
"""
from crn.cub.core import CUB, CubWarning
from crn.cub.composition import Composition
from crn.cub.synonymous import Synonymous
from crn.cub.adaptation import Adaptation
from crn.cub.enc import EffectiveCodons
from crn.cub.divergence import Divergence
from crn.cub.context import CodonContext, read_trna_counts
from crn.cub.ordination import CorrespondenceAnalysis
from crn.cub.differential import DifferentialCodonUsage
from crn.cub.report import summary_table, correlation_matrix, plot_correlations
from crn.protein import ProteinProperties

__all__ = ["CUB", "CubWarning", "Composition", "Synonymous", "Adaptation", "ProteinProperties",
           "EffectiveCodons", "Divergence", "CodonContext", "read_trna_counts",
           "CorrespondenceAnalysis", "DifferentialCodonUsage", "summary_table", "correlation_matrix", "plot_correlations"]
