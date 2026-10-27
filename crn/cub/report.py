"""Tüm genler için skaler CUB ve protein ölçümlerini tek tabloda toplar."""
from __future__ import annotations

import csv
import warnings
from pathlib import Path
from typing import Dict, List, Mapping, Optional

from crn.cub.core import CUB, CubWarning
from crn.cub.composition import Composition
from crn.cub.adaptation import Adaptation
from crn.cub.enc import EffectiveCodons
from crn.cub.divergence import Divergence
from crn.cub.context import CodonContext
from crn.protein import ProteinProperties


def summary_table(cub: CUB, reference: Optional[CUB] = None, background: Optional[CUB] = None,
                  host: Optional[CUB] = None, genome: Optional[CUB] = None,
                  trna_counts: Optional[Mapping[str, float]] = None, kingdom: str = "eukaryote",
                  protein: bool = True, file_name: Optional[str | Path] = None
                  ) -> List[Dict[str, object]]:
    """Her gen için bir satır: ``{"gene": ad, ölçüm: değer, ...}``.

    Her zaman: L, GC%, GC1, GC2, GC12, GC3, GC3s, ENC, ENC_expected, ENC_ratio,
    SCUO, CDC, RCBS. İsteğe bağlı sütunlar ilgili veri verilince eklenir:

    - ``reference`` (yüksek ifadeli genler): Fop, CAI, CBI; ``background``
      verilirse Fop_chi2 ve CBI_chi2 de (optimal kodonlar χ² ile).
    - ``host``: RCDI.  ``genome``: CPB.  ``trna_counts``: tAI.
    - ``protein`` True ise ``ProteinProperties.summary`` sütunları.

    Satır satır ``CubWarning`` uyarıları (RCDI/CPB'de referansta olmayan
    kodon/çiftler) bastırılır. ``file_name`` verilirse tablo sekmeyle ayrılmış
    (TSV) olarak da yazılır.
    """
    comp, eff, div = Composition(cub), EffectiveCodons(cub), Divergence(cub)
    ada = Adaptation(cub, reference, background) if reference is not None else None
    rcdi = Divergence(cub, host) if host is not None else None
    cpb = CodonContext(cub, genome) if genome is not None else None
    tai = CodonContext(cub) if trna_counts else None
    props = ProteinProperties(cub) if protein else None

    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", CubWarning)
        for name in cub.names:
            gc3s = comp.gc3s(name)
            enc = eff.enc(name)
            expected = eff.expected_enc(gc3s)
            row: Dict[str, object] = {"gene": name, "L": cub.total_codons(name)}
            row.update({key: value for key, value in comp.summary(name).items()
                        if key in ("GC%", "GC1", "GC2", "GC12", "GC3", "GC3s")})
            row.update({"ENC": enc, "ENC_expected": expected,
                        "ENC_ratio": (expected - enc) / expected,
                        "SCUO": div.scuo(name), "CDC": div.cdc(name), "RCBS": div.rcbs(name)})
            if ada is not None:
                row.update({"Fop": ada.synonymous.fop(name), "CAI": ada.cai(name), "CBI": ada.cbi(name)})
                if background is not None:
                    row.update({"Fop_chi2": ada.synonymous.fop(name, method="chi2"),
                                "CBI_chi2": ada.cbi(name, method="chi2")})
            if rcdi is not None:
                row["RCDI"] = rcdi.rcdi(name)
            if cpb is not None:
                row["CPB"] = cpb.cpb(name)
            if tai is not None:
                row["tAI"] = tai.tai(name, trna_counts, kingdom=kingdom)
            if props is not None:
                row.update(props.summary(name))
            rows.append(row)

    if file_name is not None:
        with open(file_name, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    return rows
