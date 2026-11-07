"""Veri kaynaklarının (Ncbi, Gencode, Hgnc, Hpo) ortak yerel depo ve ilerleme yardımcıları."""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Iterator, Optional, Tuple

from tqdm import tqdm

# Verilerin tutulduğu klasör: CRN_SYS_DIR ortam değişkeni, yoksa proje kökündeki sys/.
DEFAULT_SYS_DIR = Path(os.environ.get("CRN_SYS_DIR") or Path(__file__).resolve().parents[2] / "sys")


def atomic_write(path: Path, text: str) -> None:
    """Dosyayı önce .tmp olarak yazıp sonra yerine koyar; yarım dosya kalmaz."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    tmp.replace(path)


def split_fasta(text: str) -> Iterator[Tuple[str, str]]:
    """Çok kayıtlı FASTA metnini (ilk kelime, kaydın metni) çiftlerine böler."""
    for block in re.split(r"(?m)^(?=>)", text):
        if block.startswith(">"):
            yield block[1:].split(None, 1)[0], block.rstrip("\n") + "\n"


def read_single_fasta(path: Path) -> Tuple[str, str]:
    """Tek kayıtlı FASTA dosyasından (başlık, dizi)."""
    lines = path.read_text().splitlines()
    return lines[0][1:].strip(), "".join(line.strip() for line in lines[1:])


class LocalSource:
    """Yerel ``sys/`` klasörüyle çalışan veri kaynağının ortak tabanı.

    Alt sınıf ``SOURCE`` adını verir; çıkarılan her nesnenin features'ına
    ``source`` olarak yazılır.
    """

    SOURCE = ""

    def __init__(self, sys_dir: Optional[str | Path] = None, verbose: bool = True):
        self.sys_dir = Path(sys_dir) if sys_dir is not None else DEFAULT_SYS_DIR
        self.verbose = verbose

    def _log(self, message: str) -> None:
        if self.verbose:
            tqdm.write(message)

    def _bar(self, total: Optional[int], desc: str, unit: str = "B") -> tqdm:
        """İlerleme çubuğu; ``unit="B"`` ise bayt ölçeğiyle. verbose False ise gizli."""
        in_bytes = unit == "B"
        return tqdm(total=total or None, desc=desc, unit=unit, unit_scale=in_bytes,
                    unit_divisor=1024 if in_bytes else 1000, dynamic_ncols=True,
                    disable=not self.verbose)

    def _progress(self, items, desc: str, unit: str):
        return tqdm(items, desc=desc, unit=unit, dynamic_ncols=True, disable=not self.verbose)

    def _features(self, obj, empty=(None, ""), **values) -> None:
        """``source`` ve boş olmayan değerleri obj.features'a yazar."""
        obj.features.set("source", self.SOURCE)
        for key, value in values.items():
            if value not in empty:
                obj.features.set(key, value)
