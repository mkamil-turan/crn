from __future__ import annotations

import copy
from typing import Any, Dict, Iterator, List, Optional, Tuple


class Features:
    """Bir diziye ya da gene ait ek bilgileri anahtar -> değer listesi olarak tutar.

    Kaynak, gen, CDS koordinatı gibi bilgiler sınıflara yeni alan olarak
    değil, buraya eklenir. Her anahtar bir değer listesi tutar: ``add`` listeye
    ekler, ``set`` listeyi tek değerle değiştirir, ``first`` ilk değeri verir.
    Anahtarlar eklenme sırasını korur; dışarı verilen listeler kopyadır.

    Örnek::

        f = Features()
        f.add("gene", "BRCA1")
        f.add("gene", "BRCA2")
        f.get("gene")      # ['BRCA1', 'BRCA2']
        f.first("gene")    # 'BRCA1'
        "gene" in f        # True
    """

    __slots__ = ("_data",)

    def __init__(self):
        self._data: Dict[str, List[Any]] = {}

    @staticmethod
    def _check_key(key) -> None:
        if not isinstance(key, str):
            raise TypeError(f"Anahtar str olmalı, {type(key).__name__} verildi")
        if not key:
            raise ValueError("Anahtar boş olamaz")

    def add(self, key: str, value: Any) -> None:
        """``key`` anahtarının değer listesinin sonuna ``value`` ekler."""
        self._check_key(key)
        self._data.setdefault(key, []).append(value)

    def set(self, key: str, value: Any) -> None:
        """``key`` anahtarının tüm değerlerini silip yerine tek ``value`` koyar."""
        self._check_key(key)
        self._data[key] = [value]

    def get(self, key: str, default: Any = None) -> Optional[List[Any]]:
        """Değer listesinin kopyası; anahtar yoksa ``default``."""
        return list(self._data[key]) if key in self._data else default

    def first(self, key: str, default: Any = None) -> Any:
        """İlk değer; anahtar yoksa ``default``."""
        values = self._data.get(key)
        return values[0] if values else default

    def remove(self, key: str) -> bool:
        """Anahtarı tüm değerleriyle siler; anahtar varsa ``True`` döner."""
        return self._data.pop(key, None) is not None

    def clear(self) -> None:
        self._data = {}

    def keys(self) -> List[str]:
        return list(self._data)

    def items(self) -> Iterator[Tuple[str, List[Any]]]:
        for key, values in self._data.items():
            yield key, list(values)

    def to_dict(self) -> Dict[str, List[Any]]:
        """``{anahtar: değer listesi}``; sözlük ve listeler kopyadır."""
        return {key: list(values) for key, values in self._data.items()}

    def copy(self) -> Features:
        """Bağımsız (derin) kopya."""
        result = Features()
        result._data = copy.deepcopy(self._data)
        return result

    def __contains__(self, key) -> bool:
        return key in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self) -> Iterator[str]:
        return iter(list(self._data))

    def __repr__(self) -> str:
        return f"Keys:{','.join(self._data)}"
