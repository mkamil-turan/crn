from __future__ import annotations

import copy
from typing import (Any, Dict, Optional, Tuple, List, Iterator)

class Features:
    """Bir diziye ait tüm ek bilgileri anahtar–değer listesi olarak tutar.

    ``Sequence``, ``DNA``, ``RNA`` ve ``PRT`` sınıflarına eklenecek her türlü
    bilgi (kaynak, gen, lokus tipi, versiyon, strand vb.) sınıflara yeni bir
    property olarak değil, bu sınıf üzerinden eklenir.

    Her anahtar bir değer listesi tutar; aynı anahtara ``add`` ile birden
    fazla değer eklenebilir. Tek değerli bilgiler için ``set`` ve ``first``
    kullanılır. Anahtarlar eklenme sırasını korur.

    Örnek::

        f = Features()
        f.add("source", "NCBI")
        f.add("gene", "BRCA1")
        f.add("gene", "BRCA2")
        f.get("gene")      # ['BRCA1', 'BRCA2']
        f.first("source")  # 'NCBI'
        "gene" in f        # True
    """

    def __init__(self):
        """Boş bir özellik deposu oluşturur."""
        self.__data: Dict[str, List[Any]] = {}

    @staticmethod
    def __check_key(key) -> None:
        """Anahtarın boş olmayan bir ``str`` olduğunu doğrular.

        Raises:
            TypeError: Anahtar ``str`` değilse.
            ValueError: Anahtar boşsa.
        """
        if not isinstance(key, str):
            raise TypeError(f"Anahtar str olmalı, {type(key).__name__} verildi")
        if not key:
            raise ValueError("Anahtar boş olamaz")

    def add(self, key: str, value: Any) -> None:
        """``key`` anahtarının değer listesinin sonuna ``value`` ekler.

        Anahtar yoksa tek elemanlı yeni bir liste ile oluşturulur.

        Raises:
            TypeError: Anahtar ``str`` değilse.
            ValueError: Anahtar boşsa.
        """
        self.__check_key(key)
        if key in self.__data:
            self.__data[key].append(value)
        else:
            self.__data[key] = [value, ]

    def set(self, key: str, value: Any) -> None:
        """``key`` anahtarının tüm değerlerini silip yerine tek ``value`` koyar.

        Versiyon, strand gibi tek değer taşıması gereken bilgiler için kullanılır.

        Raises:
            TypeError: Anahtar ``str`` değilse.
            ValueError: Anahtar boşsa.
        """
        self.__check_key(key)
        self.__data[key] = [value, ]

    def get(self, key: str, default: Any = None) -> Optional[List[Any]]:
        """``key`` anahtarının değer listesini döndürür, anahtar yoksa ``default``.

        Dönen liste bir kopyadır; üzerinde yapılan değişiklikler depoyu etkilemez.
        """
        if key in self.__data:
            return list(self.__data[key])
        else:
            return default

    def first(self, key: str, default: Any = None) -> Any:
        """``key`` anahtarının ilk değerini döndürür, anahtar yoksa ``default``."""
        values = self.__data.get(key)
        return values[0] if values else default

    def keys(self) -> List[str]:
        """Anahtarları eklenme sırasıyla liste olarak döndürür."""
        return list(self.__data.keys())

    def items(self) -> Iterator[Tuple[str, List[Any]]]:
        """(anahtar, değer listesi) çiftlerini sırayla döndürür.

        Değer listeleri kopyadır; üzerinde yapılan değişiklikler depoyu etkilemez.
        """
        for key, values in self.__data.items():
            yield key, list(values)

    def to_dict(self) -> Dict[str, List[Any]]:
        """Tüm özellikleri ``{anahtar: değer listesi}`` sözlüğü olarak döndürür.

        Dönen sözlük ve listeler kopyadır; kaydetme/dışa aktarma için uygundur.
        """
        return {key: list(values) for key, values in self.__data.items()}

    def copy(self) -> Features:
        """Bağımsız bir kopya döndürür; kopyada yapılan değişiklik bu nesneyi etkilemez.

        Liste gibi değiştirilebilir değerler de kopyalanır (derin kopya).
        """
        result = Features()
        result.__data = copy.deepcopy(self.__data)
        return result

    def clear(self) -> None:
        """Tüm özellikleri siler."""
        self.__data = {}

    def remove(self, key: str) -> bool:
        """``key`` anahtarını tüm değerleriyle siler.

        Returns:
            Anahtar varsa ve silindiyse ``True``, yoksa ``False``.
        """
        if key in self.__data:
            self.__data.pop(key, None)
            return True
        else:
            return False

    def __contains__(self, key) -> bool:
        """``"gene" in features`` biçiminde anahtar kontrolü sağlar."""
        return key in self.__data

    def __len__(self) -> int:
        """Anahtar sayısını döndürür."""
        return len(self.__data)

    def __iter__(self) -> Iterator[str]:
        """Anahtarlar üzerinde eklenme sırasıyla dolaşır."""
        return iter(list(self.__data.keys()))

    def __repr__(self):
        """Anahtarları virgülle ayrılmış olarak gösterir: ``Keys:source,gene``."""
        return f"Keys:{",".join(self.__data.keys())}"
