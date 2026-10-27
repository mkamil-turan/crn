from __future__ import annotations

import gzip
import struct
import zlib
from pathlib import Path
from typing import (Any, Dict, FrozenSet, Iterator, Optional, Self)

from crn.io.crn_io_type import (Molecule, Strand, GenomicStrand)
from crn.io.crn_io_features import Features
from crn.io.crn_io_const import (Maps, GeneticCode)


# CRN ikili dosya biçimi (tüm sayılar little-endian):
#
#   magic      4 bayt   b"CRN\x00"
#   version    uint16   CRN_VERSION
#   type       str      Molecule üyesinin adı ("SEQ", "DNA" ...)
#   name       str
#   seq        str
#   n_root     uint32
#   her root alanı için (ör. DNA'nın strand'i):
#       key      str
#       value    etiketli değer
#   n_features uint32
#   her özellik için:
#       key      str
#       n_values uint32
#       values   n_values adet etiketli değer
#   crc32      uint32   magic dahil önceki tüm baytların CRC32'si
#
# "type"tan "features"a kadar olan kısım bir dizi kaydıdır (_crn_record).
#
# Gen dosyası (Gene.to_crn) aynı magic, version ve crc32'yi kullanır:
#
#   magic, version
#   type          str      "GENE"
#   name          str
#   n_features    uint32 + özellikler (yukarıdaki gibi)
#   has_dna       uint8    1 ise ardından genomik DNA'nın dizi kaydı
#   n_transcripts uint32
#   her transkript için:
#       id          str
#       rna         dizi kaydı
#   n_proteins    uint32
#   her protein için:
#       id          str
#       prt         dizi kaydı
#   crc32
#
#   str: uint32 bayt uzunluğu + UTF-8 baytları
#
# Etiketli değer: 1 baytlık etiket + içerik
#   b"N" None
#   b"T" True, b"F" False
#   b"i" int    -> uint8 uzunluk + işaretli little-endian baytlar
#   b"f" float  -> float64
#   b"s" str
#   b"b" bytes  -> uint32 uzunluk + baytlar
#   b"e" enum   -> enum sınıfının adı (str) + üye adı (str)
#                  yalnızca CRN_ENUMS içindeki enum'lar yazılır/okunur
#   b"t" tuple, b"l" list -> uint32 eleman sayısı + etiketli değerler
CRN_MAGIC = b"CRN\x00"
CRN_VERSION = 3
CRN_ENUMS = {e.__name__: e for e in (Molecule, Strand, GenomicStrand, Maps, GeneticCode)}

# Molecule -> o molekülü temsil eden sınıf. Alt sınıflar tanımlandıkça dolar.
_SEQUENCE_CLASSES: Dict[Molecule, type] = {}


def _crn_pack_str(text: str) -> bytes:
    data = text.encode("utf-8")
    return struct.pack("<I", len(data)) + data


def _crn_pack_value(value: Any) -> bytes:
    if value is None:
        return b"N"
    if isinstance(value, bool):
        return b"T" if value else b"F"
    if type(value) in CRN_ENUMS.values():
        return b"e" + _crn_pack_str(type(value).__name__) + _crn_pack_str(value.name)
    if isinstance(value, int):
        data = value.to_bytes((value.bit_length() + 8) // 8, "little", signed=True)
        return b"i" + struct.pack("<B", len(data)) + data
    if isinstance(value, float):
        return b"f" + struct.pack("<d", value)
    if isinstance(value, str):
        return b"s" + _crn_pack_str(value)
    if isinstance(value, bytes):
        return b"b" + struct.pack("<I", len(value)) + value
    if isinstance(value, (tuple, list)):
        tag = b"t" if isinstance(value, tuple) else b"l"
        return tag + struct.pack("<I", len(value)) + b"".join(_crn_pack_value(v) for v in value)
    raise TypeError(f"CRN dosyasına yazılamayan değer tipi: {type(value).__name__}")


def _crn_pack_features(features: Features) -> bytes:
    parts = [struct.pack("<I", len(features))]
    for key, values in features.items():
        parts.append(_crn_pack_str(key))
        parts.append(struct.pack("<I", len(values)))
        parts.extend(_crn_pack_value(v) for v in values)
    return b"".join(parts)


def _crn_write(file_name: str | Path, body: bytes) -> None:
    """magic + version + body + crc32 olarak dosyaya yazar."""
    data = CRN_MAGIC + struct.pack("<H", CRN_VERSION) + body
    data += struct.pack("<I", zlib.crc32(data))
    with open(file_name, "wb") as f:
        f.write(data)


def _crn_open(file_name: str | Path) -> _CrnReader:
    """Dosyayı okur; magic, CRC32 ve sürümü doğrular, okuyucuyu sürümün
    hemen sonrasına konumlanmış olarak döndürür.

    Raises:
        FileNotFoundError: Dosya yoksa.
        ValueError: Dosya CRN değilse, bozuksa ya da sürümü desteklenmiyorsa.
    """
    path = Path(file_name)
    data = path.read_bytes()

    if len(data) < len(CRN_MAGIC) + 6 or not data.startswith(CRN_MAGIC):
        raise ValueError(f"{path}: CRN dosyası değil")
    body, crc = data[:-4], struct.unpack("<I", data[-4:])[0]
    if zlib.crc32(body) != crc:
        raise ValueError(f"{path}: CRN dosyası bozuk (CRC32 uyuşmuyor)")

    reader = _CrnReader(body, str(path))
    reader.take(len(CRN_MAGIC))
    version = reader.uint("<H")
    if version != CRN_VERSION:
        raise ValueError(f"{path}: desteklenmeyen CRN sürümü: {version}")
    return reader


class _CrnReader:
    """CRN baytlarını baştan sona sırayla okur; eksik veride ValueError verir."""

    def __init__(self, data: bytes, source: str):
        self.data = data
        self.source = source
        self.pos = 0

    def take(self, n: int) -> bytes:
        if self.pos + n > len(self.data):
            raise ValueError(f"{self.source}: CRN dosyası beklenenden kısa")
        chunk = self.data[self.pos:self.pos + n]
        self.pos += n
        return chunk

    def uint(self, fmt: str) -> int:
        return struct.unpack(fmt, self.take(struct.calcsize(fmt)))[0]

    def text(self) -> str:
        return self.take(self.uint("<I")).decode("utf-8")

    def features(self) -> Features:
        """``_crn_pack_features`` ile yazılmış özellikleri okur."""
        result = Features()
        for _ in range(self.uint("<I")):
            key = self.text()
            for _ in range(self.uint("<I")):
                result.add(key, self.value())
        return result

    def finish(self) -> None:
        """Verinin tamamen okunduğunu doğrular."""
        if self.pos != len(self.data):
            raise ValueError(f"{self.source}: CRN dosyasının sonunda fazladan veri var")

    def value(self) -> Any:
        tag = self.take(1)
        if tag == b"N":
            return None
        if tag == b"T":
            return True
        if tag == b"F":
            return False
        if tag == b"i":
            return int.from_bytes(self.take(self.uint("<B")), "little", signed=True)
        if tag == b"f":
            return self.uint("<d")
        if tag == b"s":
            return self.text()
        if tag == b"b":
            return self.take(self.uint("<I"))
        if tag == b"e":
            enum_name = self.text()
            enum_cls = CRN_ENUMS.get(enum_name)
            if enum_cls is None:
                raise ValueError(f"{self.source}: bilinmeyen enum: {enum_name}")
            name = self.text()
            if name not in enum_cls.__members__:
                raise ValueError(f"{self.source}: bilinmeyen {enum_name} üyesi: {name}")
            return enum_cls[name]
        if tag in (b"t", b"l"):
            items = [self.value() for _ in range(self.uint("<I"))]
            return tuple(items) if tag == b"t" else items
        raise ValueError(f"{self.source}: bilinmeyen değer etiketi: {tag!r}")


class Sequence:

    # Alt sınıflar kendi molekül tipini ve izin verilen harfleri burada belirtir.
    # _alphabet None ise harf kontrolü yapılmaz.
    _molecule: Molecule = Molecule.SEQ
    _alphabet: Optional[FrozenSet[str]] = None

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if "_molecule" in cls.__dict__:
            _SEQUENCE_CLASSES[cls._molecule] = cls

    def __init__(self, seq:str,
                 name: str = ""):
        self._type = self._molecule
        self.__seq = self._check_seq(seq)
        self.__name = name
        self.__features = Features()

    def _check_seq(self, seq: str) -> str:
        """Diziyi büyük harfe çevirir; ``_alphabet`` tanımlıysa harfleri doğrular."""
        seq = seq.upper()
        if self._alphabet is not None:
            invalid = sorted(set(seq) - self._alphabet)
            if invalid:
                raise ValueError(
                    f"{self._molecule.value[1]} dizisinde geçersiz karakter(ler): {''.join(invalid)}"
                )
        return seq

    def _root_fields(self) -> Dict[str, Any]:
        """features dışında, sınıfın kendisine ait (root) alanlar.

        Alt sınıflar ``__init__`` parametrelerini burada döndürür (ör. DNA'nın
        strand'i). Bu alanlar CRN dosyasına yazılır, kopyalarda korunur ve
        ``cls(seq, name=..., **root_fields)`` ile geri verilir.
        """
        return {}

    def _new(self, seq: str) -> Self:
        """Aynı sınıf, ad, root alanları ve features kopyasıyla yeni bir nesne."""
        result = self.__class__(seq, name=self.__name, **self._root_fields())
        result.features = self.__features.copy()
        return result

    @classmethod
    def from_fasta_file(cls, file_name: str | Path) -> list[Self]:
        """FASTA dosyasını okur ve dizileri ``cls`` tipinde bir liste olarak döndürür.

        Dosyada tek kayıt olsa bile liste döner.

        - ``.gz`` uzantılı dosyalar gzip ile açılır.
        - Dizinin adı başlığın ilk kelimesidir.
        - Boş satırlar ve ``;`` ile başlayan yorum satırları atlanır.
        - Dizi satırlarındaki boşluklar silinir; çok satırlı diziler birleştirilir.

        Her dizinin features'ına şunlar yazılır:

        - ``fasta_header``: ``>`` işaretinden sonraki başlık satırının tamamı
        - ``source_file``: okunan dosyanın yolu

        Raises:
            FileNotFoundError: Dosya yoksa.
            ValueError: Başlıksız dizi verisi, boş başlık, dizisi olmayan başlık,
                ``cls`` için geçersiz karakter ya da hiç kayıt yoksa. Mesaj dosya
                adını ve satır numarasını içerir.
        """
        path = Path(file_name)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8-sig") as f:
            lines = f.read().splitlines()

        sequences: list[Self] = []
        header: str | None = None
        header_line = 0
        seq_chunks: list[str] = []

        def flush() -> None:
            if not seq_chunks:
                raise ValueError(f"{path}:{header_line}: '>{header}' bir seq bloğu olmadan bitti")
            try:
                seq = cls("".join(seq_chunks), name=header.split()[0])
            except ValueError as e:
                raise ValueError(f"{path}:{header_line}: {e}") from None
            seq.features.set("fasta_header", header)
            seq.features.set("source_file", str(path))
            sequences.append(seq)

        for line_no, line in enumerate(lines, start=1):
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if header is not None:
                    flush()
                header = line[1:].strip()
                header_line = line_no
                if not header:
                    raise ValueError(f"{path}:{line_no}: boş FASTA başlığı")
                seq_chunks = []
            else:
                if header is None:
                    raise ValueError(f"{path}:{line_no}: '>' header'ı olmadan sequence verisi bulundu")
                seq_chunks.append("".join(line.split()))

        if header is not None:
            flush()

        if not sequences:
            raise ValueError(f"{path}: geçerli bir FASTA kaydı bulunamadı")

        return sequences

    def to_fasta(self, file_name: str | Path) -> None:
        """Diziyi FASTA formatında ``file_name`` dosyasına yazar.

        Başlık satırı ``>NAME this sequence contains N feature(s) but not listed``
        biçimindedir. ``N``, ``from_fasta_file``'ın eklediği ``fasta_header`` ve
        ``source_file`` dışındaki özelliklerin sayısıdır. Dizi 60 karakterlik
        satırlara bölünür. Dosya varsa üzerine yazılır.

        Raises:
            ValueError: Ad boşsa ya da boşluk içeriyorsa, veya dizi boşsa. Bu
                durumlarda dosya ``from_fasta_file`` ile aynı haliyle geri
                okunamayacağı için hiçbir şey yazılmaz.
        """
        if not self.__name or len(self.__name.split()) != 1:
            raise ValueError(f"FASTA adı boş olamaz ve boşluk içeremez: {self.__name!r}")
        if not self.__seq:
            raise ValueError(f"{self.__name}: boş dizi FASTA olarak yazılamaz")

        count = sum(1 for key in self.__features if key not in ("fasta_header", "source_file"))
        word = "feature" if count == 1 else "features"
        lines = [f">{self.__name} this sequence contains {count} {word} but not listed"]
        lines += [self.__seq[i:i + 60] for i in range(0, len(self.__seq), 60)]
        with open(file_name, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    def to_crn(self, file_name: str | Path) -> None:
        """Dizinin tüm içeriğini CRN ikili dosyasına yazar.

        Yazılanlar: tip (``_type``), ad, dizi, root alanları (ör. DNA'nın
        strand'i) ve features. pickle kullanılmaz; biçim modülün başındaki
        açıklamada tanımlıdır. Değerler None, bool, int, float, str, bytes,
        CRN_ENUMS içindeki enum'lar (Molecule, Strand, GenomicStrand, Maps,
        GeneticCode) ve
        bunlardan oluşan tuple/list olabilir. Dosya varsa üzerine yazılır.

        Raises:
            TypeError: Desteklenmeyen tipte bir değer varsa. Bu durumda
                dosyaya hiçbir şey yazılmaz.
        """
        _crn_write(file_name, self._crn_record())

    def _crn_record(self) -> bytes:
        """Dizi kaydı: tip, ad, dizi, root alanları ve features (magic/crc olmadan)."""
        root = self._root_fields()
        parts = [
            _crn_pack_str(self._type.name),
            _crn_pack_str(self.__name),
            _crn_pack_str(self.__seq),
            struct.pack("<I", len(root)),
        ]
        for key, value in root.items():
            parts.append(_crn_pack_str(key))
            parts.append(_crn_pack_value(value))
        parts.append(_crn_pack_features(self.__features))
        return b"".join(parts)

    @classmethod
    def from_crn(cls, file_name: str | Path) -> Self:
        """``to_crn`` ile yazılmış CRN dosyasını okur.

        - ``Sequence.from_crn`` dosyadaki tipe uygun sınıfı döndürür: DNA
          dosyası için ``DNA``, SEQ dosyası için ``Sequence``.
        - ``DNA.from_crn`` gibi bir alt sınıftan çağrılırsa dosyadaki tip o
          sınıfla aynı olmalıdır.

        Ad, dizi, root alanları ve tüm features dosyadaki haliyle geri yüklenir.

        Raises:
            FileNotFoundError: Dosya yoksa.
            ValueError: Dosya CRN değilse, sürümü desteklenmiyorsa, bozuksa
                (CRC32 tutmuyorsa) ya da dosyadaki dizi tipi ``cls`` ile
                uyuşmuyorsa.
        """
        reader = _crn_open(file_name)
        result = cls._crn_read(reader)
        reader.finish()
        return result

    @classmethod
    def _crn_read(cls, reader: _CrnReader) -> Self:
        """Okuyucunun bulunduğu yerden bir dizi kaydı okur (``_crn_record``'un tersi).

        ``cls`` Sequence ise dosyadaki tipe uygun sınıf seçilir; alt sınıfsa
        dosyadaki tip o sınıfla aynı olmalıdır.
        """
        type_name = reader.text()
        if type_name == "GENE":
            raise ValueError(f"{reader.source}: bu bir gen dosyası; Gene.from_crn ile açın")
        if type_name not in Molecule.__members__:
            raise ValueError(f"{reader.source}: bilinmeyen dizi tipi: {type_name}")
        molecule = Molecule[type_name]

        if molecule is cls._molecule:
            target = cls
        elif cls is Sequence and molecule in _SEQUENCE_CLASSES:
            target = _SEQUENCE_CLASSES[molecule]
        else:
            raise ValueError(
                f"{reader.source}: dosyadaki dizi tipi {type_name}, "
                f"{cls.__name__} sınıfı {cls._molecule.name} tipi bekliyor"
            )

        name = reader.text()
        seq = reader.text()
        root = {}
        for _ in range(reader.uint("<I")):
            key = reader.text()
            root[key] = reader.value()

        result = target(seq, name=name, **root)
        result.features = reader.features()
        return result

    @property
    def features(self) -> Features:
        return self.__features

    @features.setter
    def features(self, value: Features) -> None:
        if not isinstance(value, Features):
            raise TypeError(f"features bir Features olmalı, {type(value).__name__} verildi")
        self.__features = value

    @property
    def seq(self) -> str:
        return self.__seq

    @seq.setter
    def seq(self, value:str) -> None:
        self.__seq = self._check_seq(value)

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, value:str) -> None:
        self.__name = value

    def __len__(self) -> int:
        return len(self.__seq)

    def __eq__(self, value: object) -> bool:
        """Aynı tip ve aynı dizi ise eşittir; ``str`` ile büyük/küçük harf
        farkı gözetmeden karşılaştırılır. Ad ve features karşılaştırılmaz."""
        if isinstance(value, Sequence):
            return self._type is value._type and self.__seq == value.__seq
        if isinstance(value, str):
            return self.__seq == value.upper()
        return NotImplemented

    # seq değiştirilebildiği için diziler hash'lenemez (set/dict anahtarı olamaz).
    __hash__ = None

    def __contains__(self, item: Sequence | str) -> bool:
        """``"ATG" in s`` veya ``other in s`` biçiminde alt dizi arar."""
        if isinstance(item, Sequence):
            return item.__seq in self.__seq
        if isinstance(item, str):
            return item.upper() in self.__seq
        raise TypeError(f"Sequence içinde aranamayan tip: {type(item).__name__}")

    def __iter__(self) -> Iterator[str]:
        return iter(self.__seq)

    def __getitem__(self, key: int | slice) -> str | Self:
        """``s[i]`` tek harfi (``str``), ``s[a:b]`` aynı sınıftan yeni bir dizi döndürür.

        Yeni dizi adı, root alanlarını ve features'ın bağımsız bir kopyasını taşır.
        Adım 0 ise ``ValueError`` verilir.
        """
        if isinstance(key, slice):
            return self._new(self.__seq[key])
        return self.__seq[key]

    def __repr__(self) -> str:
        seq = self.__seq
        if len(self) >= 80:
            seq = f"{self.__seq[:40]}...{self.__seq[-40:]}"
        return f"{self.__name}: {seq} <{len(self)}> [{self._type.value[1]}]"
