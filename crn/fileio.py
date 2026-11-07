"""Dosya biçimleri: FASTA ve CRN (pickle'sız ikili biçim).

CRN ikili dosya biçimi (tüm sayılar little-endian)::

    magic      4 bayt   b"CRN\\x00"
    version    uint16   CRN_VERSION
    kayıt               dizi kaydı, gen kaydı ya da liste kaydı
    crc32      uint32   magic dahil önceki tüm baytların CRC32'si

    dizi kaydı:
        type       str      Molecule adı ("SEQ", "DNA", "RNA", "PRT")
        name, seq  str
        n_root     uint32 + n_root × (key str, value etiketli değer)
        features

    gen kaydı:
        type       str      "GENE"
        name       str
        features
        has_dna    uint8    1 ise ardından genomik DNA'nın dizi kaydı
        n_transcripts uint32 + n × (id str, dizi kaydı)
        n_proteins    uint32 + n × (id str, dizi kaydı)

    liste kaydı:
        type       str      "LIST"
        n_items    uint32 + n × (dizi kaydı ya da gen kaydı)

    features: uint32 anahtar sayısı + her anahtar için
              key str, uint32 değer sayısı, etiketli değerler
    str:      uint32 bayt uzunluğu + UTF-8 baytları

Etiketli değer: 1 baytlık etiket + içerik::

    N None   T True   F False
    i int    uint8 uzunluk + işaretli little-endian baytlar
    f float  float64
    s str    b bytes (uint32 uzunluk + baytlar)
    e enum   enum sınıfının adı (str) + üye adı (str); yalnızca CRN_ENUMS
    t tuple, l list: uint32 eleman sayısı + etiketli değerler

Sürüm 3 dosyalarında DNA/RNA'nın root alanlarında bulunan ``maps`` (eski
``Maps`` enum'u) okunurken yok sayılır.
"""
from __future__ import annotations

import gzip
import struct
import zlib
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, List, Optional, Type

from crn.constants import GeneticCode
from crn.features import Features
from crn.seq import (Molecule, Strand, GenomicStrand, Sequence, SEQUENCE_CLASSES)

if TYPE_CHECKING:
    from crn.gene import Gene

__all__ = ["read_fasta", "write_fasta", "read_crn", "write_crn", "read_gene_crn",
           "write_gene_crn", "read_crn_list", "write_crn_list", "CRN_VERSION"]

CRN_MAGIC = b"CRN\x00"
CRN_VERSION = 3
CRN_ENUMS = {e.__name__: e for e in (Molecule, Strand, GenomicStrand, GeneticCode)}
# Eski sürümlerde yazılmış, artık kullanılmayan enum'lar; değerleri None okunur.
_LEGACY_ENUMS = frozenset({"Maps"})

FASTA_LINE_WIDTH = 60
# from_fasta_file'ın eklediği, to_fasta başlığındaki sayıma girmeyen anahtarlar.
_FASTA_KEYS = ("fasta_header", "source_file")


# --- FASTA ---------------------------------------------------------------------

def read_fasta(file_name: str | Path, cls: Type[Sequence] = Sequence) -> List[Sequence]:
    """FASTA dosyasındaki dizileri ``cls`` tipinde bir liste olarak okur.

    ``.gz`` dosyalar gzip ile açılır; dizinin adı başlığın ilk kelimesidir.
    Boş satırlar ve ``;`` yorumları atlanır, çok satırlı diziler birleştirilir.
    features'a ``fasta_header`` (başlığın tamamı) ve ``source_file`` yazılır.

    Raises:
        FileNotFoundError: Dosya yoksa.
        ValueError: Başlıksız veri, boş başlık, dizisiz başlık, ``cls`` için
            geçersiz karakter ya da hiç kayıt yoksa (mesajda dosya:satır).
    """
    path = Path(file_name)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as f:
        lines = f.read().splitlines()

    sequences: List[Sequence] = []
    header: Optional[str] = None
    header_line = 0
    chunks: List[str] = []

    def flush() -> None:
        if not chunks:
            raise ValueError(f"{path}:{header_line}: '>{header}' bir seq bloğu olmadan bitti")
        try:
            seq = cls("".join(chunks), name=header.split()[0])
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
            header, header_line, chunks = line[1:].strip(), line_no, []
            if not header:
                raise ValueError(f"{path}:{line_no}: boş FASTA başlığı")
        elif header is None:
            raise ValueError(f"{path}:{line_no}: '>' header'ı olmadan sequence verisi bulundu")
        else:
            chunks.append("".join(line.split()))
    if header is not None:
        flush()
    if not sequences:
        raise ValueError(f"{path}: geçerli bir FASTA kaydı bulunamadı")
    return sequences


def write_fasta(sequences: Iterable[Sequence], file_name: str | Path) -> None:
    """Dizileri FASTA dosyasına yazar (dosya varsa üzerine yazılır).

    Başlık ``>AD this sequence contains N feature(s) but not listed``
    biçimindedir; N, ``fasta_header`` ve ``source_file`` dışındaki anahtar
    sayısıdır. Dizi 60 karakterlik satırlara bölünür.

    Raises:
        ValueError: Bir dizinin adı boşsa/boşluk içeriyorsa ya da dizi boşsa;
            bu durumda hiçbir şey yazılmaz.
    """
    lines = []
    for s in sequences:
        if not s.name or len(s.name.split()) != 1:
            raise ValueError(f"FASTA adı boş olamaz ve boşluk içeremez: {s.name!r}")
        if not s.seq:
            raise ValueError(f"{s.name}: boş dizi FASTA olarak yazılamaz")
        count = sum(1 for key in s.features if key not in _FASTA_KEYS)
        word = "feature" if count == 1 else "features"
        lines.append(f">{s.name} this sequence contains {count} {word} but not listed")
        lines += [s.seq[i:i + FASTA_LINE_WIDTH] for i in range(0, len(s.seq), FASTA_LINE_WIDTH)]
    with open(file_name, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# --- CRN: yazma ---------------------------------------------------------------

def _pack_str(text: str) -> bytes:
    data = text.encode("utf-8")
    return struct.pack("<I", len(data)) + data


def _pack_value(value: Any) -> bytes:
    if value is None:
        return b"N"
    if isinstance(value, bool):
        return b"T" if value else b"F"
    if type(value) in CRN_ENUMS.values():
        return b"e" + _pack_str(type(value).__name__) + _pack_str(value.name)
    if isinstance(value, int):
        data = value.to_bytes((value.bit_length() + 8) // 8, "little", signed=True)
        return b"i" + struct.pack("<B", len(data)) + data
    if isinstance(value, float):
        return b"f" + struct.pack("<d", value)
    if isinstance(value, str):
        return b"s" + _pack_str(value)
    if isinstance(value, bytes):
        return b"b" + struct.pack("<I", len(value)) + value
    if isinstance(value, (tuple, list)):
        tag = b"t" if isinstance(value, tuple) else b"l"
        return tag + struct.pack("<I", len(value)) + b"".join(_pack_value(v) for v in value)
    raise TypeError(f"CRN dosyasına yazılamayan değer tipi: {type(value).__name__}")


def _pack_features(features: Features) -> bytes:
    parts = [struct.pack("<I", len(features))]
    for key, values in features.items():
        parts += [_pack_str(key), struct.pack("<I", len(values))]
        parts += [_pack_value(v) for v in values]
    return b"".join(parts)


def _pack_sequence(seq: Sequence) -> bytes:
    root = seq.root_fields()
    parts = [_pack_str(seq.molecule.name), _pack_str(seq.name), _pack_str(seq.seq),
             struct.pack("<I", len(root))]
    for key, value in root.items():
        parts += [_pack_str(key), _pack_value(value)]
    parts.append(_pack_features(seq.features))
    return b"".join(parts)


def _write(file_name: str | Path, body: bytes) -> None:
    data = CRN_MAGIC + struct.pack("<H", CRN_VERSION) + body
    data += struct.pack("<I", zlib.crc32(data))
    Path(file_name).write_bytes(data)


def write_crn(seq: Sequence, file_name: str | Path) -> None:
    """Diziyi CRN dosyasına yazar.

    Raises:
        TypeError: Desteklenmeyen tipte bir değer varsa (dosyaya hiçbir şey yazılmaz).
    """
    _write(file_name, _pack_sequence(seq))


def _pack_gene(gene: Gene) -> bytes:
    parts = [_pack_str("GENE"), _pack_str(gene.name), _pack_features(gene.features)]
    if gene.dna is None:
        parts.append(struct.pack("<B", 0))
    else:
        parts += [struct.pack("<B", 1), _pack_sequence(gene.dna)]
    for items in (gene.transcripts, gene.proteins):
        parts.append(struct.pack("<I", len(items)))
        for item_id, seq in items.items():
            parts += [_pack_str(item_id), _pack_sequence(seq)]
    return b"".join(parts)


def write_gene_crn(gene: Gene, file_name: str | Path) -> None:
    """Geni (ad, features, genomik DNA, transkriptler, proteinler) CRN dosyasına yazar."""
    _write(file_name, _pack_gene(gene))


def write_crn_list(items: Iterable[Sequence | Gene], file_name: str | Path) -> None:
    """Gene ve Sequence (DNA, RNA, PRT) nesnelerini sırayla tek CRN dosyasına yazar.

    Liste karışık olabilir; dosya varsa üzerine yazılır. Aynı addan birden çok
    kayıt olmasına izin verilir.

    Raises:
        TypeError: Bir eleman Gene ya da Sequence değilse ya da desteklenmeyen
            tipte bir değer içeriyorsa (dosyaya hiçbir şey yazılmaz).
    """
    from crn.gene import Gene

    records = []
    for item in items:
        if isinstance(item, Sequence):
            records.append(_pack_sequence(item))
        elif isinstance(item, Gene):
            records.append(_pack_gene(item))
        else:
            raise TypeError(f"CRN listesi elemanı Gene ya da Sequence olmalı, "
                            f"{type(item).__name__} verildi")
    _write(file_name, _pack_str("LIST") + struct.pack("<I", len(records)) + b"".join(records))


# --- CRN: okuma ----------------------------------------------------------------

class _Reader:
    """CRN baytlarını sırayla okur; eksik veride ValueError verir."""

    def __init__(self, file_name: str | Path):
        path = Path(file_name)
        data = path.read_bytes()
        self.source = str(path)
        if len(data) < len(CRN_MAGIC) + 6 or not data.startswith(CRN_MAGIC):
            raise ValueError(f"{path}: CRN dosyası değil")
        body, crc = data[:-4], struct.unpack("<I", data[-4:])[0]
        if zlib.crc32(body) != crc:
            raise ValueError(f"{path}: CRN dosyası bozuk (CRC32 uyuşmuyor)")
        self.data, self.pos = body, len(CRN_MAGIC)
        version = self.uint("<H")
        if version != CRN_VERSION:
            raise ValueError(f"{path}: desteklenmeyen CRN sürümü: {version}")

    def take(self, n: int) -> bytes:
        if self.pos + n > len(self.data):
            raise ValueError(f"{self.source}: CRN dosyası beklenenden kısa")
        chunk = self.data[self.pos:self.pos + n]
        self.pos += n
        return chunk

    def uint(self, fmt: str) -> Any:
        return struct.unpack(fmt, self.take(struct.calcsize(fmt)))[0]

    def text(self) -> str:
        return self.take(self.uint("<I")).decode("utf-8")

    def finish(self) -> None:
        if self.pos != len(self.data):
            raise ValueError(f"{self.source}: CRN dosyasının sonunda fazladan veri var")

    def value(self) -> Any:
        tag = self.take(1)
        if tag == b"N":
            return None
        if tag in (b"T", b"F"):
            return tag == b"T"
        if tag == b"i":
            return int.from_bytes(self.take(self.uint("<B")), "little", signed=True)
        if tag == b"f":
            return self.uint("<d")
        if tag == b"s":
            return self.text()
        if tag == b"b":
            return self.take(self.uint("<I"))
        if tag == b"e":
            enum_name, name = self.text(), self.text()
            if enum_name in _LEGACY_ENUMS:
                return None
            enum_cls = CRN_ENUMS.get(enum_name)
            if enum_cls is None:
                raise ValueError(f"{self.source}: bilinmeyen enum: {enum_name}")
            if name not in enum_cls.__members__:
                raise ValueError(f"{self.source}: bilinmeyen {enum_name} üyesi: {name}")
            return enum_cls[name]
        if tag in (b"t", b"l"):
            items = [self.value() for _ in range(self.uint("<I"))]
            return tuple(items) if tag == b"t" else items
        raise ValueError(f"{self.source}: bilinmeyen değer etiketi: {tag!r}")

    def features(self) -> Features:
        result = Features()
        for _ in range(self.uint("<I")):
            key = self.text()
            for _ in range(self.uint("<I")):
                result.add(key, self.value())
        return result

    def sequence(self, cls: Type[Sequence], type_name: Optional[str] = None) -> Sequence:
        """Bir dizi kaydı okur; ``cls`` Sequence ise dosyadaki tipin sınıfı seçilir.

        ``type_name`` verilmişse kaydın tipi zaten okunmuştur.
        """
        type_name = self.text() if type_name is None else type_name
        if type_name == "GENE":
            raise ValueError(f"{self.source}: bu bir gen dosyası; Gene.from_crn ile açın")
        if type_name == "LIST":
            raise ValueError(f"{self.source}: bu bir liste dosyası; read_crn_list ile açın")
        if type_name not in Molecule.__members__:
            raise ValueError(f"{self.source}: bilinmeyen dizi tipi: {type_name}")
        molecule = Molecule[type_name]
        if molecule is cls.molecule:
            target = cls
        elif cls is Sequence and molecule in SEQUENCE_CLASSES:
            target = SEQUENCE_CLASSES[molecule]
        else:
            raise ValueError(f"{self.source}: dosyadaki dizi tipi {type_name}, "
                             f"{cls.__name__} sınıfı {cls.molecule.name} tipi bekliyor")
        name, seq = self.text(), self.text()
        root = {}
        for _ in range(self.uint("<I")):
            key, value = self.text(), self.value()
            if key in target.root_types:
                root[key] = value
        result = target(seq, name=name, **root)
        result.features = self.features()
        return result

    def gene(self) -> Gene:
        """Tipi ("GENE") okunmuş bir gen kaydını okur."""
        from crn.gene import Gene
        from crn.seq import DNA, RNA, PRT

        gene = Gene(self.text())
        gene.features = self.features()
        if self.uint("<B"):
            gene.dna = self.sequence(DNA)
        for _ in range(self.uint("<I")):
            tid = self.text()
            gene.add_transcript(self.sequence(RNA), tid)
        for _ in range(self.uint("<I")):
            pid = self.text()
            gene.add_protein(self.sequence(PRT), pid)
        return gene


def read_crn(file_name: str | Path, cls: Type[Sequence] = Sequence) -> Sequence:
    """``write_crn`` ile yazılmış dizi dosyasını okur.

    Raises:
        FileNotFoundError: Dosya yoksa.
        ValueError: Dosya CRN değilse, bozuksa, sürümü desteklenmiyorsa ya da
            dosyadaki tip ``cls`` ile uyuşmuyorsa.
    """
    reader = _Reader(file_name)
    result = reader.sequence(cls)
    reader.finish()
    return result


def read_gene_crn(file_name: str | Path) -> Gene:
    """``write_gene_crn`` ile yazılmış gen dosyasını okur.

    Raises:
        ValueError: Dosya bir gen dosyası değilse ya da bozuksa.
    """
    reader = _Reader(file_name)
    kind = reader.text()
    if kind != "GENE":
        raise ValueError(f"{reader.source}: gen dosyası değil (dosya tipi: {kind})")
    gene = reader.gene()
    reader.finish()
    return gene


def read_crn_list(file_name: str | Path) -> List[Sequence | Gene]:
    """``write_crn_list`` ile yazılmış dosyayı yazıldığı sırayla liste olarak okur.

    Her eleman dosyadaki tipine uygun sınıftadır: Gene, DNA, RNA, PRT ya da Sequence.

    Raises:
        FileNotFoundError: Dosya yoksa.
        ValueError: Dosya bir liste dosyası değilse ya da bozuksa.
    """
    reader = _Reader(file_name)
    kind = reader.text()
    if kind != "LIST":
        raise ValueError(f"{reader.source}: liste dosyası değil (dosya tipi: {kind})")
    items: List[Sequence | Gene] = []
    for _ in range(reader.uint("<I")):
        kind = reader.text()
        items.append(reader.gene() if kind == "GENE" else reader.sequence(Sequence, kind))
    reader.finish()
    return items
