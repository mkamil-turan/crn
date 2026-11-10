"""crn komut satırı arayüzü.

Kurulumdan sonra ``crn`` komutuyla, kurulum olmadan ``python -m crn.cli`` ile çalışır::

    crn seq       dizi işlemleri (info, revcomp, transcribe, reverse-transcribe, translate, orfs)
    crn get       tek kayıt indirme (dna, rna, protein, region)
    crn gene      gen özeti, FASTA ya da CRN çıktısı
    crn cub       kodon kullanımı (table, codons, plot)
    crn dco       iki gen grubunun diferansiyel kodon karşılaştırması
    crn protein   protein özellikleri
    crn hpo       semptomdan gene sorgular (search, genes, symbols, phenotypes, diseases, summary, pack)
    crn hgnc      gen tablosu ve ortologlar (columns, symbols, query, orthologs)
    crn clinvar   klinik varyantlar (variants, genes, symbols, conditions, summary, pack)
    crn taxonomy  organizma adından NCBI Taxonomy ID'sine (ve tersi)
    crn update    yerel veri deposunu hazırlama

Çıktı kuralları: veri stdout'a (ya da ``-o`` dosyasına), ilerleme, uyarı ve hatalar
stderr'e yazılır, böylece komutlar boruyla bağlanabilir. Tablolar sekmeyle ayrılmıştır,
diziler FASTA'dır ve FASTA çıktısı her zaman 5'->3' okunuştadır. Çıkış kodları: 0 başarı,
1 işlem hatası, 2 kullanım hatası, 130 kullanıcı kesmesi.

Gen girdisi (``GENLER``) gen sembolü, virgülle ayrılmış semboller, ``.crn`` dosyası, satır
başına bir sembol içeren metin dosyası ya da stdin (``-``) olabilir. Dizi girdisi doğrudan
dizi, FASTA, ``.crn`` dosyası ya da stdin'dir.

``fetch`` ve ``cub`` paketleri yalnızca kullanan komutta yüklenir; bu sayede ``crn --help``
scipy ve matplotlib beklemez.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import warnings
from importlib import metadata
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence as Seq

from crn.constants import GeneticCode
from crn.fileio import read_crn_list, read_fasta, write_crn_list
from crn.gene import Gene
from crn.protein import ProteinProperties
from crn.seq import DNA, PRT, RNA, NucleicAcid, Sequence

__all__ = ["main", "build_parser"]

GENETIC_CODES = {"standard": GeneticCode.STANDARD, "mito": GeneticCode.VERTEBRATE_MITOCHONDRIAL}
SEQUENCE_TYPES = {"dna": DNA, "rna": RNA, "protein": PRT}
FASTA_WIDTH = 60
# Otomatik tip tanıma: IUPAC harflerinin çoğu aynı zamanda amino asit olduğundan
# (M, K, W, V ...) yalnızca bu harflerden oluşan dizi nükleik asit sayılır.
_NUCLEIC_LETTERS = frozenset("ACGTUN-")
_FILE_SUFFIXES = (".fa", ".fasta", ".fna", ".faa", ".crn", ".gz", ".txt", ".tsv")
PLOT_KINDS = ("enc", "enc-ratio", "neutrality", "pr2", "rscu", "rscu-heatmap", "cai", "coa",
              "correlations")
DCO_REPORTS = ("robust", "fisher", "significant", "permutation", "levels", "interpret")


# --- çıktı ---------------------------------------------------------------------

def _log(message: str) -> None:
    print(f"crn: {message}", file=sys.stderr)


def _emit(text: str, output: Optional[str]) -> None:
    """Metni stdout'a ya da ``output`` dosyasına yazar."""
    if output is None or output == "-":
        sys.stdout.write(text)
    else:
        Path(output).write_text(text, encoding="utf-8")
        _log(f"yazıldı: {output}")


def _cell(value: Any) -> str:
    """Tablo hücresi: None boş, ondalıklar 6 anlamlı basamak."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _write_table(rows: Seq[Dict[str, Any]], output: Optional[str]) -> None:
    """Sözlük satırlarını TSV yazar; sütunlar ilk görülme sırasındadır."""
    if not rows:
        _log("sonuç yok")
        return
    columns = list(dict.fromkeys(key for row in rows for key in row))
    lines = ["\t".join(columns)]
    lines += ["\t".join(_cell(row.get(column)) for column in columns) for row in rows]
    _emit("\n".join(lines) + "\n", output)


def _write_lines(lines: Iterable[str], output: Optional[str]) -> None:
    text = "".join(f"{line}\n" for line in lines)
    if text:
        _emit(text, output)
    else:
        _log("sonuç yok")


def _forward(seq: Sequence) -> str:
    """Dizinin 5'->3' okunuşu.

    FASTA yön bilgisi taşımaz ve okunurken 5'->3' kabul edilir; 3'->5' yazılmış bir
    tümleyici olduğu gibi yazılırsa yanlış zincir olarak okunur.
    """
    return seq.seq_5to3() if isinstance(seq, NucleicAcid) else seq.seq


def _fasta_text(sequences: Iterable[Sequence]) -> str:
    lines: List[str] = []
    for seq in sequences:
        text = _forward(seq)
        lines.append(f">{seq.name or 'seq'}")
        lines += [text[i:i + FASTA_WIDTH] for i in range(0, len(text), FASTA_WIDTH)]
    return "".join(f"{line}\n" for line in lines)


def _write_sequences(sequences: List[Sequence], output: Optional[str], raw: bool = False) -> None:
    """Dizileri FASTA, ``raw`` ile yalın satırlar, ``.crn`` uzantısında CRN olarak yazar."""
    if not sequences:
        _log("sonuç yok")
    elif output and output.endswith(".crn"):
        if len(sequences) == 1:
            sequences[0].to_crn(output)
        else:
            write_crn_list(sequences, output)
        _log(f"yazıldı: {output}")
    elif raw:
        _write_lines((_forward(seq) for seq in sequences), output)
    else:
        _emit(_fasta_text(sequences), output)


# --- dizi girdisi --------------------------------------------------------------

def _load_crn(path: str) -> List[Sequence | Gene]:
    """Dizi, gen ya da liste CRN dosyasını türüne bakmadan liste olarak okur."""
    try:
        return read_crn_list(path)
    except ValueError as error:
        message = str(error)
    if "liste dosyası değil" not in message:
        raise ValueError(message)
    return [Gene.from_crn(path) if "GENE" in message else Sequence.from_crn(path)]


def _read_stdin(what: str) -> str:
    if sys.stdin.isatty():
        raise ValueError(f"{what} verilmedi; argüman yazın ya da stdin'den yönlendirin")
    return sys.stdin.read()


def _parse_text(text: str, name: str) -> List[Sequence]:
    """stdin metni: ``>`` ile başlıyorsa FASTA, değilse yalın dizi."""
    if text.lstrip().startswith(">"):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "stdin.fa"
            path.write_text(text, encoding="utf-8")
            return read_fasta(path)
    return [Sequence("".join(text.split()).upper(), name=name)]


def _read_sequences(spec: str, index: int) -> List[Sequence]:
    """Tek bir dizi girdisini okur: stdin, dosya ya da doğrudan dizi."""
    name = f"seq{index}"
    if spec == "-":
        return _parse_text(_read_stdin("dizi"), name)
    if Path(spec).is_file():
        if spec.endswith(".crn"):
            items = _load_crn(spec)
            if any(isinstance(item, Gene) for item in items):
                raise ValueError(f"{spec}: gen dosyası; bu komut dizi bekler")
            return items
        return read_fasta(spec)
    if "/" in spec or spec.lower().endswith(_FILE_SUFFIXES):
        raise FileNotFoundError(f"Dosya yok: {spec}")
    return [Sequence(spec.upper(), name=name)]


def _coerce(seq: Sequence, kind: str) -> Sequence:
    """Diziyi istenen sınıfa çevirir.

    ``kind``: ``dna``, ``rna``, ``protein``, ``nucleic`` (U varsa RNA, yoksa DNA) ya da
    ``auto``. CRN'den gelen tipli diziler ``auto`` ve ``nucleic``te olduğu gibi kalır.
    """
    if type(seq) is not Sequence and kind in ("auto", "nucleic"):
        return seq
    if kind == "auto":
        kind = "nucleic" if set(seq.seq.upper()) <= _NUCLEIC_LETTERS else "protein"
    if kind == "nucleic":
        kind = "rna" if "U" in seq.seq.upper() else "dna"
    cls = SEQUENCE_TYPES[kind]
    return seq if isinstance(seq, cls) else cls.from_sequence(seq)


def _sequences(args: argparse.Namespace, default_kind: str) -> List[Sequence]:
    """``args.inputs`` girdilerini okur; girdi yoksa stdin kullanılır."""
    kind = getattr(args, "type", None) or default_kind
    result: List[Sequence] = []
    for spec in args.inputs or ["-"]:
        result += [_coerce(seq, kind) for seq in _read_sequences(spec, len(result) + 1)]
    return result


# --- veri kaynakları -----------------------------------------------------------

def _verbose(args: argparse.Namespace) -> bool:
    """İlerleme çubukları yalnızca terminalde gösterilir; günlük dosyasını kirletmez."""
    return not args.quiet and sys.stderr.isatty()


def _ncbi(args: argparse.Namespace):
    from crn.fetch import Ncbi, NcbiSession
    session = NcbiSession(api_key=os.environ.get("NCBI_API_KEY"))
    return Ncbi(organism=args.organism, session=session, sys_dir=args.sys_dir, verbose=_verbose(args))


def _gencode(args: argparse.Namespace, species: str = "human", release: Optional[str] = None,
             prepare: bool = False):
    """``Gencode`` nesnesi; yerel veri yoksa sessizce 1,3 GB indirmek yerine durur."""
    from crn.fetch import Gencode
    source = Gencode(species=species, release=release, sys_dir=args.sys_dir, verbose=_verbose(args))
    if not prepare and not source.is_ready():
        raise ValueError("GENCODE verisi hazır değil; 1,3 GB indirir, açılmış hâli 5 GB'tır. "
                         "Hazırlamak için: crn update gencode")
    return source


def _source(args: argparse.Namespace):
    if args.source == "gencode":
        return _gencode(args, args.organism, args.release)
    return _ncbi(args)


def _fetch_gene(source: Any, symbol: str, args: argparse.Namespace) -> Gene:
    if args.canonical:
        if args.source != "gencode":
            raise ValueError("--canonical yalnızca --source gencode ile kullanılır")
        return source.canonical(symbol)
    return source.gene(symbol)


def _symbols(spec: str) -> List[str] | None:
    """Girdi sembol listesiyse sembolleri, ``.crn`` dosyasıysa None döndürür."""
    if spec == "-":
        return _read_stdin("gen listesi").split()
    if spec.endswith(".crn") and Path(spec).is_file():
        return None
    if Path(spec).is_file():
        return Path(spec).read_text(encoding="utf-8").split()
    if spec.endswith(".crn"):
        raise FileNotFoundError(f"Dosya yok: {spec}")
    return [symbol for symbol in spec.split(",") if symbol]


def _load_genes(specs: Seq[str], args: argparse.Namespace) -> List[Gene]:
    """Gen girdilerini ``CUB``'a hazır, tek izoformlu ``Gene`` listesine çevirir.

    Listedeki tek bir sembolün bulunamaması ya da proteinle eşleşmemesi bütün işi
    durdurmaz: gen atlanır ve bildirilir. Hiç gen kalmazsa hata verilir.
    """
    code = GENETIC_CODES[args.genetic_code]
    genes: List[Gene] = []
    skipped: Dict[str, str] = {}
    source = None
    for spec in specs:
        symbols = _symbols(spec)
        if symbols is None:
            items = _load_crn(spec)
            if not all(isinstance(item, Gene) for item in items):
                raise ValueError(f"{spec}: gen içermeyen kayıtlar var")
            for gene in items:
                single = len(gene.transcripts) == 1 and len(gene.proteins) == 1
                genes += [gene] if single else gene.isoforms(genetic_code=code)
            continue
        source = source or _source(args)
        for symbol in symbols:
            try:
                isoforms = _fetch_gene(source, symbol, args).isoforms(genetic_code=code)
            except (KeyError, ValueError) as error:
                skipped[symbol] = str(error.args[0]) if error.args else type(error).__name__
                continue
            if isoforms:
                genes += isoforms
            else:
                skipped[symbol] = "proteinle eşleşen izoform yok"
    for symbol, reason in skipped.items():
        warnings.warn(f"{symbol} atlandı: {reason}", UserWarning)
    if not genes:
        raise ValueError("Hesaplanacak gen kalmadı")
    return genes


def _cub(specs: Seq[str], args: argparse.Namespace):
    from crn.cub import CUB
    return CUB(_load_genes(specs, args), genetic_code=GENETIC_CODES[args.genetic_code])


def _reference_sets(args: argparse.Namespace) -> Dict[str, Any]:
    """``summary_table`` referans argümanları; yalnızca verilenler."""
    result: Dict[str, Any] = {}
    for key in ("reference", "background", "host", "genome"):
        spec = getattr(args, key, None)
        if spec:
            result[key] = _cub([spec], args)
    if getattr(args, "trna", None):
        from crn.cub import read_trna_counts
        result["trna_counts"] = read_trna_counts(args.trna)
    return result


# --- komutlar: dizi ------------------------------------------------------------

def _cmd_seq_info(args: argparse.Namespace) -> None:
    rows = []
    for seq in _sequences(args, "auto"):
        text = seq.seq
        gc = (text.count("G") + text.count("C")) / len(text) * 100 if text else 0.0
        rows.append({"name": seq.name, "type": seq.molecule.name, "length": len(text),
                     "gc_percent": gc if isinstance(seq, NucleicAcid) else None})
    _write_table(rows, args.output)


def _cmd_seq_revcomp(args: argparse.Namespace) -> None:
    _write_sequences([s.reverse_complement() for s in _sequences(args, "nucleic")],
                     args.output, args.raw)


def _cmd_seq_transcribe(args: argparse.Namespace) -> None:
    template = args.strand == "template"
    _write_sequences([s.transcription(template=template) for s in _sequences(args, "dna")],
                     args.output, args.raw)


def _cmd_seq_reverse_transcribe(args: argparse.Namespace) -> None:
    _write_sequences([s.reverse_transcription() for s in _sequences(args, "rna")],
                     args.output, args.raw)


def _cmd_seq_translate(args: argparse.Namespace) -> None:
    code = GENETIC_CODES[args.genetic_code]
    result = []
    for seq in _sequences(args, "nucleic"):
        if isinstance(seq, PRT):
            raise ValueError(f"{seq.name}: protein çevrilemez")
        rna = seq if isinstance(seq, RNA) else seq.transcription(template=False)
        result.append(rna.translation(code))
    _write_sequences(result, args.output, args.raw)


def _cmd_seq_orfs(args: argparse.Namespace) -> None:
    code = GENETIC_CODES[args.genetic_code]
    orfs: List[DNA] = []
    for seq in _sequences(args, "dna"):
        orfs += seq.find_orfs(args.min_length, code, args.alternative_starts)
    if args.tsv:
        _write_table([{"name": orf.name, "start": orf.features.first("orf_start"),
                       "end": orf.features.first("orf_end"),
                       "strand": orf.features.first("orf_strand"),
                       "frame": orf.features.first("orf_frame"), "length": len(orf.seq)}
                      for orf in orfs], args.output)
    else:
        _write_sequences(orfs, args.output, args.raw)


def _cmd_protein(args: argparse.Namespace) -> None:
    rows = []
    for seq in _sequences(args, "protein"):
        rows.append({"name": seq.name, **ProteinProperties(seq).summary()})
    _write_table(rows, args.output)


# --- komutlar: veri erişimi ----------------------------------------------------

def _cmd_get(args: argparse.Namespace) -> None:
    source = _source(args)
    fetch = {"dna": source.dna, "rna": source.rna, "protein": source.protein}[args.kind]
    _write_sequences([fetch(identifier) for identifier in args.ids], args.output, args.raw)


def _cmd_taxonomy(args: argparse.Namespace) -> None:
    """Taxid'lerin kaydını ya da addan bulunan taxid'leri yazar."""
    from crn.fetch._store import parse_taxid
    source = _ncbi(args)
    rows = []
    for query in args.queries:
        found = [source.taxonomy(query)] if parse_taxid(query) else source.find_taxid(query, args.limit)
        if not found:
            _log(f"uyarı: '{query}' NCBI Taxonomy'de bulunamadı")
        rows += [{"query": query, **record} for record in found]
    if not rows:
        raise ValueError("Hiçbir organizma bulunamadı")
    _write_table(rows, args.output)


def _cmd_region(args: argparse.Namespace) -> None:
    region = _source(args).region(args.accession, args.start, args.end, args.strand)
    _write_sequences([region], args.output, args.raw)


def _gene_sequences(gene: Gene, kind: str) -> List[Sequence]:
    """Genin istenen dizileri; ``cds`` transkriptlerin ``cds`` koordinatından kesilir."""
    if kind == "dna":
        return [gene.dna] if gene.dna is not None else []
    if kind == "rna":
        return list(gene.transcripts.values())
    if kind == "protein":
        return list(gene.proteins.values())
    result = []
    for rna in gene.transcripts.values():
        location = rna.features.first("cds")
        if location:
            result.append(rna[location[0] - 1:location[1]])
    return result


def _cmd_gene(args: argparse.Namespace) -> None:
    source = _source(args)
    genes = [_fetch_gene(source, symbol, args) for symbol in args.symbols]
    if args.isoforms:
        code = GENETIC_CODES[args.genetic_code]
        genes = [isoform for gene in genes for isoform in gene.isoforms(genetic_code=code)]
    if args.fasta:
        _write_sequences([s for gene in genes for s in _gene_sequences(gene, args.fasta)],
                         args.output, args.raw)
    elif args.output and args.output.endswith(".crn"):
        if not genes:
            raise ValueError("Yazılacak gen yok")
        if len(genes) == 1:
            genes[0].to_crn(args.output)
        else:
            write_crn_list(genes, args.output)
        _log(f"yazıldı: {args.output}")
    else:
        rows: List[Dict[str, Any]] = []
        for gene in genes:
            if gene.dna is not None:
                rows.append({"gene": gene.name, "kind": "dna", "id": gene.dna.name,
                             "length": len(gene.dna.seq)})
            rows += [{"gene": gene.name, "kind": "rna", "id": key, "length": len(rna.seq)}
                     for key, rna in gene.transcripts.items()]
            rows += [{"gene": gene.name, "kind": "protein", "id": key, "length": len(prt.seq)}
                     for key, prt in gene.proteins.items()]
        _write_table(rows, args.output)


def _cmd_update(args: argparse.Namespace) -> None:
    from crn.fetch import Clinvar, Hgnc, Hpo
    if args.target == "ncbi":
        if not args.genes:
            raise ValueError("crn update ncbi için en az bir gen gerekir")
        source = _ncbi(args)
        for gene in args.genes:
            _log(f"{gene}: NCBI Gene ID {source.update(gene, force=args.force)}")
    elif args.target == "gencode":
        _log(f"GENCODE {_gencode(args, args.organism, args.release, prepare=True).update(force=args.force)} hazır")
    elif args.target == "hpo":
        source = Hpo(release=args.release, sys_dir=args.sys_dir, verbose=_verbose(args))
        _log(f"HPO {source.update(force=args.force)} hazır")
    elif args.target == "clinvar":
        source = Clinvar(release=args.release, assembly=args.assembly, sys_dir=args.sys_dir,
                         verbose=_verbose(args))
        _log(f"ClinVar {source.update(force=args.force)} ({args.assembly}) hazır")
    else:
        species = args.species.split(",") if args.species else None
        source = Hgnc(sys_dir=args.sys_dir, verbose=_verbose(args))
        path = source.update(force=args.force, species=species)
        _log(f"HGNC hazır: {path} ({len(source.prepared_species())} ortolog türü)")


# --- komutlar: kodon kullanımı -------------------------------------------------

def _cmd_cub_table(args: argparse.Namespace) -> None:
    from crn.cub import summary_table
    rows = summary_table(_cub(args.genes, args), kingdom=args.kingdom,
                         protein=not args.no_protein, **_reference_sets(args))
    _write_table(rows, args.output)


def _cmd_cub_codons(args: argparse.Namespace) -> None:
    from crn.cub import Synonymous
    cub = _cub(args.genes, args)
    counts = cub.codon_counts(args.gene)
    per_thousand = cub.codon_per_thousand(args.gene)
    rscu = Synonymous(cub).rscu(args.gene)
    _write_table([{"codon": codon, "amino_acid": cub.genetic_code.table[codon], "count": count,
                   "per_thousand": per_thousand.get(codon), "rscu": rscu.get(codon)}
                  for codon, count in counts.items()], args.output)


def _cmd_cub_plot(args: argparse.Namespace) -> None:
    import matplotlib
    matplotlib.use("Agg")
    from crn.cub import (Adaptation, CorrespondenceAnalysis, EffectiveCodons, Synonymous,
                         plot_correlations)
    cub = _cub(args.genes, args)
    references = _reference_sets(args)
    file_name = args.output
    if args.kind in ("enc", "enc-ratio", "neutrality", "pr2"):
        effective = EffectiveCodons(cub)
        if args.kind == "enc":
            effective.plot_enc(labels=args.labels, file_name=file_name)
        elif args.kind == "enc-ratio":
            effective.plot_enc_ratio(file_name=file_name)
        elif args.kind == "neutrality":
            effective.plot_neutrality(labels=args.labels, file_name=file_name)
        else:
            effective.plot_pr2(labels=args.labels, file_name=file_name)
    elif args.kind == "rscu":
        Synonymous(cub).plot_rscu(args.gene, file_name=file_name)
    elif args.kind == "rscu-heatmap":
        Synonymous(cub).plot_rscu_heatmap(file_name=file_name)
    elif args.kind == "cai":
        Adaptation(cub, reference=references.get("reference")).plot_cai(file_name=file_name)
    elif args.kind == "coa":
        CorrespondenceAnalysis(cub).plot_coa(labels=args.labels, file_name=file_name)
    else:
        plot_correlations(cub, file_name=file_name, **references)
    _log(f"yazıldı: {file_name}")


def _cmd_dco(args: argparse.Namespace) -> None:
    from crn.cub import CUB, DifferentialCodonUsage
    code = GENETIC_CODES[args.genetic_code]
    if args.file:
        dcu = DifferentialCodonUsage.from_crn(args.file, args.group_a, args.group_b,
                                              key=args.key, genetic_code=code)
    else:
        dcu = DifferentialCodonUsage(CUB(_load_genes([args.group_a], args), genetic_code=code),
                                     CUB(_load_genes([args.group_b], args), genetic_code=code),
                                     name_a=args.name_a or "A", name_b=args.name_b or "B")
    if args.report == "interpret":
        _emit(dcu.interpret(q=args.fdr, permutations=args.permutations, seed=args.seed) + "\n",
              args.output)
    else:
        rows = {
            "robust": lambda: dcu.robust(args.level, args.fdr, args.min_log2_or,
                                         args.permutations, args.seed),
            "fisher": lambda: dcu.table(args.level, args.fdr),
            "significant": lambda: dcu.significant(args.level, args.fdr, args.min_log2_or),
            "permutation": lambda: dcu.permutation(args.level, args.fdr, args.permutations,
                                                   args.seed),
            "levels": lambda: dcu.compare_levels(args.fdr, args.permutations, args.seed),
        }[args.report]()
        _write_table(rows, args.output)
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        dcu.plot_volcano(args.level, args.fdr, file_name=args.plot)
        _log(f"yazıldı: {args.plot}")


# --- komutlar: HPO ve HGNC -----------------------------------------------------

def _hpo(args: argparse.Namespace):
    from crn.fetch import Hpo
    return Hpo(release=args.release, sys_dir=args.sys_dir, verbose=_verbose(args))


def _hpo_filters(args: argparse.Namespace) -> Dict[str, Any]:
    """Komut satırında verilen ``Hpo.annotations`` süzgeçleri."""
    names = ("exclude", "disease", "min_frequency", "association", "evidence", "onset",
             "inheritance")
    filters = {name: getattr(args, name) for name in names if getattr(args, name) is not None}
    if args.no_descendants:
        filters["descendants"] = False
    return filters


def _cmd_hpo_search(args: argparse.Namespace) -> None:
    _write_table(_hpo(args).search(args.text, limit=args.limit), args.output)


def _cmd_hpo_genes(args: argparse.Namespace) -> None:
    _write_table(_hpo(args).genes(args.terms, **_hpo_filters(args)), args.output)


def _cmd_hpo_symbols(args: argparse.Namespace) -> None:
    _write_lines(_hpo(args).symbols(args.terms, **_hpo_filters(args)), args.output)


def _cmd_hpo_phenotypes(args: argparse.Namespace) -> None:
    _write_table(_hpo(args).phenotypes(args.gene, **_hpo_filters(args)), args.output)


def _cmd_hpo_diseases(args: argparse.Namespace) -> None:
    _write_table(_hpo(args).diseases(args.name, **_hpo_filters(args)), args.output)


def _cmd_hpo_summary(args: argparse.Namespace) -> None:
    summary = _hpo(args).summary()
    _write_table([{"field": key, "value": value} for key, value in summary.items()], args.output)


def _cmd_hpo_pack(args: argparse.Namespace) -> None:
    report = _hpo(args).pack(_groups(args.groups, "TERİM"), args.output, source=_gencode(args),
                             key=args.key, overlap=args.overlap, **_hpo_filters(args))
    _pack_report(report, args.output)


def _groups(items: Seq[str], what: str) -> Dict[str, str]:
    """``AD=DEĞER`` argümanlarını grup sözlüğüne çevirir."""
    groups = {}
    for item in items:
        name, separator, value = item.partition("=")
        if not separator or not name or not value:
            raise ValueError(f"Grup AD={what} biçiminde olmalı: {item!r}")
        groups[name] = value
    return groups


def _pack_report(report: Dict[str, Any], output: str) -> None:
    for name, count in report["groups"].items():
        _log(f"{name}: {count} gen")
    if report["overlap"]:
        _log(f"gruplar arası çakışan gen: {len(report['overlap'])}")
    if report["skipped"]:
        _log(f"atlanan gen: {len(report['skipped'])}")
    _log(f"yazıldı: {output}")


def _clinvar(args: argparse.Namespace):
    from crn.fetch import Clinvar
    return Clinvar(release=args.release, assembly=args.assembly, sys_dir=args.sys_dir,
                   verbose=_verbose(args))


def _clinvar_filters(args: argparse.Namespace) -> Dict[str, Any]:
    """Komut satırında verilen ``Clinvar.variants`` süzgeçleri."""
    names = ("condition", "significance", "min_stars", "type", "origin", "effect", "chromosome",
             "start", "end")
    filters = {name: getattr(args, name) for name in names if getattr(args, name) is not None}
    filters["max_genes"] = args.max_genes or None
    return filters


def _cmd_clinvar_variants(args: argparse.Namespace) -> None:
    rows = _clinvar(args).variants(args.genes or None, limit=args.limit, **_clinvar_filters(args))
    _write_table(rows, args.output)


def _cmd_clinvar_genes(args: argparse.Namespace) -> None:
    _write_table(_clinvar(args).genes(**_clinvar_filters(args)), args.output)


def _cmd_clinvar_symbols(args: argparse.Namespace) -> None:
    _write_lines(_clinvar(args).symbols(**_clinvar_filters(args)), args.output)


def _cmd_clinvar_conditions(args: argparse.Namespace) -> None:
    filters = _clinvar_filters(args)
    if args.gene:
        filters["gene"] = args.gene
    rows = _clinvar(args).conditions(args.name, informative=not args.all, limit=args.limit, **filters)
    _write_table(rows, args.output)


def _cmd_clinvar_summary(args: argparse.Namespace) -> None:
    summary = _clinvar(args).summary()
    _write_table([{"field": key, "value": value} for key, value in summary.items()], args.output)


def _cmd_clinvar_pack(args: argparse.Namespace) -> None:
    report = _clinvar(args).pack(_groups(args.groups, "KOŞUL"), args.output, source=_gencode(args),
                                 key=args.key, overlap=args.overlap, **_clinvar_filters(args))
    _pack_report(report, args.output)


def _hgnc_filters(items: Seq[str]) -> Dict[str, Any]:
    """``sütun=değer`` çiftleri; virgüllü değer öğelerinden biriyle eşleşir."""
    filters: Dict[str, Any] = {}
    for item in items:
        column, separator, value = item.partition("=")
        if not separator or not column:
            raise ValueError(f"Süzgeç sütun=değer biçiminde olmalı: {item!r}")
        filters[column] = value.split(",") if "," in value else value
    return filters


def _hgnc(args: argparse.Namespace):
    from crn.fetch import Hgnc
    return Hgnc(sys_dir=args.sys_dir, verbose=_verbose(args))


def _fields(args: argparse.Namespace) -> Optional[List[str]]:
    return args.fields.split(",") if args.fields else None


def _cmd_hgnc_columns(args: argparse.Namespace) -> None:
    _write_lines(_hgnc(args).columns, args.output)


def _cmd_hgnc_symbols(args: argparse.Namespace) -> None:
    _write_lines(_hgnc(args).symbols(**_hgnc_filters(args.filters)), args.output)


def _cmd_hgnc_query(args: argparse.Namespace) -> None:
    _write_table(_hgnc(args).query(_fields(args), **_hgnc_filters(args.filters)), args.output)


def _cmd_hgnc_orthologs(args: argparse.Namespace) -> None:
    _write_table(_hgnc(args).orthologs(_fields(args), **_hgnc_filters(args.filters)), args.output)


# --- argümanlar ----------------------------------------------------------------

def _common_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("genel seçenekler")
    group.add_argument("--sys-dir", metavar="KLASÖR", help="veri deposu (varsayılan: CRN_SYS_DIR ya da sys/)")
    group.add_argument("--quiet", action="store_true", help="ilerleme çubuklarını kapatır")
    group.add_argument("--debug", action="store_true", help="hatada tam izleme yığınını gösterir")
    return parser


def _output_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-o", "--output", metavar="DOSYA", help="çıktı dosyası (varsayılan: stdout)")
    return parser


def _raw_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--raw", action="store_true", help="FASTA başlığı olmadan, kayıt başına bir satır")
    return parser


def _code_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--genetic-code", choices=GENETIC_CODES, default="standard",
                        help="genetik kod (varsayılan: standard)")
    return parser


def _input_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("inputs", nargs="*", metavar="GİRDİ",
                        help="dizi, FASTA ya da .crn dosyası, '-' stdin (varsayılan: stdin)")
    parser.add_argument("--type", choices=SEQUENCE_TYPES, help="dizi tipi (varsayılan: otomatik)")
    return parser


def _source_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("veri kaynağı")
    group.add_argument("--source", choices=("ncbi", "gencode"), default="ncbi",
                       help="veri kaynağı (varsayılan: ncbi)")
    group.add_argument("--organism", "--taxid", dest="organism", default="human", metavar="ORGANİZMA",
                       help="organizma adı ya da NCBI Taxonomy ID, ör. mouse, 9031 "
                            "(varsayılan: human); taxid bulmak için: crn taxonomy AD")
    group.add_argument("--release", help="GENCODE sürümü (ör. v50)")
    group.add_argument("--canonical", action="store_true",
                       help="yalnızca ana transkript (--source gencode)")
    return parser


def _hpo_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--release", help="HPO sürümü (ör. v2026-09-01)")
    group = parser.add_argument_group("HPO süzgeçleri")
    group.add_argument("--exclude", action="append", metavar="TERİM", help="dışlanacak terim (yinelenebilir)")
    group.add_argument("--disease", metavar="AD", help="hastalık adı ya da kimliği")
    group.add_argument("--min-frequency", type=float, metavar="ORAN",
                       help="en küçük sıklık; sıklığı bilinmeyen satırları da eler")
    group.add_argument("--association", action="append", help="ilişki türü, ör. MENDELIAN (yinelenebilir)")
    group.add_argument("--evidence", action="append", help="kanıt kodu, ör. PCS (yinelenebilir)")
    group.add_argument("--onset", metavar="TERİM", help="başlangıç yaşı terimi")
    group.add_argument("--inheritance", metavar="TERİM", help="kalıtım terimi")
    group.add_argument("--no-descendants", action="store_true", help="alt terimleri katmaz")
    return parser


class _LazyChoices:
    """``crn.fetch.clinvar`` sabitini ilk denetimde yükleyen seçenek kümesi.

    Ayrıştırıcı kurulurken ``crn.fetch`` (requests, tqdm) yüklenmez; liste yalnızca
    seçenek denetlenirken ya da yardım metni yazılırken okunur.
    """

    def __init__(self, name: str):
        self.name = name

    def _values(self) -> Seq[str]:
        from crn.fetch import clinvar
        return getattr(clinvar, self.name)

    def __contains__(self, value: object) -> bool:
        return value in self._values()

    def __iter__(self):
        return iter(self._values())


def _clinvar_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--release", help="ClinVar arşiv ayı (ör. 2026-11)")
    parser.add_argument("--assembly", choices=_LazyChoices("ASSEMBLIES"), default="GRCh38",
                        metavar="DERLEME", help="%(choices)s (varsayılan: GRCh38)")
    group = parser.add_argument_group("ClinVar süzgeçleri")
    group.add_argument("--condition", action="append", metavar="KOŞUL",
                       help="hastalık adı parçası ya da kimlik, ör. OMIM:249100 (yinelenebilir)")
    group.add_argument("--significance", action="append", choices=_LazyChoices("SIGNIFICANCES"),
                       metavar="SINIF", help="klinik anlam sınıfı (yinelenebilir): %(choices)s")
    group.add_argument("--min-stars", type=int, choices=range(5), metavar="0-4", help="en az inceleme yıldızı")
    group.add_argument("--type", action="append", metavar="TİP", help="varyant tipi (yinelenebilir)")
    group.add_argument("--origin", action="append", metavar="KÖKEN", help="germline, somatic ... (yinelenebilir)")
    group.add_argument("--effect", action="append", choices=_LazyChoices("EFFECTS"), metavar="ETKİ",
                       help="protein etkisi (yinelenebilir): %(choices)s")
    group.add_argument("--chromosome", metavar="KROMOZOM")
    group.add_argument("--start", type=int)
    group.add_argument("--end", type=int)
    group.add_argument("--max-genes", type=int, default=1, metavar="N",
                       help="en çok N gen içeren varyantlar; 0 sınırsız (varsayılan: 1)")
    return parser


def _reference_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group("referans setleri (GENLER biçiminde)")
    group.add_argument("--reference", metavar="GENLER", help="yüksek ifadeli set: Fop, CAI, CBI")
    group.add_argument("--background", metavar="GENLER", help="düşük ifadeli set: χ² optimal kodonlar")
    group.add_argument("--host", metavar="GENLER", help="konak seti: RCDI")
    group.add_argument("--genome", metavar="GENLER", help="genom seti: CPB")
    group.add_argument("--trna", metavar="DOSYA", help="GtRNAdb FASTA dosyası: tAI")
    return parser


def build_parser() -> argparse.ArgumentParser:
    """Bütün komutları içeren ayrıştırıcıyı kurar."""
    try:
        version = metadata.version("crn")
    except metadata.PackageNotFoundError:
        version = "bilinmiyor"
    parser = argparse.ArgumentParser(
        prog="crn", description="DNA, RNA ve protein dizileri, genomik veri erişimi ve kodon kullanımı.",
        epilog="Komut ayrıntısı için: crn KOMUT --help")
    parser.add_argument("--version", action="version", version=f"crn {version}")
    commands = parser.add_subparsers(dest="command", metavar="KOMUT", required=True)

    common, output, raw, code = _common_parser(), _output_parser(), _raw_parser(), _code_parser()
    inputs, source = _input_parser(), _source_parser()
    hpo_filters, references = _hpo_parser(), _reference_parser()

    def add(subs: Any, name: str, func: Any, text: str, *parents: argparse.ArgumentParser
            ) -> argparse.ArgumentParser:
        sub = subs.add_parser(name, help=text, description=text, parents=[common, *parents])
        sub.set_defaults(func=func)
        return sub

    def group(name: str, text: str) -> Any:
        sub = commands.add_parser(name, help=text, description=text)
        return sub.add_subparsers(dest=f"{name}_command", metavar="İŞLEM", required=True)

    seq = group("seq", "dizi işlemleri")
    add(seq, "info", _cmd_seq_info, "ad, tip, uzunluk ve GC yüzdesi", inputs, output)
    add(seq, "revcomp", _cmd_seq_revcomp, "ters tümleyici", inputs, output, raw)
    sub = add(seq, "transcribe", _cmd_seq_transcribe, "DNA'dan RNA", inputs, output, raw)
    sub.add_argument("--strand", choices=("coding", "template"), default="coding",
                     help="girdinin zinciri: coding yalnızca T->U yapar, template tümleyici "
                          "RNA üretir (varsayılan: coding)")
    add(seq, "reverse-transcribe", _cmd_seq_reverse_transcribe, "RNA'dan tümleyici cDNA",
        inputs, output, raw)
    add(seq, "translate", _cmd_seq_translate,
        "ilk nükleotitten itibaren çeviri; DNA kodlayan zincir sayılır", inputs, output, raw, code)
    sub = add(seq, "orfs", _cmd_seq_orfs, "altı çerçevede ORF taraması", inputs, output, raw, code)
    sub.add_argument("--min-length", type=int, default=-1, metavar="NT", help="stop dahil en küçük uzunluk")
    sub.add_argument("--alternative-starts", action="store_true", help="kodun tüm başlangıç kodonları")
    sub.add_argument("--tsv", action="store_true", help="dizi yerine koordinat tablosu")

    add(commands, "protein", _cmd_protein, "protein özellikleri tablosu", inputs, output)

    get = group("get", "tek kayıt indirme")
    for kind, text in (("dna", "genin genomik DNA'sı"), ("rna", "transkript"), ("protein", "protein")):
        sub = add(get, kind, _cmd_get, text, source, output, raw)
        sub.add_argument("ids", nargs="+", metavar="KİMLİK")
        sub.set_defaults(kind=kind)
    sub = add(get, "region", _cmd_region, "genomik bölge", source, output, raw)
    sub.add_argument("accession", metavar="DİZİ", help="kromozom ya da erişim numarası")
    sub.add_argument("start", type=int)
    sub.add_argument("end", type=int)
    sub.add_argument("--strand", choices=("+", "-"), default="+", help="zincir (varsayılan: +)")

    sub = add(commands, "taxonomy", _cmd_taxonomy,
              "organizma adından NCBI Taxonomy ID'sine ya da taxid'den organizmaya", output)
    sub.add_argument("queries", nargs="+", metavar="AD|TAXID",
                     help="bilimsel ya da yaygın ad (\"Gallus gallus\", chicken) veya taxid")
    sub.add_argument("--limit", type=int, default=10, help="ad başına en çok kayıt (varsayılan: 10)")
    sub.set_defaults(organism="human")

    sub = add(commands, "gene", _cmd_gene, "gen özeti, FASTA ya da CRN çıktısı", source, output, raw, code)
    sub.add_argument("symbols", nargs="+", metavar="GEN")
    sub.add_argument("--isoforms", action="store_true", help="izoformlara ayırır")
    sub.add_argument("--fasta", choices=("dna", "rna", "cds", "protein"), help="dizileri FASTA yazar")

    cub = group("cub", "kodon kullanımı")
    sub = add(cub, "table", _cmd_cub_table, "gen başına ölçüm tablosu", source, code, references, output)
    sub.add_argument("genes", nargs="+", metavar="GENLER")
    sub.add_argument("--kingdom", choices=("eukaryote", "prokaryote"), default="eukaryote", help="tAI kuralı")
    sub.add_argument("--no-protein", action="store_true", help="protein sütunlarını katmaz")
    sub = add(cub, "codons", _cmd_cub_codons, "kodon sayıları ve RSCU", source, code, output)
    sub.add_argument("genes", nargs="+", metavar="GENLER")
    sub.add_argument("--gene", metavar="AD", help="yalnızca bu izoform (varsayılan: hepsi birlikte)")
    sub = add(cub, "plot", _cmd_cub_plot, "grafik üretir", source, code, references)
    sub.add_argument("kind", choices=PLOT_KINDS, metavar="TÜR", help=", ".join(PLOT_KINDS))
    sub.add_argument("genes", nargs="+", metavar="GENLER")
    sub.add_argument("-o", "--output", required=True, metavar="DOSYA", help="görüntü dosyası")
    sub.add_argument("--gene", metavar="AD", help="rscu grafiğinde izoform")
    sub.add_argument("--labels", action="store_true", help="gen adlarını yazar")

    sub = add(commands, "dco", _cmd_dco, "iki gen grubunun diferansiyel kodon karşılaştırması",
              source, code, output)
    sub.add_argument("group_a", metavar="A", help="GENLER; --file ile grup adı")
    sub.add_argument("group_b", metavar="B", help="GENLER; --file ile grup adı")
    sub.add_argument("--file", metavar="DOSYA", help="iki grubu içeren CRN listesi (crn hpo pack)")
    sub.add_argument("--key", default="group", help="grup anahtarı (varsayılan: group)")
    sub.add_argument("--name-a", metavar="AD", help="A grubunun adı (--file ile grup adı kullanılır)")
    sub.add_argument("--name-b", metavar="AD", help="B grubunun adı (--file ile grup adı kullanılır)")
    sub.add_argument("--report", choices=DCO_REPORTS, default="robust",
                     help="rapor; robust Fisher ve permütasyondan geçenlerdir (varsayılan)")
    sub.add_argument("--level", choices=("family", "subfamily", "global"), default="family",
                     help="karşılaştırma düzeyi (varsayılan: family)")
    sub.add_argument("--fdr", type=float, default=0.05, metavar="Q", help="q eşiği (varsayılan: 0.05)")
    sub.add_argument("--min-log2-or", type=float, default=0.0, metavar="X", help="en küçük |log2 OR|")
    sub.add_argument("--permutations", type=int, default=5000, metavar="N",
                     help="karıştırma sayısı (varsayılan: 5000)")
    sub.add_argument("--seed", type=int, default=0, help="permütasyon tohumu")
    sub.add_argument("--plot", metavar="DOSYA", help="volkan grafiği dosyası")

    hpo = group("hpo", "semptomdan gene sorgular")
    sub = add(hpo, "search", _cmd_hpo_search, "terim arama", output)
    sub.add_argument("text", metavar="METİN")
    sub.add_argument("--limit", type=int)
    sub.add_argument("--release", help="HPO sürümü")
    for name, func, text in (("genes", _cmd_hpo_genes, "terime bağlı genler"),
                             ("symbols", _cmd_hpo_symbols, "gen sembolleri, satır başına bir")):
        sub = add(hpo, name, func, text, hpo_filters, output)
        sub.add_argument("terms", nargs="+", metavar="TERİM")
    sub = add(hpo, "phenotypes", _cmd_hpo_phenotypes, "genin fenotipleri", hpo_filters, output)
    sub.add_argument("gene", metavar="GEN")
    sub = add(hpo, "diseases", _cmd_hpo_diseases, "hastalık arama", hpo_filters, output)
    sub.add_argument("name", metavar="AD")
    sub = add(hpo, "summary", _cmd_hpo_summary, "sürüm ve kayıt sayıları", output)
    sub.add_argument("--release", help="HPO sürümü")
    sub = add(hpo, "pack", _cmd_hpo_pack, "gen gruplarını tek CRN listesine yazar (GENCODE gerekir)",
              hpo_filters)
    sub.add_argument("groups", nargs="+", metavar="AD=TERİM")
    sub.add_argument("-o", "--output", required=True, metavar="DOSYA", help=".crn dosyası")
    sub.add_argument("--key", default="group", help="grup anahtarı (varsayılan: group)")
    sub.add_argument("--overlap", choices=("exclude", "first", "error"), default="exclude",
                     help="iki gruba giren gen (varsayılan: exclude)")

    clinvar = group("clinvar", "klinik varyantlar: gen ve hastalık sorguları")
    clinvar_filters = _clinvar_parser()
    sub = add(clinvar, "variants", _cmd_clinvar_variants, "süzgece uyan varyantlar", clinvar_filters, output)
    sub.add_argument("genes", nargs="*", metavar="GEN")
    sub.add_argument("--limit", type=int)
    add(clinvar, "genes", _cmd_clinvar_genes, "gen başına varyant sayıları", clinvar_filters, output)
    add(clinvar, "symbols", _cmd_clinvar_symbols, "gen sembolleri, satır başına bir", clinvar_filters, output)
    sub = add(clinvar, "conditions", _cmd_clinvar_conditions, "hastalık arama", clinvar_filters, output)
    sub.add_argument("name", nargs="?", metavar="AD")
    sub.add_argument("--gene", action="append", metavar="GEN", help="yalnızca bu genin hastalıkları")
    sub.add_argument("--all", action="store_true", help="'not provided' gibi koşulları da gösterir")
    sub.add_argument("--limit", type=int)
    add(clinvar, "summary", _cmd_clinvar_summary, "sürüm, derleme ve kayıt sayıları", clinvar_filters, output)
    sub = add(clinvar, "pack", _cmd_clinvar_pack,
              "gen gruplarını tek CRN listesine yazar (GENCODE gerekir)", clinvar_filters)
    sub.add_argument("groups", nargs="+", metavar="AD=KOŞUL")
    sub.add_argument("-o", "--output", required=True, metavar="DOSYA", help=".crn dosyası")
    sub.add_argument("--key", default="group", help="grup anahtarı (varsayılan: group)")
    sub.add_argument("--overlap", choices=("exclude", "first", "error"), default="exclude",
                     help="iki gruba giren gen (varsayılan: exclude)")

    hgnc = group("hgnc", "gen tablosu ve ortologlar")
    add(hgnc, "columns", _cmd_hgnc_columns, "gen tablosu sütunları", output)
    for name, func, text in (("symbols", _cmd_hgnc_symbols, "süzgece uyan semboller"),
                             ("query", _cmd_hgnc_query, "gen tablosu sorgusu"),
                             ("orthologs", _cmd_hgnc_orthologs, "HCOP ortolog sorgusu")):
        sub = add(hgnc, name, func, text, output)
        sub.add_argument("filters", nargs="*", metavar="SÜTUN=DEĞER")
        if name != "symbols":
            sub.add_argument("--fields", metavar="A,B", help="gösterilecek sütunlar")

    sub = add(commands, "update", _cmd_update, "yerel veri deposunu hazırlar")
    sub.add_argument("target", choices=("ncbi", "gencode", "hpo", "hgnc", "clinvar"))
    sub.add_argument("genes", nargs="*", metavar="GEN", help="ncbi için gen adları")
    sub.add_argument("--organism", "--taxid", dest="organism", default="human", metavar="ORGANİZMA",
                     help="ncbi ve gencode için organizma adı ya da NCBI Taxonomy ID (varsayılan: human)")
    sub.add_argument("--species", metavar="TÜR,TÜR",
                     help="hgnc için ayrıca hazırlanacak HCOP türleri; ad ya da taxid, ör. chicken,7955")
    sub.add_argument("--release", help="GENCODE, HPO ya da ClinVar sürümü")
    sub.add_argument("--assembly", choices=_LazyChoices("ASSEMBLIES"), default="GRCh38",
                     metavar="DERLEME", help="ClinVar derlemesi: %(choices)s (varsayılan: GRCh38)")
    sub.add_argument("--force", action="store_true", help="yeniden indirir")
    return parser


# --- giriş noktası -------------------------------------------------------------

def _show_warning(message: Any, category: Any, filename: str, lineno: int,
                  file: Any = None, line: Optional[str] = None) -> None:
    _log(f"uyarı: {message}")


def main(argv: Optional[Seq[str]] = None) -> int:
    """Komutu çalıştırır ve çıkış kodunu döndürür."""
    args = build_parser().parse_args(argv)
    try:
        with warnings.catch_warnings():
            warnings.showwarning = _show_warning
            args.func(args)
        sys.stdout.flush()
    except BrokenPipeError:
        # Okuyucu erken kapandı (ör. `| head`); yorumlayıcı kapanışta yeniden denemesin.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
    except KeyboardInterrupt:
        return 130
    except (KeyError, ValueError, OSError) as error:
        if args.debug:
            raise
        _log(f"hata: {error.args[0] if isinstance(error, KeyError) and error.args else error}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
