"""HGVS nomenclature (version 21.1) for variants in one coding sequence of a genome.

A variant is always described on public reference sequences:

* ``NC_012920.1:m.8993T>G`` on the mitochondrial genome (``m.``, plus strand),
* ``NC_012920.1(MT-ATP6):c.467T>G`` on the coding DNA annotated in that genome - by
  exception HGVS lets a gene symbol specify the annotated segment of a mitochondrial
  genome, ``c.1`` being the first nucleotide of the start codon,
* ``YP_003024031.1:p.(Leu156Arg)`` for the predicted protein change.

Both DNA coordinate systems come from the GenBank CDS feature. A ``complement(a..b)``
location is read in transcription order: ``c.1`` is genome position ``b`` and ``c.``
numbering runs towards ``a``, while the ``m.`` description stays on the plus strand and
spells the complementary bases.

Every description is normalised as HGVS requires:

* deletions, duplications and insertions are shifted to the most 3' position of their
  own reference sequence - the plus strand for ``m.``, the coding strand for ``c.``;
* an insertion copying the sequence directly 5' of it is a duplication, two or more
  extra tandem copies are a repeat (in a ``c.`` coding sequence only for 3n units);
* changes touching each other, or two changes one nucleotide apart inside one codon,
  are merged into one deletion-insertion; a replaced nucleotide is a substitution and a
  stretch replaced by its reverse complement an inversion.

Mitochondrial genomes are circular, so the coordinates are kept on a rotated view of the
molecule centred on the gene: a 3' shift can then cross the origin like on the real
molecule (``m.16569_1insT``).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from typing import Callable, Iterable

from Bio.SeqFeature import ExactPosition

from .cds_annotation import TranslationException, get_codon_table, read_cds_annotation
from .hgvs_protein import (
    FRAMESHIFT,
    START,
    CdsEdit,
    ProteinConsequence,
    predict_consequence,
    three_letter,
    translation_end,
)

HGVS_VERSION = "21.1"
GENOME_PREFIX = "m"
CODING_PREFIX = "c"
SUPPORTED_PREFIXES = (CODING_PREFIX, GENOME_PREFIX)

SUBSTITUTION = "substitution"
DELETION = "deletion"
DUPLICATION = "duplication"
INSERTION = "insertion"
DELINS = "deletion-insertion"
INVERSION = "inversion"
REPEAT = "repeat"

_COMPLEMENT = str.maketrans("ACGT", "TGCA")

_SEQUENCE = r"[ACGT]+"
_M_POSITION = r"[1-9][0-9]*"
_C_POSITION = r"[-*]?[1-9][0-9]*"
_REFERENCE = re.compile(
    r"(?:(?P<accession>[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)?\.[0-9]+)"
    r"(?:\((?P<gene>[^()\s]+)\))?:)?(?P<prefix>[A-Za-z])\.(?P<body>.*)",
    re.DOTALL,
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


class HgvsError(ValueError):
    """A description that is not valid HGVS or does not fit the reference sequence."""


@dataclass(frozen=True)
class Edit:
    """Replace ``seq[start:end]`` of a frame with ``alt``; ``start == end`` inserts."""

    start: int
    end: int
    alt: str


@dataclass(frozen=True)
class Frame:
    """A linear, oriented view of the genome.

    ``seq[i]`` is the base at genome index ``to_genome(i)`` (0-based), complemented in a
    reversed frame. A circular genome is rotated so the gene sits in the middle.
    """

    seq: str
    origin: int
    reverse: bool
    circular: bool

    def to_genome(self, index: int) -> int:
        position = self.origin - index if self.reverse else self.origin + index
        return position % len(self.seq) if self.circular else position

    def from_genome(self, position: int) -> int:
        index = self.origin - position if self.reverse else position - self.origin
        return index % len(self.seq) if self.circular else index


def _build_frame(genome: str, circular: bool, reverse: bool, anchor: int, anchor_index: int) -> Frame:
    """Build a frame where genome index ``anchor`` sits at frame index ``anchor_index``."""
    size = len(genome)
    if not circular:
        if reverse:
            return Frame(reverse_complement(genome), size - 1, True, False)
        return Frame(genome, 0, False, False)

    if reverse:
        origin = (anchor + anchor_index) % size
        cut = (origin + 1) % size
        return Frame(reverse_complement(genome[cut:] + genome[:cut]), origin, True, True)

    origin = (anchor - anchor_index) % size
    return Frame(genome[origin:] + genome[:origin], origin, False, True)


def _convert(edit: Edit, source: Frame, target: Frame) -> Edit:
    """Express an edit of one frame in another frame of the same genome."""
    if source.reverse == target.reverse:
        start = target.from_genome(source.to_genome(edit.start))
        return Edit(start, start + edit.end - edit.start, edit.alt)

    if edit.end > edit.start:
        start = target.from_genome(source.to_genome(edit.end - 1))
        return Edit(start, start + edit.end - edit.start, reverse_complement(edit.alt))

    # An insertion between source bases start-1 and start sits, read backwards, between
    # the target bases of start and start-1.
    point = target.from_genome(source.to_genome(edit.start - 1))
    return Edit(point, point, reverse_complement(edit.alt))


@dataclass(frozen=True)
class ReferenceContext:
    """One coding sequence of a GenBank record together with its HGVS coordinates."""

    accession: str
    gene: str
    symbol: str
    protein_id: str | None
    aliases: tuple[str, ...]
    strand: int
    circular: bool
    genome_length: int
    is_coding: bool
    cds_complete: bool
    transl_table: int
    codon_start: int
    transl_except: tuple[TranslationException, ...]
    plus: Frame = field(repr=False)
    oriented: Frame = field(repr=False)
    cds_offset: int
    cds_length: int

    @classmethod
    def from_genbank(cls, record, feature) -> "ReferenceContext":
        """Read the context of ``feature`` (normally a CDS) of a GenBank ``record``."""
        genome = str(record.seq).upper()
        size = len(genome)
        circular = str(record.annotations.get("topology", "")).lower() == "circular"
        gene = (feature.qualifiers.get("gene") or [feature.qualifiers.get("locus_tag", ["?"])[0]])[0]
        positions = [int(position) for position in feature.location]
        if not positions:
            raise HgvsError(f"The {gene} feature has an empty location")

        strands = {part.strand for part in feature.location.parts}
        if len(strands) > 1:
            raise HgvsError(f"The {gene} feature spans both strands, which is not supported")
        strand = -1 if strands == {-1} else 1

        for previous, current in zip(positions, positions[1:]):
            step = current - previous
            contiguous = step % size == strand % size if circular else step == strand
            if not contiguous:
                raise HgvsError(
                    f"The {gene} location {feature.location} is not one contiguous stretch "
                    "of the genome; joined (spliced) locations are not supported"
                )

        length = len(positions)
        padding = (size - length) // 2
        plus_start = positions[0] if strand == 1 else positions[-1]
        plus = _build_frame(genome, circular, False, plus_start, padding)
        oriented = _build_frame(genome, circular, strand == -1, positions[0], padding)
        cds_offset = oriented.from_genome(positions[0])

        expected = str(feature.extract(record.seq)).upper()
        if oriented.seq[cds_offset:cds_offset + length] != expected:
            raise HgvsError(f"Could not map the {gene} location {feature.location} onto the genome")

        annotation = read_cds_annotation(feature)
        return cls(
            accession=record.id,
            gene=gene,
            symbol=_official_symbol(record, gene),
            protein_id=(feature.qualifiers.get("protein_id") or [None])[0],
            aliases=_gene_aliases(record, feature, gene),
            strand=strand,
            circular=circular,
            genome_length=size,
            is_coding=feature.type == "CDS",
            cds_complete=all(
                type(part.start) is ExactPosition and type(part.end) is ExactPosition
                for part in feature.location.parts
            ),
            transl_table=annotation.transl_table,
            codon_start=annotation.codon_start,
            transl_except=annotation.transl_except,
            plus=plus,
            oriented=oriented,
            cds_offset=cds_offset,
            cds_length=length,
        )

    @property
    def cds_sequence(self) -> str:
        return self.oriented.seq[self.cds_offset:self.cds_offset + self.cds_length]

    @property
    def c_available(self) -> bool:
        """``c.`` numbering needs a complete CDS whose first base starts the start codon."""
        return self.is_coding and self.cds_complete and self.codon_start == 1

    @property
    def c_unavailable_reason(self) -> str:
        if not self.is_coding:
            return f"{self.gene} is not a protein coding sequence"
        if not self.cds_complete:
            return f"the {self.gene} CDS is annotated as incomplete (partial 5' or 3' end)"
        return f"the {self.gene} CDS does not start with a complete codon"

    @property
    def first_genome_position(self) -> int:
        """1-based genome position of the first CDS base (``c.1``)."""
        return self.oriented.to_genome(self.cds_offset) + 1

    @property
    def last_genome_position(self) -> int:
        return self.oriented.to_genome(self.cds_offset + self.cds_length - 1) + 1

    @property
    def m_reference(self) -> str:
        return f"{self.accession}:{GENOME_PREFIX}."

    @property
    def c_reference(self) -> str:
        return f"{self.accession}({self.symbol}):{CODING_PREFIX}."

    def genome_position(self, cds_index: int) -> int:
        """1-based genome position of a 0-based CDS index."""
        return self.oriented.to_genome(self.cds_offset + cds_index) + 1

    def m_label(self, index: int) -> str:
        return str(self.plus.to_genome(index) + 1)

    def c_label(self, index: int) -> str:
        number = index - self.cds_offset + 1
        if number > self.cds_length:
            return f"*{number - self.cds_length}"
        if number < 1:
            return f"-{1 - number}"
        return str(number)

    def codon_number(self, index: int) -> int | None:
        """0-based codon of an oriented frame index inside the CDS, ``None`` outside."""
        cds_index = index - self.cds_offset - (self.codon_start - 1)
        if cds_index < 0 or index >= self.cds_offset + self.cds_length:
            return None
        return cds_index // 3


def _gene_aliases(record, feature, gene: str) -> tuple[str, ...]:
    """Names a description may use for the gene: GenBank name, synonyms, official symbol."""
    names = [gene]
    for synonyms in feature.qualifiers.get("gene_synonym", []):
        names.extend(name.strip() for name in synonyms.split(";"))
    for candidate in record.features:
        if candidate.type == "gene" and (candidate.qualifiers.get("gene") or [None])[0] == gene:
            for note in candidate.qualifiers.get("nomenclature", []):
                match = re.search(r"Official Symbol:\s*([^\s|]+)", note)
                if match:
                    names.append(match.group(1))
    names.append(_official_symbol(record, gene))
    return tuple(dict.fromkeys(name for name in names if name))


def _official_symbol(record, gene: str) -> str:
    """The symbol of the nomenclature authority (HGNC, MGI, ZFIN ...) when annotated.

    A symbol containing characters with a meaning in HGVS (FlyBase uses ``mt:ND1``) would
    break the ``reference:description`` structure, so the GenBank gene name is used then.
    """
    for feature in record.features:
        if feature.type != "gene" or (feature.qualifiers.get("gene") or [None])[0] != gene:
            continue
        for note in feature.qualifiers.get("nomenclature", []):
            match = re.search(r"Official Symbol:\s*([^\s|]+)", note)
            if match and re.fullmatch(r"[A-Za-z0-9._-]+", match.group(1)):
                return match.group(1)
    return gene


# --------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class VariantDescription:
    """One normalised variant with its HGVS descriptions and its place in the CDS."""

    inputs: tuple[int, ...]
    texts: tuple[str, ...]
    kind: str
    m_edits: tuple[str, ...]
    c_edits: tuple[str, ...] | None
    m: str
    c: str | None
    cds_start: int
    cds_end: int
    ref: str
    alt: str
    protein: ProteinConsequence | None = None
    p: str | None = None

    @property
    def length_change(self) -> int:
        return len(self.alt) - len(self.ref)

    @property
    def frameshift(self) -> bool:
        return self.length_change % 3 != 0

    @property
    def consequence(self) -> str:
        return self.protein.label if self.protein else "Not analysed"

    def to_dict(self) -> dict:
        return {
            "inputs": list(self.inputs),
            "type": self.kind,
            "m": self.m,
            "c": self.c,
            "p": self.p,
            "consequence": self.consequence,
            "proteinKind": self.protein.kind if self.protein else None,
            "cdsStart": self.cds_start,
            "cdsEnd": self.cds_end,
            "ref": self.ref,
            "alt": self.alt,
            "lengthChange": self.length_change,
            "frameshift": self.frameshift,
        }


@dataclass(frozen=True)
class DescriptionError:
    inputs: tuple[int, ...]
    texts: tuple[str, ...]
    message: str

    def to_dict(self) -> dict:
        return {"inputs": list(self.inputs), "texts": list(self.texts), "message": self.message}


@dataclass(frozen=True)
class DescriptionResult:
    variants: tuple[VariantDescription, ...]
    errors: tuple[DescriptionError, ...]
    allele_m: str | None = None
    allele_c: str | None = None
    allele_p: str | None = None

    def raise_for_errors(self) -> None:
        if self.errors:
            raise HgvsError("; ".join(error.message for error in self.errors))

    def to_dict(self) -> dict:
        return {
            "variants": [variant.to_dict() for variant in self.variants],
            "errors": [error.to_dict() for error in self.errors],
            "allele": {"m": self.allele_m, "c": self.allele_c, "p": self.allele_p},
            "hgvsVersion": HGVS_VERSION,
        }


# --------------------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Parsed:
    index: int
    text: str
    edit: Edit  # in the oriented (coding) frame


def _edit_patterns(position: str) -> list[tuple[str, re.Pattern]]:
    single = rf"(?P<start>{position})"
    pair = rf"(?P<start>{position})_(?P<end>{position})"
    either = rf"(?P<start>{position})(?:_(?P<end>{position}))?"
    return [
        ("sub", re.compile(rf"{single}(?P<ref>[ACGT])>(?P<alt>[ACGT])")),
        ("delins", re.compile(rf"{either}delins(?P<alt>{_SEQUENCE})")),
        ("del", re.compile(rf"{either}del(?P<ref>{_SEQUENCE})?")),
        ("dup", re.compile(rf"{either}dup(?P<ref>{_SEQUENCE})?")),
        ("ins", re.compile(rf"{pair}ins(?P<alt>{_SEQUENCE})")),
        ("inv", re.compile(rf"{pair}inv")),
        ("repeat", re.compile(rf"{either}(?P<unit>{_SEQUENCE})\[(?P<copies>[1-9][0-9]*)\]")),
    ]


_PATTERNS = {GENOME_PREFIX: _edit_patterns(_M_POSITION), CODING_PREFIX: _edit_patterns(_C_POSITION)}
_LENIENT = [
    (kind, re.compile(pattern.pattern.replace("[ACGT]", "[A-Za-z]").replace("[1-9]", "[0-9]"), re.I))
    for kind, pattern in _edit_patterns(r"[-*]?[0-9]+")
]
_EXAMPLES = (
    "substitution c.10A>G, deletion c.10del or c.10_12del, duplication c.10dup, "
    "insertion c.10_11insA, deletion-insertion c.10_12delinsTT, inversion c.10_13inv, "
    "repeat c.10_12CAG[4]"
)


def _split_allele(body: str) -> list[str]:
    """Split ``[a;b]`` (variants in cis) into its edits; other forms stay as they are."""
    if "(;)" in body:
        raise HgvsError(
            "Variants with an unknown phase '(;)' cannot be analysed together; "
            "list them as one allele, e.g. m.[100A>G;200del]"
        )
    if not body.startswith("["):
        return [body]

    depth = 0
    for index, character in enumerate(body):
        depth += {"[": 1, "]": -1}.get(character, 0)
        if depth == 0:
            if index != len(body) - 1:
                raise HgvsError(
                    "Only one allele (variants in cis, e.g. m.[100A>G;200del]) can be "
                    "analysed; descriptions of several alleles '[...];[...]' are not supported"
                )
            break
    else:
        raise HgvsError("Unbalanced '[' in the allele description")

    edits, depth, current = [], 0, ""
    for character in body[1:-1]:
        depth += {"[": 1, "]": -1}.get(character, 0)
        if character == ";" and depth == 0:
            edits.append(current)
            current = ""
        else:
            current += character
    edits.append(current)
    if any(not edit for edit in edits):
        raise HgvsError("An allele description contains an empty variant")
    return edits


def _explain_invalid_edit(edit: str, prefix: str) -> str:
    label = f"{prefix}.{edit}"
    if not edit:
        return "The variant description is empty"
    if re.search(r"\s", edit):
        return "Spaces are not permitted in HGVS descriptions"
    if "=" in edit:
        return f"'{label}' describes a sequence found unchanged ('='), there is no variant to analyse"
    if "?" in edit or "(" in edit:
        return f"'{label}': uncertain positions or sequences are not supported"
    if prefix == CODING_PREFIX and re.search(r"[0-9][+-][0-9]", edit):
        return f"'{label}': intronic positions do not exist in a mitochondrial coding sequence"
    if re.search(r"(?<![0-9])0", edit):
        return f"'{label}': positions start at 1 and are written without leading zeros"
    if re.search(r"[0-9]ins", edit) and not re.search(r"_[-*]?[0-9]+ins", edit):
        return f"'{label}': an insertion lies between two flanking positions, e.g. {prefix}.10_11insA"
    if re.search(r"(del|dup)[0-9]", edit):
        return f"'{label}': give the first and last position of the range, not its length"
    if "con" in edit:
        return f"'{label}': 'con' is no longer used, describe a conversion as a deletion-insertion"
    for _, pattern in _LENIENT:
        match = pattern.fullmatch(edit)
        if not match:
            continue
        sequences = "".join(
            value for value in match.groupdict().values() if value and value[0].isalpha()
        )
        if re.search(r"[BDHKMNRSVWY]", sequences):
            return f"'{label}': ambiguous nucleotides (N, R, Y, ...) are not supported"
        return (
            f"'{label}': HGVS descriptions are case-sensitive - use 'del', 'dup', 'ins', "
            "'delins' and 'inv' in lower case and nucleotides in upper case (A, C, G, T)"
        )
    return f"'{label}' is not a supported HGVS DNA variant; supported are {_EXAMPLES}"


class _Parser:
    def __init__(self, context: ReferenceContext):
        self.context = context

    def parse(self, index: int, text: str) -> list[_Parsed]:
        description = text.strip()
        if not description:
            raise HgvsError("The variant description is empty")
        if re.search(r"\s", description):
            raise HgvsError(f"'{description}': spaces are not permitted in HGVS descriptions")

        match = _REFERENCE.fullmatch(description)
        if match is None:
            if re.match(r"[A-Z]+_?[0-9]+[:(]", description):
                raise HgvsError(
                    f"'{description}': the reference sequence needs its version number, "
                    f"e.g. {self.context.accession}"
                )
            raise HgvsError(
                f"Unsupported mutation notation: '{description}'. Every variant needs the "
                f"'{CODING_PREFIX}.' (coding DNA) or '{GENOME_PREFIX}.' (mitochondrial genome) "
                "prefix, e.g. 'm.8993T>G' or 'c.10_12del'"
            )

        prefix = match.group("prefix")
        if prefix.lower() in SUPPORTED_PREFIXES and prefix.isupper():
            raise HgvsError(f"'{description}': the '{prefix.lower()}.' prefix is written in lower case")
        if prefix not in SUPPORTED_PREFIXES:
            raise HgvsError(
                f"Unsupported coordinate prefix '{prefix}.' in '{description}'. Only "
                f"'{GENOME_PREFIX}.' (mitochondrial genome) and '{CODING_PREFIX}.' (coding DNA) "
                "descriptions can be analysed"
            )

        self._check_reference(description, prefix, match.group("accession"), match.group("gene"))
        if prefix == CODING_PREFIX and not self.context.c_available:
            raise HgvsError(
                f"'{description}': c. numbering is not available because "
                f"{self.context.c_unavailable_reason}; use m. positions instead"
            )

        return [
            _Parsed(index, description, self._parse_edit(edit, prefix))
            for edit in _split_allele(match.group("body"))
        ]

    def _check_reference(self, description: str, prefix: str, accession, gene) -> None:
        context = self.context
        if accession is None:
            return
        if accession != context.accession:
            raise HgvsError(
                f"'{description}' refers to {accession}, but the selected genome is "
                f"{context.accession}"
            )
        if gene is not None and gene not in context.aliases:
            raise HgvsError(
                f"'{description}' refers to the gene {gene}, but the selected gene is "
                f"{context.symbol}"
            )
        if prefix == CODING_PREFIX and gene is None:
            raise HgvsError(
                f"'{description}': {accession} is a genomic reference, so a c. description "
                f"has to name the gene, e.g. {context.c_reference}1A>G"
            )

    def _index(self, token: str, prefix: str) -> int:
        context = self.context
        # No genome has positions with this many digits; it also keeps int() cheap.
        if len(token.lstrip("-*")) > 12:
            raise HgvsError(f"{prefix}.{token[:15]}... is outside of the reference sequence {context.accession}")
        if prefix == GENOME_PREFIX:
            number = int(token)
            if number > context.genome_length:
                raise HgvsError(
                    f"m.{number} is outside of the {context.genome_length} nt reference "
                    f"sequence {context.accession}"
                )
            index = context.plus.from_genome(number - 1)
            if not 0 <= index < len(context.plus.seq):
                raise HgvsError(f"m.{number} is outside of the reference sequence")
            return index

        if token.startswith("*"):
            index = context.cds_offset + context.cds_length - 1 + int(token[1:])
        elif token.startswith("-"):
            index = context.cds_offset - int(token[1:])
        else:
            number = int(token)
            if number > context.cds_length:
                raise HgvsError(
                    f"c.{number} is outside of the {context.symbol} coding sequence "
                    f"(c.1-c.{context.cds_length}); positions after the stop codon are "
                    "written c.*1, c.*2, ..."
                )
            index = context.cds_offset + number - 1
        if not 0 <= index < len(context.oriented.seq):
            raise HgvsError(f"c.{token} is outside of the reference sequence {context.accession}")
        return index

    def _parse_edit(self, edit: str, prefix: str) -> Edit:
        for kind, pattern in _PATTERNS[prefix]:
            match = pattern.fullmatch(edit)
            if match:
                break
        else:
            raise HgvsError(_explain_invalid_edit(edit, prefix))

        frame = self.context.plus if prefix == GENOME_PREFIX else self.context.oriented
        seq = frame.seq
        label = f"{prefix}.{edit}"
        start = self._index(match.group("start"), prefix)
        end = self._index(match.group("end"), prefix) if match.groupdict().get("end") else start

        if kind == "ins":
            if end != start + 1:
                raise HgvsError(
                    f"'{label}': an insertion lies between two flanking (adjacent) positions "
                    f"listed 5' to 3', e.g. {prefix}.10_11insA"
                )
        elif match.groupdict().get("end") is not None and end <= start:
            raise HgvsError(
                f"'{label}': a range lists two different positions from 5' to 3', "
                f"e.g. {prefix}.10_12del"
            )

        inserted = match.groupdict().get("alt") or ""
        if len(inserted) > self.context.genome_length:
            raise HgvsError(
                f"'{label}': the inserted sequence is longer than the whole reference sequence"
            )

        reference = seq[start:end + 1]
        if kind == "sub":
            if seq[start] != match.group("ref"):
                raise HgvsError(
                    f"'{label}': the reference nucleotide at {prefix}.{match.group('start')} "
                    f"is {seq[start]}, not {match.group('ref')}"
                )
            if match.group("ref") == match.group("alt"):
                raise HgvsError(f"'{label}' does not change the sequence")
            parsed = Edit(start, start + 1, match.group("alt"))
        elif kind in {"del", "dup"}:
            given = match.group("ref")
            if given is not None and given != reference:
                raise HgvsError(
                    f"'{label}': the reference sequence there is {reference}, not {given} "
                    "(HGVS recommends not to list the sequence at all)"
                )
            parsed = Edit(start, end + 1, "") if kind == "del" else Edit(end + 1, end + 1, reference)
        elif kind == "ins":
            parsed = Edit(end, end, match.group("alt"))
        elif kind == "delins":
            parsed = Edit(start, end + 1, match.group("alt"))
        elif kind == "inv":
            parsed = Edit(start, end + 1, reverse_complement(reference))
        else:
            parsed = self._repeat(label, seq, start, end, match, self.context.genome_length)
            unit = match.group("unit")
            if prefix == CODING_PREFIX and len(unit) % 3 and self._in_cds(start, parsed.end):
                raise HgvsError(
                    f"'{label}': in a coding sequence a repeat is only described for units "
                    "of 3, 6, 9 ... nucleotides; describe it as a duplication or an insertion"
                )

        if prefix == GENOME_PREFIX:
            return _convert(parsed, self.context.plus, self.context.oriented)
        return parsed

    def _in_cds(self, start: int, end: int) -> bool:
        """Whether the oriented frame stretch ``[start, end)`` overlaps the CDS."""
        first = self.context.cds_offset
        return start < first + self.context.cds_length and end > first

    @staticmethod
    def _repeat(label: str, seq: str, start: int, end: int, match, limit: int) -> Edit:
        unit, copies_text = match.group("unit"), match.group("copies")
        size = len(unit)
        # Checked before building anything: [N] could otherwise ask for gigabytes.
        if len(copies_text) > len(str(limit)) or int(copies_text) * size > limit:
            raise HgvsError(
                f"'{label}': {copies_text} copies of {unit} are longer than the whole "
                "reference sequence"
            )
        copies = int(copies_text)
        if match.group("end") is not None:
            tract = seq[start:end + 1]
            if len(tract) % size or tract != unit * (len(tract) // size):
                raise HgvsError(
                    f"'{label}': the reference sequence there is {tract}, which is not a "
                    f"whole number of {unit} units"
                )
            present = len(tract) // size
        else:
            present = 0
            while seq[start + present * size:start + (present + 1) * size] == unit:
                present += 1
            if not present:
                raise HgvsError(f"'{label}': the reference sequence does not start with {unit}")

        tract_end = start + present * size
        if copies == present:
            raise HgvsError(f"'{label}' describes the reference number of {unit} copies")
        if copies > present:
            return Edit(tract_end, tract_end, unit * (copies - present))
        return Edit(tract_end - (present - copies) * size, tract_end, "")


# --------------------------------------------------------------------------------------
# Normalisation
# --------------------------------------------------------------------------------------


def _trim(seq: str, edit: Edit) -> Edit:
    reference, alt, start = seq[edit.start:edit.end], edit.alt, edit.start
    prefix = 0
    while prefix < len(reference) and prefix < len(alt) and reference[prefix] == alt[prefix]:
        prefix += 1
    reference, alt, start = reference[prefix:], alt[prefix:], start + prefix
    suffix = 0
    while suffix < len(reference) and suffix < len(alt) and reference[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    if suffix:
        reference, alt = reference[:-suffix], alt[:-suffix]
    return Edit(start, start + len(reference), alt)


def _roll_right(seq: str, edit: Edit, limit: int | None = None) -> Edit:
    """Shift a deletion or an insertion 3' (right) while it describes the same change.

    The edit never reaches past ``limit``, the start of the next change of the allele.
    """
    limit = len(seq) if limit is None else limit
    if edit.start == edit.end:
        point, inserted = edit.start, edit.alt
        while point < limit and seq[point] == inserted[0]:
            inserted = inserted[1:] + seq[point]
            point += 1
        return Edit(point, point, inserted)

    start, end = edit.start, edit.end
    while end < limit and seq[start] == seq[end]:
        start, end = start + 1, end + 1
    return Edit(start, end, "")


def _roll_left(seq: str, edit: Edit, limit: int = 0) -> Edit:
    """Shift a deletion or an insertion 5' (left), never before ``limit``."""
    if edit.start == edit.end:
        point, inserted = edit.start, edit.alt
        while point > limit and seq[point - 1] == inserted[-1]:
            inserted = seq[point - 1] + inserted[:-1]
            point -= 1
        return Edit(point, point, inserted)

    start, end = edit.start, edit.end
    while start > limit and seq[start - 1] == seq[end - 1]:
        start, end = start - 1, end - 1
    return Edit(start, end, "")


def _smallest_unit(sequence: str) -> str:
    for size in range(1, len(sequence) + 1):
        if len(sequence) % size == 0 and sequence == sequence[:size] * (len(sequence) // size):
            return sequence[:size]
    return sequence


@dataclass(frozen=True)
class _Normalised:
    kind: str
    start: int
    end: int
    alt: str = ""
    copies: int = 0


def _classify(
    seq: str, edit: Edit, lower: int, repeat_allowed: Callable[[str, int, int], bool]
) -> _Normalised:
    """HGVS type of a trimmed edit already shifted 3'.

    A duplicated or repeated copy is a stretch of the reference, so it has to lie at or
    after ``lower``, the end of the previous change of the allele.
    """
    reference = seq[edit.start:edit.end]
    if reference and edit.alt:
        if len(reference) == 1 and len(edit.alt) == 1:
            return _Normalised(SUBSTITUTION, edit.start, edit.end, edit.alt)
        if len(reference) == len(edit.alt) and edit.alt == reverse_complement(reference):
            return _Normalised(INVERSION, edit.start, edit.end)
        return _Normalised(DELINS, edit.start, edit.end, edit.alt)

    if reference:
        return _Normalised(DELETION, edit.start, edit.end)

    point, inserted = edit.start, edit.alt
    size = len(inserted)
    if point - size >= lower and seq[point - size:point] == inserted:
        return _Normalised(DUPLICATION, point - size, point)

    unit = _smallest_unit(inserted)
    extra = size // len(unit)
    if extra >= 2 and point - len(unit) >= lower and seq[point - len(unit):point] == unit:
        tract_start = point - len(unit)
        while (
            tract_start - len(unit) >= lower
            and seq[tract_start - len(unit):tract_start] == unit
        ):
            tract_start -= len(unit)
        if repeat_allowed(unit, tract_start, point):
            present = (point - tract_start) // len(unit)
            return _Normalised(REPEAT, tract_start, point, unit, present + extra)

    return _Normalised(INSERTION, point, point, inserted)


def _format(normalised: _Normalised, seq: str, label: Callable[[int], str]) -> str:
    start, end = normalised.start, normalised.end

    def span(first: int, last: int) -> str:
        return label(first) if first == last else f"{label(first)}_{label(last)}"

    if normalised.kind == SUBSTITUTION:
        return f"{label(start)}{seq[start]}>{normalised.alt}"
    if normalised.kind == DELETION:
        return f"{span(start, end - 1)}del"
    if normalised.kind == DUPLICATION:
        return f"{span(start, end - 1)}dup"
    if normalised.kind == INVERSION:
        return f"{span(start, end - 1)}inv"
    if normalised.kind == DELINS:
        return f"{span(start, end - 1)}delins{normalised.alt}"
    if normalised.kind == INSERTION:
        return f"{label(start - 1)}_{label(start)}ins{normalised.alt}"
    return f"{span(start, end - 1)}{normalised.alt}[{normalised.copies}]"


# --------------------------------------------------------------------------------------
# Variant sets
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Piece:
    """A change in one frame and the input pieces (by number) it stands for."""

    edit: Edit
    members: frozenset[int]


class _AlleleError(HgvsError):
    """An error of some input pieces; the other pieces can still be described."""

    def __init__(self, members, message: str):
        super().__init__(message)
        self.members = frozenset(members)


def _conflict(first: Edit, second: Edit) -> bool:
    """Whether two edits of the reference touch the same bases in an ambiguous way."""
    first_insertion, second_insertion = first.start == first.end, second.start == second.end
    if first_insertion and second_insertion:
        return first.start == second.start
    if first_insertion:
        return second.start < first.start < second.end
    if second_insertion:
        return first.start < second.start < first.end
    return first.start < second.end and second.start < first.end


def _combine(seq: str, edits: Iterable[Edit]) -> Edit:
    """One edit replacing the whole stretch that non-overlapping edits cover.

    ``edits`` are given 5' to 3'; two insertions at one position keep that order.
    """
    edits = list(edits)
    start = min(edit.start for edit in edits)
    end = max(edit.end for edit in edits)
    region = seq[start:end]
    # Right to left, a deletion before an insertion at the same position and the later of
    # two insertions at one position first, so the earlier one ends up in front of it.
    ordered = sorted(
        enumerate(edits), key=lambda item: (item[1].start, item[1].end, item[0]), reverse=True
    )
    for _, edit in ordered:
        region = region[: edit.start - start] + edit.alt + region[edit.end - start :]
    return Edit(start, end, region)


def _is_fixed(edit: Edit) -> bool:
    """A substitution, deletion-insertion or inversion cannot be shifted."""
    return edit.end > edit.start and bool(edit.alt)


def _same_codon(first: Edit, second: Edit, codon_of: Callable[[int], int | None]) -> bool:
    """Two changes one nucleotide apart that together change one amino acid."""
    if not (_is_fixed(first) and _is_fixed(second)) or second.start - first.end != 1:
        return False
    codon = codon_of(first.start)
    return codon is not None and codon == codon_of(second.end - 1)


def _decompose(seq: str, edit: Edit, codon_of: Callable[[int], int | None]) -> list[Edit]:
    """Split a deletion-insertion into the changes HGVS describes individually.

    Replaced stretches of equal length are split at unchanged nucleotides: two variants
    separated by one or more nucleotides are described individually, except two changes
    one nucleotide apart inside one codon. A reverse-complemented stretch is one
    inversion, and stretches of different length stay one deletion-insertion, the form
    HGVS recommends for complex changes.
    """
    reference = seq[edit.start:edit.end]
    if not (reference and edit.alt) or len(reference) != len(edit.alt):
        return [edit]
    if len(reference) > 1 and edit.alt == reverse_complement(reference):
        return [edit]

    runs: list[list[int]] = []
    for offset, (old, new) in enumerate(zip(reference, edit.alt)):
        if old == new:
            continue
        if runs and runs[-1][1] == offset:
            runs[-1][1] = offset + 1
        else:
            runs.append([offset, offset + 1])

    joined: list[list[int]] = []
    for run in runs:
        if joined:
            gap_start, gap_end = edit.start + joined[-1][1], edit.start + run[0]
            first_codon = codon_of(edit.start + joined[-1][0])
            if (
                gap_end - gap_start == 1
                and first_codon is not None
                and first_codon == codon_of(edit.start + run[1] - 1)
            ):
                joined[-1][1] = run[1]
                continue
        joined.append(list(run))

    return [
        Edit(edit.start + first, edit.start + last, edit.alt[first:last]) for first, last in joined
    ]


def _can_touch(seq: str, first: Edit, second: Edit) -> bool:
    """Whether two changes, in this order, touch in some equivalent placement.

    A deletion or insertion describes the same change anywhere along its repeat, so it is
    shifted 3' towards the next change, and the next one 5' towards it.
    """
    if not _is_fixed(first):
        first = _roll_right(seq, first, second.start)
    if not _is_fixed(second):
        second = _roll_left(seq, second, first.end)
    return first.end >= second.start


def _normalise_allele(
    seq: str, pieces: list[_Piece], codon_of: Callable[[int], int | None]
) -> list[_Piece]:
    """Describe the changes of one allele the way HGVS does, in one frame.

    Changes that touch in some equivalent placement change consecutive nucleotides, so
    they are one change (deletion-insertion, or simpler once trimmed): ``AAA`` to ``GA``
    is ``c.4_5delinsG`` however the input placed the deletion. So are two substitutions
    one nucleotide apart inside one codon. Pairs are merged from 5' to 3', so the result
    depends only on the variant sequence. Changes separated by at least one nucleotide in
    every placement stay apart: HGVS describes ``g.[3C>G;7dup]``, not ``g.3_4delinsGGT``
    (SVD-WG010, which proposed merging nearby changes, was rejected). Finally every
    deletion and insertion is shifted 3' as far as it goes before the next change.
    """
    pieces = sorted(pieces, key=lambda piece: (piece.edit.start, piece.edit.end))
    merging = True
    while merging:
        merging = False
        for index in range(len(pieces) - 1):
            first, second = pieces[index].edit, pieces[index + 1].edit
            if _can_touch(seq, first, second) or _same_codon(first, second, codon_of):
                merged = _trim(seq, _combine(seq, [first, second]))
                members = pieces[index].members | pieces[index + 1].members
                if merged.start == merged.end and not merged.alt:
                    raise _AlleleError(members, "The variants cancel each other out")
                pieces[index:index + 2] = [
                    _Piece(part, members) for part in _decompose(seq, merged, codon_of)
                ]
                merging = True
                break

    for index in range(len(pieces) - 1, -1, -1):
        edit = pieces[index].edit
        if not _is_fixed(edit):
            limit = pieces[index + 1].edit.start if index + 1 < len(pieces) else None
            pieces[index] = _Piece(_roll_right(seq, edit, limit), pieces[index].members)
    return pieces


def _components(groups: Iterable[frozenset[int]]) -> list[frozenset[int]]:
    """Join the groups sharing a member (union-find)."""
    parent: dict[int, int] = {}

    def root(item: int) -> int:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    for group in groups:
        items = sorted(group)
        for item in items:
            parent.setdefault(item, item)
        for item in items[1:]:
            parent[root(item)] = root(items[0])

    components: dict[int, set[int]] = {}
    for item in parent:
        components.setdefault(root(item), set()).add(item)
    return [frozenset(component) for component in components.values()]


def _row_pieces(pieces: list[_Piece], rows: list[frozenset]) -> list[list[Edit]]:
    """The normalised changes of one frame, split into the rows of the report."""
    row_of = {key: index for index, row in enumerate(rows) for key in row}
    changes: list[list[Edit]] = [[] for _ in rows]
    for piece in sorted(pieces, key=lambda piece: (piece.edit.start, piece.edit.end)):
        changes[row_of[next(iter(piece.members))]].append(piece.edit)
    return changes


def _format_changes(seq, changes: list[list[Edit]], label, repeat_allowed) -> list[list[tuple[str, str]]]:
    """HGVS (type, edit) of every change of every row, in one frame."""
    flat = sorted(
        ((row, position, edit) for row, edits in enumerate(changes) for position, edit in enumerate(edits)),
        key=lambda item: (item[2].start, item[2].end),
    )
    formatted: list[list[tuple[str, str]]] = [[None] * len(edits) for edits in changes]
    lower = 0
    for row, position, edit in flat:
        normalised = _classify(seq, edit, lower, repeat_allowed)
        formatted[row][position] = (normalised.kind, _format(normalised, seq, label))
        lower = edit.end
    return formatted


def _allele_edits(changes: list[list[Edit]], formatted: list[list[tuple[str, str]]]) -> list[str]:
    flat = sorted(
        ((edit, formatted[row][position][1]) for row, edits in enumerate(changes) for position, edit in enumerate(edits)),
        key=lambda item: (item[0].start, item[0].end),
    )
    return [text for _, text in flat]


class _Describer:
    def __init__(self, context: ReferenceContext):
        self.context = context
        self.seq = context.oriented.seq
        self.codon_table = get_codon_table(context.transl_table) if context.is_coding else None

    def plus_codon(self, index: int) -> int | None:
        context = self.context
        return context.codon_number(context.oriented.from_genome(context.plus.to_genome(index)))

    def coding_repeat_allowed(self, unit: str, start: int, end: int) -> bool:
        # HGVS allows c. repeat descriptions in the coding sequence only for 3n units.
        first, last = self.context.cds_offset, self.context.cds_offset + self.context.cds_length
        return len(unit) % 3 == 0 or end <= first or start >= last

    def placement(self, edit: Edit, lower: int, upper: int, members) -> Edit:
        """An equivalent placement of the edit inside the CDS, used for the analysis.

        ``lower`` and ``upper`` are the neighbouring changes the edit must not pass.
        """
        context = self.context
        first, last = context.cds_offset, context.cds_offset + context.cds_length
        size = edit.end - edit.start
        outside = _AlleleError(
            members,
            f"The variant lies outside of the {context.symbol} coding sequence "
            f"(c.1-c.{context.cds_length}, m.{context.first_genome_position}-"
            f"m.{context.last_genome_position})",
        )

        if size and edit.alt:
            if first <= edit.start and edit.end <= last:
                return edit
            raise outside

        left = _roll_left(self.seq, edit, lower)
        right = _roll_right(self.seq, edit, upper)
        if size:
            lowest, highest = max(left.start, first), min(right.start, last - size)
        else:
            # Between two CDS bases or right after the last one, which can complete a
            # partial stop codon; an insertion in front of c.1 leaves the CDS untouched.
            lowest, highest = max(left.start, first + 1), min(right.start, last)
        if lowest > highest:
            raise outside

        chosen = min(max(edit.start, lowest), highest)
        if size:
            return Edit(chosen, chosen + size, "")
        point, inserted = left.start, left.alt
        while point < chosen:
            inserted = inserted[1:] + self.seq[point]
            point += 1
        return Edit(chosen, chosen, inserted)

    def protein(self, variants: Iterable["VariantDescription"]) -> ProteinConsequence | None:
        if not self.context.is_coding:
            return None
        return predict_consequence(
            self.context.cds_sequence,
            [CdsEdit(item.cds_start, item.cds_end, item.alt) for item in variants],
            self.codon_table,
            codon_start=self.context.codon_start,
            transl_except=self.context.transl_except,
            complete_start=self.context.cds_complete,
        )

    def translation_end(self, variants: Iterable["VariantDescription"]) -> int | None:
        return translation_end(
            self.context.cds_sequence,
            [CdsEdit(item.cds_start, item.cds_end, item.alt) for item in variants],
            self.codon_table,
            codon_start=self.context.codon_start,
            transl_except=self.context.transl_except,
        )

    def describe(self, parsed: list[_Parsed]) -> tuple[list[VariantDescription], list[str], list[str]]:
        """Rows of the report and the m. and c. edits of the whole allele."""
        context = self.context
        oriented, plus = context.oriented, context.plus

        given: dict[int, Edit] = {}
        for number, item in enumerate(parsed):
            edit = _trim(oriented.seq, item.edit)
            if edit.start == edit.end and not edit.alt:
                raise _AlleleError({number}, f"'{item.text}' does not change the sequence")
            given[number] = edit

        numbers = sorted(given)
        for position, first in enumerate(numbers):
            for second in numbers[position + 1:]:
                if given[first] == given[second]:
                    raise _AlleleError({first, second}, "The same variant is listed more than once")
                if _conflict(given[first], given[second]):
                    raise _AlleleError(
                        {first, second},
                        "The variants overlap, so it is not clear which change happens first",
                    )

        # Every part is identified by (input piece number, part number).
        coding_parts = {
            (number, index): part
            for number, edit in given.items()
            for index, part in enumerate(_decompose(oriented.seq, edit, context.codon_number))
        }
        genomic_parts = {
            key: _trim(plus.seq, _convert(edit, oriented, plus)) for key, edit in coding_parts.items()
        }

        try:
            coding_pieces = _normalise_allele(
                oriented.seq,
                [_Piece(edit, frozenset({key})) for key, edit in coding_parts.items()],
                context.codon_number,
            )
            genomic_pieces = _normalise_allele(
                plus.seq,
                [_Piece(edit, frozenset({key})) for key, edit in genomic_parts.items()],
                self.plus_codon,
            )
        except _AlleleError as exc:
            raise _AlleleError({key[0] for key in exc.members}, str(exc)) from exc

        # Changes merged in either frame are reported in one row; each frame describes
        # them with its own pieces, since the 3' rule follows the strand.
        rows = sorted(
            _components(piece.members for piece in coding_pieces + genomic_pieces),
            key=lambda component: min(coding_parts[key].start for key in component),
        )
        coding = _row_pieces(coding_pieces, rows)
        genomic = _row_pieces(genomic_pieces, rows)

        m_texts = _format_changes(plus.seq, genomic, context.m_label, lambda *_: True)
        c_texts = (
            _format_changes(oriented.seq, coding, context.c_label, self.coding_repeat_allowed)
            if context.c_available
            else None
        )

        described = []
        for position, component in enumerate(rows):
            lower = max(coding_parts[key].end for key in rows[position - 1]) if position else 0
            upper = (
                min(coding_parts[key].start for key in rows[position + 1])
                if position + 1 < len(rows)
                else len(oriented.seq)
            )
            members = {key[0] for key in component}
            combined = _trim(
                oriented.seq, _combine(oriented.seq, [coding_parts[key] for key in component])
            )
            placed = self.placement(combined, lower, upper, members)
            described.append(
                self._row(
                    parsed,
                    members,
                    placed,
                    m_texts[position],
                    c_texts[position] if c_texts is not None else None,
                )
            )

        allele_m = _allele_edits(genomic, m_texts)
        allele_c = _allele_edits(coding, c_texts) if c_texts is not None else []
        return described, allele_m, allele_c

    def _row(self, parsed, members, placed, m_pieces, c_pieces) -> VariantDescription:
        context = self.context
        inputs = tuple(sorted({parsed[number].index for number in members}))
        texts = tuple(dict.fromkeys(parsed[number].text for number in sorted(members)))
        kinds = list(dict.fromkeys(kind for kind, _ in m_pieces))
        m_edits = tuple(text for _, text in m_pieces)
        c_edits = tuple(text for _, text in c_pieces) if c_pieces is not None else None
        row = VariantDescription(
            inputs=inputs,
            texts=texts,
            kind=" + ".join(kinds),
            m_edits=m_edits,
            c_edits=c_edits,
            m=_allele(context.m_reference, list(m_edits)),
            c=_allele(context.c_reference, list(c_edits)) if c_edits else None,
            cds_start=placed.start - context.cds_offset,
            cds_end=placed.end - context.cds_offset,
            ref=self.seq[placed.start:placed.end],
            alt=placed.alt,
        )
        consequence = self.protein([row])
        if consequence is None:
            return row
        return replace(row, protein=consequence, p=consequence.description(context.protein_id))


def describe_variants(context: ReferenceContext, descriptions: Iterable[str]) -> DescriptionResult:
    """Parse, validate, merge and normalise HGVS descriptions of one gene.

    Every input may be a ``c.`` or ``m.`` description with or without its reference
    sequence (``NC_012920.1:m.``, ``NC_012920.1(MT-ATP6):c.``) or an allele of variants
    in cis (``m.[8993T>G;9000del]``). All inputs are one allele, so every change is
    described in the context of the others. Invalid inputs are reported in ``errors``;
    the remaining ones are described.
    """
    parser, describer = _Parser(context), _Describer(context)
    parsed: list[_Parsed] = []
    errors: list[DescriptionError] = []
    for index, text in enumerate(descriptions):
        try:
            parsed.extend(parser.parse(index, text))
        except HgvsError as exc:
            errors.append(DescriptionError((index,), (str(text).strip(),), str(exc)))

    failed = {index for error in errors for index in error.inputs}
    parsed = [item for item in parsed if item.index not in failed]
    while True:
        try:
            rows, allele_m, allele_c = describer.describe(parsed) if parsed else ([], [], [])
            break
        except _AlleleError as exc:
            # Drop every input involved and describe the remaining allele again.
            inputs = tuple(sorted({parsed[number].index for number in exc.members}))
            texts = tuple(dict.fromkeys(item.text for item in parsed if item.index in inputs))
            errors.append(DescriptionError(inputs, texts, str(exc)))
            parsed = [item for item in parsed if item.index not in inputs]

    rows.sort(key=lambda item: item.cds_start)
    errors.sort(key=lambda error: error.inputs)
    return DescriptionResult(
        variants=tuple(rows),
        errors=tuple(errors),
        allele_m=_allele(context.m_reference, allele_m),
        allele_c=_allele(context.c_reference, allele_c) if context.c_available else None,
        allele_p=_protein_allele(describer, rows),
    )


def _allele(reference: str, edits: list[str]) -> str | None:
    if not edits:
        return None
    if len(edits) == 1:
        return f"{reference}{edits[0]}"
    return f"{reference}[{';'.join(edits)}]"


def _touching(first: ProteinConsequence, second: ProteinConsequence) -> bool:
    """Changes of consecutive amino acids, which HGVS describes as one deletion-insertion."""
    if first.span is None or second.span is None:
        return False
    return first.span[1] >= second.span[0]


def _protein_allele(describer: _Describer, variants: list[VariantDescription]) -> str | None:
    """The predicted protein of all variants together, ``p.[(Leu156Arg;Asn4del)]``.

    Changes of consecutive amino acids are one deletion-insertion. Nothing after a
    premature stop codon is described, and a frameshift is predicted together with the
    changes its shifted frame still reads (a later change can restore the frame).
    """
    if not variants or variants[0].protein is None:
        return None

    accession = describer.context.protein_id
    units = [[variant] for variant in variants]
    while True:
        consequences = [
            unit[0].protein if len(unit) == 1 else describer.protein(unit) for unit in units
        ]
        # A changed start codon makes the whole protein unpredictable.
        if any(consequence.kind == START for consequence in consequences):
            return ProteinConsequence("Met1?", START).description(accession)

        changed = False
        for index, consequence in enumerate(consequences):
            if consequence.kind == FRAMESHIFT:
                end = describer.translation_end(units[index])
                following = index + 1
                while following < len(units) and (end is None or units[following][0].cds_start < end):
                    following += 1
                if following > index + 1:
                    units[index:following] = [[item for unit in units[index:following] for item in unit]]
                    changed = True
                    break
            if consequence.terminal:
                if index + 1 < len(units):
                    del units[index + 1:]
                    changed = True
                break
            if index + 1 < len(units) and _touching(consequence, consequences[index + 1]):
                units[index:index + 2] = [units[index] + units[index + 1]]
                changed = True
                break
        if not changed:
            break

    edits = [consequence.edit for consequence in consequences if consequence.edit != "="]
    if not edits:
        edits = ["="]
    if len(edits) == 1:
        body = f"p.({edits[0]})"
    else:
        body = f"p.[({';'.join(edits)})]"
    return f"{accession}:{body}" if accession else body


__all__ = [
    "HGVS_VERSION",
    "DescriptionError",
    "DescriptionResult",
    "HgvsError",
    "ReferenceContext",
    "VariantDescription",
    "describe_variants",
    "reverse_complement",
    "three_letter",
]
