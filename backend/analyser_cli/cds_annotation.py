"""Translation related qualifiers read from a GenBank CDS feature.

GenBank annotates every coding sequence with the information needed to reproduce its
``/translation`` exactly:

* ``/transl_table`` - the NCBI genetic code id. It is omitted when the standard code
  (id 1) applies, therefore the default here is 1 and never a "guessed" organism group.
* ``/codon_start`` - 1-based offset of the first complete codon inside the feature.
* ``/transl_except`` - single codons whose amino acid differs from the genetic code,
  used in mitochondrial genomes for 3' partial stop codons completed by polyadenylation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from Bio.Data import CodonTable

DEFAULT_TRANSL_TABLE = 1
DEFAULT_CODON_START = 1

# /transl_except spells the amino acid with a three letter code.
AMINO_ACID_CODES = {
    "TERM": "*",
    "SEC": "U",
    "PYL": "O",
    "OTHER": "X",
}

_TRANSL_EXCEPT_PATTERN = re.compile(
    r"\(\s*pos\s*:\s*(?P<pos>.+?)\s*,\s*aa\s*:\s*(?P<aa>[A-Za-z]+)", re.IGNORECASE
)
_POSITION_PATTERN = re.compile(r"(?P<start>\d+)(?:\.\.(?P<end>\d+))?")


@dataclass(frozen=True)
class TranslationException:
    """A ``/transl_except`` entry mapped onto gene relative coordinates.

    Attributes:
        start: 0-based offset of the first affected base inside the extracted gene.
        end: exclusive 0-based offset of the last affected base.
        amino_acid: single letter amino acid the codon translates to.
    """

    start: int
    end: int
    amino_acid: str


@dataclass(frozen=True)
class CdsAnnotation:
    """Everything needed to translate a gene the way GenBank does."""

    transl_table: int = DEFAULT_TRANSL_TABLE
    codon_start: int = DEFAULT_CODON_START
    transl_except: tuple[TranslationException, ...] = field(default_factory=tuple)
    translation: str | None = None


def get_codon_table(transl_table: int = DEFAULT_TRANSL_TABLE):
    """Return the NCBI genetic code with the given ``/transl_table`` id.

    Args:
        transl_table: NCBI genetic code id as annotated on a GenBank CDS feature.

    Raises:
        ValueError: If the id is not a known NCBI genetic code.
    """
    try:
        return CodonTable.unambiguous_dna_by_id[int(transl_table)]
    except (KeyError, TypeError, ValueError):
        raise ValueError(
            f"Unknown NCBI genetic code id {transl_table!r}, expected one of "
            f"{sorted(CodonTable.unambiguous_dna_by_id)}"
        )


def genetic_code_name(transl_table: int = DEFAULT_TRANSL_TABLE) -> str:
    """Human-readable name of a genetic code, e.g. 'Vertebrate Mitochondrial'."""
    return get_codon_table(transl_table).names[0]


def read_cds_annotation(feature) -> CdsAnnotation:
    """Read the translation qualifiers of a GenBank feature.

    Features without the qualifiers (a bare ``gene`` feature for example) fall back to
    the standard genetic code, which is exactly what GenBank means by their absence.
    """
    qualifiers = getattr(feature, "qualifiers", {}) or {}

    return CdsAnnotation(
        transl_table=_read_int(qualifiers, "transl_table", DEFAULT_TRANSL_TABLE),
        codon_start=_read_int(qualifiers, "codon_start", DEFAULT_CODON_START),
        transl_except=read_transl_except(feature),
        translation=(qualifiers.get("translation") or [None])[0],
    )


def read_transl_except(feature) -> tuple[TranslationException, ...]:
    """Translate ``/transl_except`` genome positions into gene relative offsets."""
    qualifiers = getattr(feature, "qualifiers", {}) or {}
    raw_entries = qualifiers.get("transl_except") or []
    if not raw_entries:
        return ()

    # Iterating a location yields genome positions in transcription order, so this also
    # covers genes annotated on the complement strand and joined locations.
    genome_to_gene = {position: offset for offset, position in enumerate(feature.location)}

    exceptions = []
    for raw_entry in raw_entries:
        parsed = _TRANSL_EXCEPT_PATTERN.search(str(raw_entry))
        if parsed is None:
            continue

        amino_acid = AMINO_ACID_CODES.get(parsed.group("aa").upper())
        position = _POSITION_PATTERN.search(parsed.group("pos"))
        if amino_acid is None or position is None:
            continue

        first = int(position.group("start")) - 1
        last = int(position.group("end") or position.group("start")) - 1
        offsets = [genome_to_gene[base] for base in range(first, last + 1) if base in genome_to_gene]
        if not offsets:
            continue

        exceptions.append(TranslationException(min(offsets), max(offsets) + 1, amino_acid))

    return tuple(sorted(exceptions, key=lambda exception: exception.start))


def _read_int(qualifiers: dict, name: str, default: int) -> int:
    try:
        return int(qualifiers[name][0])
    except (KeyError, IndexError, TypeError, ValueError):
        return default
