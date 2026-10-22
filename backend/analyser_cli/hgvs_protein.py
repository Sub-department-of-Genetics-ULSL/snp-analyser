"""Predicted protein consequences of DNA variants, described as HGVS ``p.`` edits.

HGVS derives a protein description by comparing the variant protein with the reference
protein, never from the DNA change itself, so both are translated here:

* the reference protein reproduces the GenBank ``/translation`` (an alternative start
  codon initiates with methionine, ``/transl_except`` codons are honoured) and ends with
  its stop codon;
* the variant protein is read from the variant CDS followed by a poly(A) tail, which is
  how a processed mitochondrial mRNA ends. A partial 3' stop codon is completed by that
  tail exactly like GenBank annotates it, and a frame that runs into the tail never
  meets a stop codon again (``Ter?``).

The rules applied are the ones of HGVS version 21.1: ``Met1?`` for a changed initiation
codon, ``=`` for a silent change, ``Ter`` substitutions for an immediate stop codon,
``ext`` for a lost stop codon, ``fs`` counted from the first changed amino acid, and
in-frame deletions, insertions and duplications shifted to the most C-terminal position.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "U": "Sec",
    "O": "Pyl", "X": "Xaa", "*": "Ter",
}

START = "start"
SYNONYMOUS = "synonymous"
MISSENSE = "missense"
NONSENSE = "nonsense"
EXTENSION = "extension"
FRAMESHIFT = "frameshift"
INFRAME_DELETION = "inframe_deletion"
INFRAME_INSERTION = "inframe_insertion"
INFRAME_DUPLICATION = "inframe_duplication"
INFRAME_DELINS = "inframe_delins"
UNKNOWN = "unknown"

CONSEQUENCE_LABELS = {
    START: "Start codon changed (consequence unknown)",
    SYNONYMOUS: "Synonymous (silent)",
    MISSENSE: "Missense",
    NONSENSE: "Nonsense (premature stop codon)",
    EXTENSION: "Stop codon lost (extension)",
    FRAMESHIFT: "Frameshift",
    INFRAME_DELETION: "In-frame deletion",
    INFRAME_INSERTION: "In-frame insertion",
    INFRAME_DUPLICATION: "In-frame duplication",
    INFRAME_DELINS: "In-frame deletion-insertion",
    UNKNOWN: "Unknown",
}

# A frame running into the poly(A) tail reads AAA codons, which no genetic code stops on.
POLY_A_CODON = "AAA"


def three_letter(residues: str) -> str:
    """Spell residues with the three-letter amino acid code HGVS recommends."""
    return "".join(THREE_LETTER.get(residue, "Xaa") for residue in residues)


@dataclass(frozen=True)
class CdsEdit:
    """Replace ``cds[start:end]`` with ``alt`` (0-based, coding strand)."""

    start: int
    end: int
    alt: str

    @property
    def length_change(self) -> int:
        return len(self.alt) - (self.end - self.start)


@dataclass(frozen=True)
class ProteinConsequence:
    """One predicted protein change, e.g. ``Leu156Arg`` with the kind ``missense``.

    ``span`` holds the boundaries (0-based, between residues) of the reference residues
    the change touches; ``None`` when no residue changes.
    """

    edit: str
    kind: str
    first_residue: int | None = None
    span: tuple[int, int] | None = None

    @property
    def label(self) -> str:
        return CONSEQUENCE_LABELS[self.kind]

    @property
    def terminal(self) -> bool:
        """Whether the protein ends (or leaves the reference frame) with this change."""
        if self.kind in {NONSENSE, EXTENSION, FRAMESHIFT}:
            return True
        return self.kind in {INFRAME_INSERTION, INFRAME_DELINS} and "Ter" in self.edit

    def description(self, accession: str | None = None) -> str:
        """HGVS description of a predicted consequence, which goes in parentheses."""
        body = f"p.({self.edit})"
        return f"{accession}:{body}" if accession else body


def apply_edits(sequence: str, edits: Iterable[CdsEdit]) -> str:
    """Apply non-overlapping edits given in reference coordinates, 5' to 3'."""
    result = sequence
    # Right to left, so every edit still finds its reference coordinates; an insertion
    # sharing its position with a deletion goes in front of it, and of two insertions at
    # one position the later one is applied first so the earlier one ends up in front.
    ordered = sorted(
        enumerate(edits), key=lambda item: (item[1].start, item[1].end, item[0]), reverse=True
    )
    for _, edit in ordered:
        result = result[: edit.start] + edit.alt + result[edit.end :]
    return result


def translate_protein(
    sequence: str,
    codon_table,
    *,
    codon_start: int = 1,
    exceptions: dict[int, str] | None = None,
) -> tuple[str, bool]:
    """Translate up to and including the first stop codon.

    A partial last codon is completed with the poly(A) tail. Returns the protein and
    whether it reached a stop codon.
    """
    coding = sequence[max(codon_start - 1, 0):]
    if len(coding) % 3:
        coding += "A" * (3 - len(coding) % 3)

    exceptions = exceptions or {}
    residues = []
    for index in range(0, len(coding), 3):
        codon = coding[index:index + 3]
        codon_number = index // 3
        if codon_number in exceptions:
            residue = exceptions[codon_number]
        elif codon in codon_table.stop_codons:
            residue = "*"
        else:
            residue = codon_table.forward_table.get(codon, "X")
            # Translation initiates with methionine even on an alternative start codon.
            if codon_number == 0 and codon in codon_table.start_codons:
                residue = "M"

        residues.append(residue)
        if residue == "*":
            return "".join(residues), True

    return "".join(residues), False


def predict_consequence(
    cds: str,
    edits: Sequence[CdsEdit],
    codon_table,
    *,
    codon_start: int = 1,
    transl_except: Iterable = (),
    complete_start: bool = True,
) -> ProteinConsequence:
    """Describe the protein change caused by ``edits`` applied together to ``cds``.

    ``edits`` holds one variant, or a frameshift with the variants downstream of it in
    the shifted frame. ``transl_except`` takes
    :class:`~analyser_cli.cds_annotation.TranslationException` objects.
    """
    if not edits:
        raise ValueError("At least one edit is required")

    offset = max(codon_start - 1, 0)
    edits = sorted(edits, key=lambda item: (item.start, item.end))
    reference_exceptions = {
        (exception.start - offset) // 3: exception.amino_acid
        for exception in transl_except
        if exception.start >= offset
    }
    reference, reference_terminated = translate_protein(
        cds, codon_table, codon_start=codon_start, exceptions=reference_exceptions
    )
    variant_cds = apply_edits(cds, edits)

    if (
        complete_start
        and reference.startswith("M")
        and variant_cds[offset:offset + 3] != cds[offset:offset + 3]
    ):
        return ProteinConsequence("Met1?", START, 1, (0, 1))

    variant, variant_terminated = translate_protein(
        variant_cds,
        codon_table,
        codon_start=codon_start,
        exceptions=_unaffected_exceptions(reference_exceptions, cds, edits, offset),
    )

    if variant == reference and variant_terminated == reference_terminated:
        return _synonymous(edits, offset, reference)

    padding = codon_table.forward_table.get(POLY_A_CODON, "X")
    view = variant if variant_terminated else variant + padding * (len(reference) + 1)
    first = _first_difference(reference, view)
    if first >= len(reference):
        return ProteinConsequence("?", UNKNOWN)

    reference_residue, variant_residue = reference[first], view[first]
    if variant_residue == "*":
        return ProteinConsequence(
            f"{three_letter(reference_residue)}{first + 1}Ter", NONSENSE, first + 1, (first, first + 1)
        )

    if reference_residue == "*":
        length = str(len(variant) - 1 - first) if variant_terminated else "?"
        return ProteinConsequence(
            f"Ter{first + 1}{three_letter(variant_residue)}extTer{length}",
            EXTENSION,
            first + 1,
            (first, first + 1),
        )

    stop_index = len(reference) - 1 if reference_terminated else None
    # HGVS describes deletions starting 5' of and including the stop codon as frameshifts.
    touches_stop = stop_index is not None and any(
        _overlaps_codon(edit, offset + stop_index * 3) for edit in edits
    )
    if sum(edit.length_change for edit in edits) % 3 or touches_stop:
        length = str(len(variant) - first) if variant_terminated else "?"
        return ProteinConsequence(
            f"{three_letter(reference_residue)}{first + 1}"
            f"{three_letter(variant_residue)}fsTer{length}",
            FRAMESHIFT,
            first + 1,
            (first, first + 1),
        )

    # A frameshift compensated by a later change only replaces the stretch between them.
    start, end = edits[0].start, max(edit.end for edit in edits)
    combined = CdsEdit(
        start,
        end,
        apply_edits(
            cds[start:end],
            [CdsEdit(edit.start - start, edit.end - start, edit.alt) for edit in edits],
        ),
    )
    return _in_frame(combined, offset, reference, view)


def translation_end(
    cds: str,
    edits: Sequence[CdsEdit],
    codon_table,
    *,
    codon_start: int = 1,
    transl_except: Iterable = (),
) -> int | None:
    """Reference CDS index right after the stop codon ending the variant protein.

    Variants from there on cannot change the protein. ``None`` when the frame runs into
    the poly(A) tail without meeting a stop codon.
    """
    offset = max(codon_start - 1, 0)
    edits = sorted(edits, key=lambda item: (item.start, item.end))
    reference_exceptions = {
        (exception.start - offset) // 3: exception.amino_acid
        for exception in transl_except
        if exception.start >= offset
    }
    protein, terminated = translate_protein(
        apply_edits(cds, edits),
        codon_table,
        codon_start=codon_start,
        exceptions=_unaffected_exceptions(reference_exceptions, cds, edits, offset),
    )
    if not terminated:
        return None

    variant_position = offset + 3 * len(protein)
    shift = 0
    for edit in edits:
        start = edit.start + shift
        if variant_position <= start:
            break
        if variant_position <= start + len(edit.alt):
            # The stop codon ends inside the inserted bases.
            return edit.end
        shift += edit.length_change
    return variant_position - shift


def _unaffected_exceptions(
    exceptions: dict[int, str], cds: str, edits: Sequence[CdsEdit], offset: int
) -> dict[int, str]:
    """Keep the internal ``/transl_except`` codons no edit touches or moves.

    The terminal one completes a partial stop codon; the poly(A) tail does that for the
    variant, so it is read from the variant bases instead.
    """
    complete_codons = (len(cds) - offset) // 3
    kept = {}
    for codon_number, residue in exceptions.items():
        if codon_number >= complete_codons:
            continue
        codon_start = offset + codon_number * 3
        before_all = all(codon_start + 3 <= edit.start for edit in edits)
        untouched = all(not _overlaps_codon(edit, codon_start) for edit in edits)
        if before_all or (untouched and sum(edit.length_change for edit in edits) == 0):
            kept[codon_number] = residue
    return kept


def _overlaps_codon(edit: CdsEdit, codon_start: int) -> bool:
    codon_end = codon_start + 3
    if edit.end > edit.start:
        return edit.start < codon_end and codon_start < edit.end
    return codon_start < edit.start < codon_end


def _first_difference(reference: str, variant: str) -> int:
    for index, (left, right) in enumerate(zip(reference, variant)):
        if left != right:
            return index
    return min(len(reference), len(variant))


def _synonymous(edits: Sequence[CdsEdit], offset: int, reference: str) -> ProteinConsequence:
    codons = set()
    for edit in edits:
        start, end = edit.start - offset, edit.end - offset
        if end > start:
            codons.update(range(start // 3, (end - 1) // 3 + 1))
        elif start % 3:
            codons.add(start // 3)

    codons = {codon for codon in codons if 0 <= codon < len(reference)}
    if len(codons) == 1:
        codon = codons.pop()
        return ProteinConsequence(
            f"{three_letter(reference[codon])}{codon + 1}=", SYNONYMOUS, codon + 1
        )
    return ProteinConsequence("=", SYNONYMOUS)


def _in_frame(edit: CdsEdit, offset: int, reference: str, variant: str) -> ProteinConsequence:
    """Describe an in-frame change from the codons the DNA edit touches.

    Working on those codons keeps a stop codon inside an inserted sequence an insertion
    (``Met3_His4insGlyTer``) instead of a deletion-insertion of the whole C-terminus.
    """
    start, end = edit.start - offset, edit.end - offset
    window_start = (start // 3) * 3
    if end > start:
        window_end = -(-end // 3) * 3
    else:
        window_end = window_start if start % 3 == 0 else window_start + 3

    position = window_start // 3
    removed = reference[position:window_end // 3]
    added = variant[position:(window_end + edit.length_change) // 3]
    if "*" in added:
        added = added[: added.index("*") + 1]

    while removed and added and removed[0] == added[0]:
        removed, added, position = removed[1:], added[1:], position + 1
    if "*" not in added:
        while removed and added and removed[-1] == added[-1]:
            removed, added = removed[:-1], added[:-1]

    if not removed and not added:
        return _synonymous([edit], offset, reference)

    if added.startswith("*"):
        return ProteinConsequence(
            f"{three_letter(reference[position])}{position + 1}Ter",
            NONSENSE,
            position + 1,
            (position, position + 1),
        )

    # The protein keeps its stop codon at the end; nothing is shifted across it.
    last_sense = len(reference) - 1 if reference.endswith("*") else len(reference)

    if not removed:
        if "*" not in added:
            while position < last_sense and reference[position] == added[0]:
                added = added[1:] + reference[position]
                position += 1
            size = len(added)
            if position >= size and reference[position - size:position] == added:
                first = position - size + 1
                return ProteinConsequence(
                    f"{_residue_range(reference, first, position)}dup",
                    INFRAME_DUPLICATION,
                    first,
                    (position, position),
                )
        return ProteinConsequence(
            f"{three_letter(reference[position - 1])}{position}_"
            f"{three_letter(reference[position])}{position + 1}ins{three_letter(added)}",
            INFRAME_INSERTION,
            position + 1,
            (position, position),
        )

    if not added:
        size = len(removed)
        while position + size < last_sense and reference[position] == reference[position + size]:
            position += 1
        return ProteinConsequence(
            f"{_residue_range(reference, position + 1, position + size)}del",
            INFRAME_DELETION,
            position + 1,
            (position, position + size),
        )

    if len(removed) == 1 and len(added) == 1:
        return ProteinConsequence(
            f"{three_letter(removed)}{position + 1}{three_letter(added)}",
            MISSENSE,
            position + 1,
            (position, position + 1),
        )

    return ProteinConsequence(
        f"{_residue_range(reference, position + 1, position + len(removed))}"
        f"delins{three_letter(added)}",
        INFRAME_DELINS,
        position + 1,
        (position, position + len(removed)),
    )


def _residue_range(reference: str, first: int, last: int) -> str:
    """``Lys23`` or ``Lys23_Val25`` for 1-based residue numbers."""
    start = f"{three_letter(reference[first - 1])}{first}"
    if first == last:
        return start
    return f"{start}_{three_letter(reference[last - 1])}{last}"
