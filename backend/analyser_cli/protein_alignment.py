"""Side by side comparison of the original and the mutated protein translation.

The two translations rarely have the same length - an indel shifts everything that
follows it - so they are aligned before being compared. That keeps the residues that
did not change on top of each other and shows an in-frame deletion as a gap instead of
marking the whole tail of the protein as different.

Pairwise alignment only tells the truth while both translations stay homologous, which
is the case for an in-frame indel. A frameshift reads the rest of the gene in another
frame, so everything past the indel is a different peptide rather than a shifted one,
and a premature stop codon simply cuts the protein short. Forcing an aligner onto
either case makes it chase the ~6% of residue pairs that match by chance and scatter
gaps through the tail, which reads as a series of small edits instead of the single
catastrophic one that actually happened. The residues before the frameshift are still
homologous and are aligned normally; only the tail after it is anchored and padded.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio import Align

GAP = "-"
STOP = "*"
DEFAULT_LINE_LENGTH = 60
DEFAULT_GROUP_SIZE = 5


@dataclass(frozen=True)
class AlignmentLine:
    """One wrapped line of the alignment, ready to be rendered."""

    original: str
    mutated: str
    differences: tuple[bool, ...]
    original_start: int | None
    original_end: int | None
    mutated_start: int | None
    mutated_end: int | None
    divergence: int | None = None
    """Offset inside this line where the two translations stop being homologous."""


@dataclass(frozen=True)
class ProteinAlignment:
    """The original and mutated translation padded with gaps to the same length."""

    original: str
    mutated: str
    divergence_start: int | None = None
    """1-based column from which the mutated translation is no longer homologous."""
    divergence_reason: str | None = None
    """``'frameshift'`` or ``'truncation'`` when the tail is not comparable."""

    @property
    def has_mutated(self) -> bool:
        return bool(self.mutated.replace(GAP, ""))

    @property
    def is_homologous(self) -> bool:
        """Whether every column can be read as a substitution, an insertion or a deletion."""
        return self.divergence_start is None

    @property
    def differences(self) -> tuple[bool, ...]:
        # Past a frameshift or a premature stop codon the two translations no longer
        # describe the same residues, so the columns that happen to carry the same letter
        # are a coincidence rather than a conserved residue and are reported as changed.
        divergence = (
            len(self.original) if self.divergence_start is None else self.divergence_start - 1
        )

        return tuple(
            left != right or index >= divergence
            for index, (left, right) in enumerate(zip(self.original, self.mutated))
        )

    @property
    def identical_residues(self) -> int:
        return sum(1 for different in self.differences if not different)

    @property
    def changed_residues(self) -> int:
        return sum(1 for different in self.differences if different)

    @property
    def original_length(self) -> int:
        return len(self.original.replace(GAP, ""))

    @property
    def mutated_length(self) -> int:
        return len(self.mutated.replace(GAP, ""))

    def lines(
        self,
        line_length: int = DEFAULT_LINE_LENGTH,
        group_size: int = DEFAULT_GROUP_SIZE,
    ) -> list[AlignmentLine]:
        """Split the alignment into fixed width lines so it never has to scroll."""
        if line_length < 1:
            raise ValueError("line_length must be at least 1")

        # Keep whole groups on a line, otherwise the grouping spaces would not line up.
        if group_size > 0:
            line_length = max(line_length - line_length % group_size, group_size)

        differences = self.differences
        original_position = 0
        mutated_position = 0
        lines = []

        for start in range(0, len(self.original), line_length):
            end = start + line_length
            original_chunk = self.original[start:end]
            mutated_chunk = self.mutated[start:end]

            original_first, original_position = _advance(original_chunk, original_position)
            mutated_first, mutated_position = _advance(mutated_chunk, mutated_position)

            lines.append(
                AlignmentLine(
                    original=original_chunk,
                    mutated=mutated_chunk,
                    differences=differences[start:end],
                    original_start=original_first,
                    original_end=original_position if original_first else None,
                    mutated_start=mutated_first,
                    mutated_end=mutated_position if mutated_first else None,
                    divergence=_divergence_offset(
                        self.divergence_start, start, len(original_chunk)
                    ),
                )
            )

        return lines


def align_translations(
    original: object,
    mutated: object,
    *,
    frameshift_residue: int | None = None,
) -> ProteinAlignment:
    """Align two protein translations, padding both with gaps to the same length.

    ``frameshift_residue`` is the 1-based residue of the original translation whose codon
    is the first one read in another frame. The residues before it are still homologous
    and are aligned normally, while the tail after it is anchored and padded, because the
    two tails are different peptides rather than shifted copies of each other. Pass ``1``
    when a frameshift is known but its position is not. A translation cut short by a
    premature stop codon is detected from the sequences themselves.
    """
    original = str(original or "")
    mutated = str(mutated or "")

    if not mutated:
        return ProteinAlignment(original, GAP * len(original))

    if not original:
        return ProteinAlignment(GAP * len(mutated), mutated)

    if original == mutated:
        return ProteinAlignment(original, mutated)

    if frameshift_residue is not None and frameshift_residue <= len(original):
        return _split_at_frameshift(original, mutated, max(frameshift_residue, 1))

    if _is_truncation(original, mutated):
        return _anchor_tails(original, mutated, "", "", "truncation")

    aligned_original, aligned_mutated = _normalise_gaps(*_align(original, mutated))

    return ProteinAlignment(aligned_original, aligned_mutated)


def group_residues(residues: str, group_size: int = DEFAULT_GROUP_SIZE) -> str:
    """Insert a space every ``group_size`` residues to make counting them easier."""
    if group_size < 1:
        return residues

    return " ".join(
        residues[start:start + group_size] for start in range(0, len(residues), group_size)
    )


def _advance(chunk: str, position: int) -> tuple[int | None, int]:
    """Return the 1-based position of the first residue of the chunk and the last one."""
    first = None

    for residue in chunk:
        if residue == GAP:
            continue
        position += 1
        if first is None:
            first = position

    return first, position


def _divergence_offset(divergence_start: int | None, start: int, length: int) -> int | None:
    """Return where the divergence falls inside a wrapped line, if it falls in it at all."""
    if divergence_start is None:
        return None

    offset = divergence_start - 1 - start

    return offset if 0 <= offset < length else None


def _is_truncation(original: str, mutated: str) -> bool:
    """Whether the mutated translation is the original one cut short by a stop codon."""
    if len(mutated) >= len(original):
        return False

    return original.startswith(mutated[:-1] if mutated.endswith(STOP) else mutated)


def _split_at_frameshift(original: str, mutated: str, residue: int) -> ProteinAlignment:
    """Align the homologous residues before the frameshift and anchor the tails after it."""
    original_head = original[:residue - 1]
    mutated_head = mutated[:_matching_prefix_length(original_head, mutated)]
    aligned_head = _normalise_gaps(*_align(original_head, mutated_head))

    return _anchor_tails(
        original[len(original_head):], mutated[len(mutated_head):], *aligned_head, "frameshift"
    )


def _anchor_tails(
    original_tail: str,
    mutated_tail: str,
    aligned_original_head: str,
    aligned_mutated_head: str,
    reason: str,
) -> ProteinAlignment:
    """Put the two tails side by side and pad the shorter one at the end."""
    if reason == "truncation":
        # A premature stop codon leaves the start of the protein untouched, so whatever the
        # two translations still share is the homologous head.
        shared = _common_prefix_length(original_tail, mutated_tail)
        aligned_original_head = aligned_mutated_head = original_tail[:shared]
        original_tail, mutated_tail = original_tail[shared:], mutated_tail[shared:]

    width = max(len(original_tail), len(mutated_tail))

    return ProteinAlignment(
        aligned_original_head + original_tail.ljust(width, GAP),
        aligned_mutated_head + mutated_tail.ljust(width, GAP),
        divergence_start=len(aligned_original_head) + 1,
        divergence_reason=reason,
    )


def _matching_prefix_length(original_head: str, mutated: str) -> int:
    """Return how much of ``mutated`` the still homologous head of the original covers.

    An in-frame indel before the frameshift moves the boundary, so the two heads do not
    have to be the same length and the aligner has to say where the mutated one ends.
    """
    if not original_head:
        return 0

    aligner = _aligner()
    # The mutated tail is a different peptide, so leaving it out has to be free.
    aligner.open_right_insertion_score = 0.0
    aligner.extend_right_insertion_score = 0.0

    return int(aligner.align(mutated, original_head)[0].aligned[0][-1][-1])


def _common_prefix_length(original: str, mutated: str) -> int:
    for index, (left, right) in enumerate(zip(original, mutated)):
        if left != right:
            return index

    return min(len(original), len(mutated))


def _align(original: str, mutated: str) -> tuple[str, str]:
    if original == mutated:
        return original, mutated

    alignment = _aligner().align(original, mutated)[0]

    return str(alignment[0]), str(alignment[1])


def _normalise_gaps(original: str, mutated: str) -> tuple[str, str]:
    """Push every gap to the rightmost position that describes the same indel.

    A residue deleted from a repeat can be placed anywhere inside it, and the aligner
    picks the leftmost spot. HGVS names the rightmost one, so the gaps are moved to match
    the notation the rest of the report uses.
    """
    original = _shift_gaps_right(original, mutated)

    return original, _shift_gaps_right(mutated, original)


def _shift_gaps_right(sequence: str, reference: str) -> str:
    """Slide the gap runs of ``sequence`` right while that pairs up the same residues."""
    residues = list(sequence)
    start = 0

    while start < len(residues):
        if residues[start] != GAP:
            start += 1
            continue

        end = start
        while end < len(residues) and residues[end] == GAP:
            end += 1

        # Swapping the gap run with the residue that follows it leaves the sequence itself
        # untouched, so it is only worth doing while that residue keeps its partner.
        while end < len(residues) and reference[start] == reference[end] != GAP:
            residues[start], residues[end] = residues[end], residues[start]
            start += 1
            end += 1

        start = end

    return "".join(residues)


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    # A plain identity scoring keeps every alphabet working, including the '*' of a stop
    # codon and the 'U'/'O' of an annotated /transl_except. Gaps are expensive so that a
    # single changed residue is shown as a substitution rather than as a deletion plus an
    # insertion, and end gaps are penalised like any other so a residue lost in the middle
    # is not pushed to the end of the sequence.
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    return aligner
