"""Side by side comparison of the original and the mutated protein translation.

The two translations rarely have the same length - an indel shifts everything that
follows it - so they are aligned before being compared. That keeps the residues that
did not change on top of each other and shows an in-frame deletion as a gap instead of
marking the whole tail of the protein as different.
"""

from __future__ import annotations

from dataclasses import dataclass

from Bio import Align

GAP = "-"
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


@dataclass(frozen=True)
class ProteinAlignment:
    """The original and mutated translation padded with gaps to the same length."""

    original: str
    mutated: str

    @property
    def has_mutated(self) -> bool:
        return bool(self.mutated.replace(GAP, ""))

    @property
    def differences(self) -> tuple[bool, ...]:
        return tuple(left != right for left, right in zip(self.original, self.mutated))

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
                )
            )

        return lines


def align_translations(original: object, mutated: object) -> ProteinAlignment:
    """Align two protein translations, padding both with gaps to the same length."""
    original = str(original or "")
    mutated = str(mutated or "")

    if not mutated:
        return ProteinAlignment(original, GAP * len(original))

    if not original:
        return ProteinAlignment(GAP * len(mutated), mutated)

    if original == mutated:
        return ProteinAlignment(original, mutated)

    alignment = _aligner().align(original, mutated)[0]

    return ProteinAlignment(str(alignment[0]), str(alignment[1]))


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
