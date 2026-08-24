import Bio
from Bio.Seq import Seq, MutableSeq

from .cds_annotation import (
    DEFAULT_CODON_START,
    DEFAULT_TRANSL_TABLE,
    TranslationException,
    get_codon_table,
)

VALID_NUCLEOTIDES = frozenset("ACGT")
SUPPORTED_MUTATION_TYPES = frozenset({"sub", "del", "ins"})


def _as_translation_exceptions(entries) -> tuple[TranslationException, ...]:
    """Accept both ``TranslationException`` objects and their dict form."""
    return tuple(
        entry
        if isinstance(entry, TranslationException)
        else TranslationException(
            start=int(entry["start"]),
            end=int(entry["end"]),
            amino_acid=str(entry["aminoAcid"]),
        )
        for entry in (entries or ())
    )


class Analyser:
    """Analyser class for processing DNA sequences and their amino acid translations.

    The genetic code is never guessed from the organism: it is taken from the
    ``/transl_table`` qualifier of the GenBank CDS feature the sequence comes from,
    together with ``/codon_start`` and ``/transl_except``. That makes the translation
    reproduce the ``/translation`` qualifier of the record exactly.

    Attributes:
        dna_type (str): Descriptive label of the sequence ('mitochondrial' or 'nuclear').
            It documents the report only, the genetic code comes from ``transl_table``.
        transl_table (int): NCBI genetic code id used for every translation.
        codon_table (Bio.Data.CodonTable.CodonTable): Resolved codon table.
        codon_start (int): 1-based offset of the first complete codon in the sequence.
        transl_except (tuple[TranslationException, ...]): Annotated codon exceptions.
        input_sequence (Bio.Seq.Seq): Original DNA sequence.
        amino_acid_translation (Bio.Seq.Seq): Amino acid translation of original sequence.
        mutated_sequence (Bio.Seq.Seq or None): DNA sequence after applying mutations.
        mutated_amino_acid_translation (Bio.Seq.Seq or None): Amino acid translation of mutated sequence.

    Example:
        >>> analyser = Analyser("ATGCGATCG", transl_table=2)
        >>> mutations = [{"type": "sub", "position": 3, "alt": "T"}]
        >>> analyser.apply_mutations(mutations)
        >>> print(analyser.amino_acid_translation)
        >>> print(analyser.mutated_amino_acid_translation)
    """

    def __init__(
        self,
        input_sequence: str,
        dna_type: str = "mitochondrial",
        transl_table: int = DEFAULT_TRANSL_TABLE,
        codon_start: int = DEFAULT_CODON_START,
        transl_except: tuple[TranslationException, ...] = (),
    ):
        self.dna_type = dna_type
        self.transl_table = int(transl_table)
        self.codon_table = get_codon_table(self.transl_table)
        self.codon_start = int(codon_start)
        self.transl_except = _as_translation_exceptions(transl_except)
        self.input_sequence = Seq(input_sequence)
        self.amino_acid_translation = self._translate(self.input_sequence, self.transl_except)

        self.mutated_sequence = None
        self.mutated_amino_acid_translation = None

    @property
    def codon_table_name(self) -> str:
        """Human-readable name of the codon table in use."""
        return self.codon_table.names[0]

    @property
    def genetic_code_description(self) -> str:
        """Codon table name together with the NCBI id it was resolved from."""
        return f"{self.codon_table_name} (NCBI transl_table={self.transl_table})"

    @property
    def coding_offset(self) -> int:
        """0-based index of the first base belonging to a complete codon."""
        return max(self.codon_start - 1, 0)

    def _translate(self, sequence: Seq, exceptions: tuple[TranslationException, ...]) -> Seq:
        """Translate a coding sequence the way GenBank builds its ``/translation``."""
        coding_sequence = sequence[self.coding_offset:]
        complete_length = (len(coding_sequence) // 3) * 3

        try:
            residues = list(str(coding_sequence[:complete_length].translate(table=self.codon_table)))
        except Exception:
            return Seq("")

        # Translation initiates with methionine even when the start codon is an
        # alternative one (ATT, GTG, TTG ...), which is what GenBank reports.
        if complete_length >= 3 and str(coding_sequence[:3]).upper() in self.codon_table.start_codons:
            residues[0] = "M"

        for exception in exceptions:
            codon_index = (exception.start - self.coding_offset) // 3
            if 0 <= codon_index < len(residues):
                residues[codon_index] = exception.amino_acid
            elif codon_index == len(residues) and complete_length < len(coding_sequence):
                # A 1-2 base partial codon at the 3' end, completed by polyadenylation.
                residues.append(exception.amino_acid)

        return Seq("".join(residues))

    def _exceptions_for_mutated(self, mutated_sequence: Seq) -> tuple[TranslationException, ...]:
        """Re-anchor the annotated codon exceptions onto the mutated sequence.

        A 3' terminal exception marks a stop codon completed by polyadenylation of the
        mRNA, so it follows the new end of the sequence - but only while that end is
        still a partial codon. Once an indel makes the last codon complete it is read
        from the genetic code again. An internal exception is only kept when no indel
        moved it.
        """
        length_change = len(mutated_sequence) - len(self.input_sequence)
        coding_length = max(len(mutated_sequence) - self.coding_offset, 0)
        partial_bases = coding_length % 3
        exceptions = []

        for exception in self.transl_except:
            if exception.end < len(self.input_sequence):
                if length_change == 0:
                    exceptions.append(exception)
                continue

            if partial_bases:
                exceptions.append(
                    TranslationException(
                        len(mutated_sequence) - partial_bases,
                        len(mutated_sequence),
                        exception.amino_acid,
                    )
                )

        return tuple(exceptions)

    def apply_mutations(self, mutations: list[dict]):
        """Apply mutations to input sequence and store result.

        Args:
            mutations: List of dictionaries containing 'position' (1-based), 'type'
                      ('sub', 'del' or 'ins'), 'ref' (optional) and 'alt'. A deletion
                      does not need 'alt'.

        Raises:
            ValueError: If the list is empty or a mutation is malformed.
            TypeError: If a mutation, its position or its base has the wrong type.
        """
        self._validate_mutations(mutations)

        sorted_mutations = sorted(mutations, key=lambda x: x['position'], reverse=True)
        mutated_seq = MutableSeq(str(self.input_sequence))
        
        for mutation in sorted_mutations:
            position_idx = mutation['position'] - 1
            mutation_type = mutation.get('type', 'sub')
            
            if position_idx < 0:
                continue

            if mutation_type == 'sub':
                if position_idx < len(mutated_seq):
                    mutated_seq[position_idx] = mutation['alt'].upper()
            
            elif mutation_type == 'del':
                if position_idx < len(mutated_seq):
                    del mutated_seq[position_idx]
            
            elif mutation_type == 'ins':
                insert_idx = position_idx + 1
                seq_to_insert = mutation['alt'].upper()
                mutated_seq[insert_idx:insert_idx] = seq_to_insert

        self.mutated_sequence = Seq(str(mutated_seq))
        self.mutated_amino_acid_translation = self._translate(
            self.mutated_sequence, self._exceptions_for_mutated(self.mutated_sequence)
        )

    def _validate_mutations(self, mutations: list[dict]):
        """Reject malformed mutations before any of them is applied."""
        if not mutations:
            raise ValueError("Mutations list cannot be empty")

        sequence_length = len(self.input_sequence)

        for number, mutation in enumerate(mutations, start=1):
            if not isinstance(mutation, dict):
                raise TypeError(
                    f"Mutation {number} must be a dictionary, got {type(mutation).__name__}"
                )

            mutation_type = mutation.get("type", "sub")
            if mutation_type not in SUPPORTED_MUTATION_TYPES:
                raise ValueError(
                    f"Mutation {number} has an unsupported type {mutation_type!r}, "
                    f"expected one of {sorted(SUPPORTED_MUTATION_TYPES)}"
                )

            # A deletion removes a base, it never carries a replacement for it.
            required_keys = ("position",) if mutation_type == "del" else ("position", "alt")
            for key in required_keys:
                if mutation.get(key) is None:
                    raise ValueError(f"Mutation {number} is missing the {key!r} key")

            position = mutation["position"]
            if isinstance(position, bool) or not isinstance(position, int):
                raise TypeError(
                    f"Position in mutation {number} must be an integer, "
                    f"got {type(position).__name__}"
                )

            if not 1 <= position <= sequence_length:
                raise ValueError(
                    f"Position {position} in mutation {number} is out of bounds, "
                    f"the sequence has {sequence_length} nucleotides"
                )

            if mutation_type == "del":
                continue

            base = mutation["alt"]
            if not isinstance(base, str):
                raise TypeError(
                    f"Base in mutation {number} must be a string, got {type(base).__name__}"
                )

            if mutation_type == "sub" and len(base) != 1:
                raise ValueError(
                    f"Base in mutation {number} must be a single character, got {base!r}"
                )

            if not base:
                raise ValueError(f"Base in mutation {number} must not be empty")

            for nucleotide in base.upper():
                if nucleotide not in VALID_NUCLEOTIDES:
                    raise ValueError(
                        f"Base {nucleotide!r} in mutation {number} is not a valid DNA "
                        f"nucleotide, expected one of {sorted(VALID_NUCLEOTIDES)}"
                    )
