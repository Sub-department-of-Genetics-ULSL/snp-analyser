import Bio
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.Data import CodonTable


MITOCHONDRIAL_TABLE = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
STANDARD_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]


class Analyser:
    """Analyser class for processing DNA sequences and their amino acid translations.
    
    This class provides functionality to analyze DNA sequences (both mitochondrial and nuclear)
    and translate them to amino acid sequences using appropriate codon tables. It supports
    applying mutations to sequences and comparing original vs mutated amino acid translations.
    
    Attributes:
        dna_type (str): Type of DNA sequence ('mitochondrial' or 'nuclear').
        input_sequence (str): Original DNA sequence.
        amino_acid_translation (Bio.Seq.Seq): Amino acid translation of original sequence.
        mutated_sequence (Bio.Seq.Seq or None): DNA sequence after applying mutations.
        mutated_amino_acid_translation (Bio.Seq.Seq or None): Amino acid translation of mutated sequence.
    
    Example:
        >>> analyser = Analyser("ATGCGATCG", dna_type="mitochondrial")
        >>> mutations = [{"position": 3, "base": "T"}]
        >>> analyser.apply_mutations(mutations)
        >>> print(analyser.amino_acid_translation)
        >>> print(analyser.mutated_amino_acid_translation)
    """

    def __init__(self, input_sequence: str, dna_type: str = "mitochondrial"):
        self.dna_type = dna_type
        self.input_sequence = Seq(input_sequence)
        self.amino_acid_translation = self.input_sequence.translate(
            table=MITOCHONDRIAL_TABLE if dna_type == "mitochondrial" else STANDARD_TABLE
        )
        self.mutated_sequence = None
        self.mutated_amino_acid_translation = None

    def apply_mutations(self, mutations: list[dict]):
        """Apply mutations to input sequence and store result.
        
        Args:
            mutations: List of dictionaries containing 'position' (1-based) and 'base' keys.
            
        Raises:
            ValueError: If mutations list is empty, positions are invalid, or bases are invalid.
            TypeError: If mutation format is incorrect.
        """
        if not mutations:
            raise ValueError("Mutations list cannot be empty")
        
        # Validate mutations format and content
        valid_bases = set('ATGC')
        sequence_length = len(self.input_sequence)
        
        for i, mutation in enumerate(mutations):
            if not isinstance(mutation, dict):
                raise TypeError(f"Mutation {i+1} must be a dictionary")
            
            if 'position' not in mutation or 'base' not in mutation:
                raise ValueError(f"Mutation {i+1} must contain 'position' and 'base' keys")
            
            position = mutation['position']
            if not isinstance(position, int):
                raise TypeError(f"Position in mutation {i+1} must be an integer")
            if position < 1 or position > sequence_length:
                raise ValueError(f"Position {position} in mutation {i+1} is out of bounds. "
                               f"Sequence length is {sequence_length}")
            
            base = mutation['base']
            if not isinstance(base, str):
                raise TypeError(f"Base in mutation {i+1} must be a string")
            if len(base) != 1:
                raise ValueError(f"Base in mutation {i+1} must be a single character")
            if base.upper() not in valid_bases:
                raise ValueError(f"Base '{base}' in mutation {i+1} is not a valid DNA nucleotide (A, T, G, C)")
        
        # Apply mutations if all validations pass
        mutated_seq = MutableSeq(str(self.input_sequence))
        for mutation in mutations:
            pos = mutation["position"] - 1  # Given positions start at 1
            new_base = mutation["base"].upper()
            mutated_seq[pos] = new_base
        self.mutated_sequence = Seq(str(mutated_seq))
        self.mutated_amino_acid_translation = self.mutated_sequence.translate(
            table=MITOCHONDRIAL_TABLE if self.dna_type == "mitochondrial" else STANDARD_TABLE
        )
