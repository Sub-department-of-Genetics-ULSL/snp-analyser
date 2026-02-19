import Bio
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
        >>> mutations = [{"type": "sub", "position": 3, "alt": "T"}]
        >>> analyser.apply_mutations(mutations)
        >>> print(analyser.amino_acid_translation)
        >>> print(analyser.mutated_amino_acid_translation)
    """

    def __init__(self, input_sequence: str, dna_type: str = "mitochondrial"):
        self.dna_type = dna_type
        self.input_sequence = Seq(input_sequence)
        try:
            self.amino_acid_translation = self.input_sequence.translate(
                table=MITOCHONDRIAL_TABLE if dna_type == "mitochondrial" else STANDARD_TABLE
            )
        except:
            self.amino_acid_translation = Seq("")
            
        self.mutated_sequence = None
        self.mutated_amino_acid_translation = None

    def apply_mutations(self, mutations: list[dict]):
        """Apply mutations to input sequence and store result.
        
        Args:
            mutations: List of dictionaries containing 'position' (1-based), 'type' (sub, del, ins), 
                      'ref' (optional), and 'alt' keys.
            
        Raises:
            ValueError: If mutations list is empty.
        """
        if not mutations:
            raise ValueError("Mutations list cannot be empty")
        
        sorted_mutations = sorted(mutations, key=lambda x: x['position'], reverse=True)
        mutated_seq = MutableSeq(str(self.input_sequence))
        
        for mutation in sorted_mutations:
            position_idx = mutation['position'] - 1
            mutation_type = mutation.get('type', 'sub')
            
            if position_idx < 0:
                continue

            if mutation_type == 'sub':
                if position_idx < len(mutated_seq):
                    mutated_seq[position_idx] = mutation['alt']
            
            elif mutation_type == 'del':
                if position_idx < len(mutated_seq):
                    del mutated_seq[position_idx]
            
            elif mutation_type == 'ins':
                insert_idx = position_idx + 1
                seq_to_insert = mutation['alt']
                mutated_seq[insert_idx:insert_idx] = seq_to_insert

        self.mutated_sequence = Seq(str(mutated_seq))
        
        valid_len = (len(self.mutated_sequence) // 3) * 3
        self.mutated_amino_acid_translation = self.mutated_sequence[:valid_len].translate(
            table=MITOCHONDRIAL_TABLE if self.dna_type == "mitochondrial" else STANDARD_TABLE
        )
