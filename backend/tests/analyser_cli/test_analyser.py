import pytest
from Bio.Seq import Seq
from analyser_cli.analyser import Analyser, get_codon_table

MITOCHONDRIAL_TABLE = get_codon_table(2)
STANDARD_TABLE = get_codon_table(1)


@pytest.fixture
def simple_sequence():
    """Fixture providing a simple DNA sequence for testing."""
    return "ATGCGATCG"


@pytest.fixture
def long_sequence():
    """Fixture providing a longer DNA sequence for testing."""
    return "ATGCGATCGTATGCGATCGTATGCGATCGT"


@pytest.fixture
def sequence_with_stop():
    """Fixture providing a sequence with stop codon."""
    return "ATGTAA"


@pytest.fixture
def mitochondrial_analyser(simple_sequence):
    """Fixture providing a mitochondrial Analyser instance."""
    return Analyser(simple_sequence, dna_type="mitochondrial", transl_table=2)


@pytest.fixture
def nuclear_analyser(simple_sequence):
    """Fixture providing a nuclear Analyser instance."""
    return Analyser(simple_sequence, dna_type="nuclear", transl_table=1)


@pytest.fixture
def extended_analyser():
    """Fixture providing an extended Analyser instance for testing."""
    return Analyser("ATGCGATCGTAA", dna_type="mitochondrial")


@pytest.fixture
def multiple_mutations():
    """Fixture providing multiple mutations for testing."""
    return [
        {"position": 3, "alt": "T"},
        {"position": 6, "alt": "C"},
        {"position": 9, "alt": "A"}
    ]


@pytest.fixture
def single_mutation():
    """Fixture providing a single mutation for testing."""
    return [{"position": 3, "alt": "T"}]


@pytest.fixture
def invalid_mutations():
    """Fixture providing various invalid mutation formats for testing."""
    return {
        "empty_list": [],
        "non_dict": ["invalid", {"position": 1, "alt": "T"}],
        "missing_position": [{"alt": "T"}],
        "missing_base": [{"position": 1}],
        "invalid_position_type": [{"position": "1", "alt": "T"}],
        "position_too_low": [{"position": 0, "alt": "T"}],
        "position_too_high": [{"position": 10, "alt": "T"}],
        "invalid_base_type": [{"position": 1, "alt": 123}],
        "invalid_base_length": [{"position": 1, "alt": "AT"}],
        "invalid_base_character": [{"position": 1, "alt": "X"}]
    }


class TestAnalyserInitialization:
    """Test cases for Analyser class initialization."""
    
    def test_init_mitochondrial_default(self, simple_sequence):
        """Test initialization with mitochondrial DNA (default)."""
        analyser = Analyser(simple_sequence)
        
        assert analyser.dna_type == "mitochondrial"
        assert str(analyser.input_sequence) == simple_sequence
        assert analyser.mutated_sequence is None
        assert analyser.mutated_amino_acid_translation is None
        assert analyser.amino_acid_translation is not None
    
    def test_init_mitochondrial_explicit(self, simple_sequence):
        """Test initialization with explicitly specified mitochondrial DNA."""
        analyser = Analyser(simple_sequence, dna_type="mitochondrial")
        
        assert analyser.dna_type == "mitochondrial"
        assert str(analyser.input_sequence) == simple_sequence
    
    def test_init_nuclear(self, simple_sequence):
        """Test initialization with nuclear DNA."""
        analyser = Analyser(simple_sequence, dna_type="nuclear")
        
        assert analyser.dna_type == "nuclear"
        assert str(analyser.input_sequence) == simple_sequence
    
    def test_amino_acid_translation_mitochondrial(self, mitochondrial_analyser, simple_sequence):
        """Test amino acid translation using mitochondrial codon table."""
        expected = Seq(simple_sequence).translate(table=MITOCHONDRIAL_TABLE)
        assert str(mitochondrial_analyser.amino_acid_translation) == str(expected)
    
    def test_amino_acid_translation_nuclear(self, nuclear_analyser, simple_sequence):
        """Test amino acid translation using standard codon table."""
        expected = Seq(simple_sequence).translate(table=STANDARD_TABLE)
        assert str(nuclear_analyser.amino_acid_translation) == str(expected)


class TestApplyMutations:
    """Test cases for the apply_mutations method."""
    
    def test_apply_single_mutation(self, mitochondrial_analyser, single_mutation):
        """Test applying a single mutation."""
        mitochondrial_analyser.apply_mutations(single_mutation)
        
        assert str(mitochondrial_analyser.mutated_sequence) == "ATTCGATCG"
        assert mitochondrial_analyser.mutated_amino_acid_translation is not None
    
    def test_apply_multiple_mutations(self, mitochondrial_analyser, multiple_mutations):
        """Test applying multiple mutations."""
        mitochondrial_analyser.apply_mutations(multiple_mutations)
        
        assert str(mitochondrial_analyser.mutated_sequence) == "ATTCGCTCA"
    
    def test_mutations_case_insensitive(self, mitochondrial_analyser):
        """Test that mutations handle lowercase bases correctly."""
        mutations = [{"position": 3, "alt": "t"}]
        mitochondrial_analyser.apply_mutations(mutations)
        
        assert str(mitochondrial_analyser.mutated_sequence) == "ATTCGATCG"
    
    def test_mutation_amino_acid_translation(self):
        """Test that mutations properly update amino acid translation."""
        sequence = "ATGAAATAA"  # M-K-* (mitochondrial stop)
        analyser = Analyser(sequence, dna_type="mitochondrial")
        original_translation = str(analyser.amino_acid_translation)
        
        # Change middle codon AAA to TTT
        mutations = [
            {"position": 4, "alt": "T"},
            {"position": 5, "alt": "T"},
            {"position": 6, "alt": "T"}
        ]
        
        analyser.apply_mutations(mutations)
        mutated_translation = str(analyser.mutated_amino_acid_translation)
        
        assert original_translation != mutated_translation
        assert "F" in mutated_translation  # TTT codes for Phenylalanine


class TestMutationValidation:
    """Test cases for mutation validation."""
    
    def test_empty_mutations_list(self, mitochondrial_analyser, invalid_mutations):
        """Test that empty mutations list raises ValueError."""
        with pytest.raises(ValueError, match="Mutations list cannot be empty"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["empty_list"])
    
    def test_none_mutations_list(self, mitochondrial_analyser, invalid_mutations):
        """Test that None mutations list raises ValueError."""
        with pytest.raises(ValueError, match="Mutations list cannot be empty"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["empty_list"])
    
    def test_invalid_mutation_format(self, mitochondrial_analyser, invalid_mutations):
        """Test that non-dictionary mutations raise TypeError."""
        with pytest.raises(TypeError, match="Mutation 1 must be a dictionary"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["non_dict"])
    
    def test_missing_position_key(self, mitochondrial_analyser, invalid_mutations):
        """Test that mutations missing 'position' key raise ValueError."""
        with pytest.raises(ValueError, match="Mutation 1 is missing the 'position' key"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["missing_position"])
    
    def test_missing_alt_key(self, mitochondrial_analyser, invalid_mutations):
        """Test that mutations missing the 'alt' key raise ValueError."""
        with pytest.raises(ValueError, match="Mutation 1 is missing the 'alt' key"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["missing_base"])
    
    def test_invalid_position_type(self, mitochondrial_analyser, invalid_mutations):
        """Test that non-integer positions raise TypeError."""
        with pytest.raises(TypeError, match="Position in mutation 1 must be an integer"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["invalid_position_type"])
    
    def test_position_out_of_bounds_low(self, mitochondrial_analyser, invalid_mutations):
        """Test that positions < 1 raise ValueError."""
        with pytest.raises(ValueError, match="Position 0 in mutation 1 is out of bounds"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["position_too_low"])
    
    def test_position_out_of_bounds_high(self, mitochondrial_analyser, invalid_mutations):
        """Test that positions > sequence length raise ValueError."""
        with pytest.raises(ValueError, match="Position 10 in mutation 1 is out of bounds"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["position_too_high"])
    
    def test_invalid_base_type(self, mitochondrial_analyser, invalid_mutations):
        """Test that non-string bases raise TypeError."""
        with pytest.raises(TypeError, match="Base in mutation 1 must be a string"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["invalid_base_type"])
    
    def test_invalid_base_length(self, mitochondrial_analyser, invalid_mutations):
        """Test that multi-character bases raise ValueError."""
        with pytest.raises(ValueError, match="Base in mutation 1 must be a single character"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["invalid_base_length"])
    
    def test_invalid_base_character(self, mitochondrial_analyser, invalid_mutations):
        """Test that invalid DNA bases raise ValueError."""
        with pytest.raises(ValueError, match="Base 'X' in mutation 1 is not a valid DNA nucleotide"):
            mitochondrial_analyser.apply_mutations(invalid_mutations["invalid_base_character"])

    def test_deletion_does_not_need_a_base(self, mitochondrial_analyser):
        """A deletion removes a base, so it carries no replacement for it."""
        mitochondrial_analyser.apply_mutations([{"type": "del", "position": 3}])

        assert str(mitochondrial_analyser.mutated_sequence) == "ATCGATCG"

    def test_insertion_accepts_several_nucleotides(self, mitochondrial_analyser):
        """An insertion is the only type whose base may be longer than one character."""
        mitochondrial_analyser.apply_mutations([{"type": "ins", "position": 3, "alt": "TTT"}])

        assert str(mitochondrial_analyser.mutated_sequence) == "ATGTTTCGATCG"

    def test_unsupported_mutation_type_raises(self, mitochondrial_analyser):
        """Test that a type the analyser cannot apply is rejected instead of ignored."""
        with pytest.raises(ValueError, match="Mutation 1 has an unsupported type 'dup'"):
            mitochondrial_analyser.apply_mutations([{"type": "dup", "position": 1, "alt": "A"}])

    def test_nothing_is_applied_when_a_later_mutation_is_invalid(self, mitochondrial_analyser):
        """The whole batch is validated before the first mutation is applied."""
        with pytest.raises(ValueError):
            mitochondrial_analyser.apply_mutations(
                [{"position": 1, "alt": "T"}, {"position": 99, "alt": "G"}]
            )

        assert mitochondrial_analyser.mutated_sequence is None


class TestDNATypes:
    """Test cases for different DNA types."""
    
    def test_mitochondrial_vs_nuclear_translation(self):
        """Test that mitochondrial and nuclear DNA types produce different translations."""
        # UGA is a stop codon in standard code but codes for Trp in mitochondrial
        sequence = "TGATGATGA"
        
        mito_analyser = Analyser(sequence, dna_type="mitochondrial", transl_table=2)
        nuclear_analyser = Analyser(sequence, dna_type="nuclear", transl_table=1)
        
        mito_translation = str(mito_analyser.amino_acid_translation)
        nuclear_translation = str(nuclear_analyser.amino_acid_translation)
        
        # Should be different due to different codon tables
        assert mito_translation != nuclear_translation


class TestEdgeCases:
    """Test cases for edge cases and boundary conditions."""
    
    def test_single_nucleotide_sequence(self, single_mutation):
        """Test with a single nucleotide sequence."""
        analyser = Analyser("A")
        mutations = [{"position": 1, "alt": "T"}]
        
        analyser.apply_mutations(mutations)
        
        assert str(analyser.mutated_sequence) == "T"
    
    def test_sequence_with_stop_codon(self, sequence_with_stop):
        """Test sequence that contains stop codons."""
        analyser = Analyser(sequence_with_stop)
        
        assert "*" in str(analyser.amino_acid_translation)
    
    def test_mutation_creating_stop_codon(self):
        """Test mutation that creates a stop codon."""
        sequence = "ATGAAA"  # M-K
        analyser = Analyser(sequence, dna_type="nuclear")
        
        # Change AAA to TAA (stop codon)
        mutations = [{"position": 4, "alt": "T"}]
        analyser.apply_mutations(mutations)
        
        assert "*" in str(analyser.mutated_amino_acid_translation)
    
    def test_multiple_mutations_same_position(self, mitochondrial_analyser):
        """Test multiple mutations at the same position (last one should win)."""
        mutations = [
            {"position": 3, "alt": "A"},
            {"position": 3, "alt": "T"}  # This should be the final base
        ]
        
        mitochondrial_analyser.apply_mutations(mutations)
        
        assert str(mitochondrial_analyser.mutated_sequence) == "ATTCGATCG"


class TestComplexScenarios:
    """Test cases for complex mutation scenarios."""
    
    def test_all_possible_bases(self):
        """Test mutations with all valid DNA bases."""
        sequence = "AAAAAAAA"
        analyser = Analyser(sequence)
        mutations = [
            {"position": 1, "alt": "A"},
            {"position": 2, "alt": "T"},
            {"position": 3, "alt": "G"},
            {"position": 4, "alt": "C"}
        ]
        
        analyser.apply_mutations(mutations)
        
        assert str(analyser.mutated_sequence)[:4] == "ATGC"
    
    def test_long_sequence_many_mutations(self, long_sequence):
        """Test with a longer sequence and many mutations."""
        analyser = Analyser(long_sequence)
        
        # Apply mutations at every 3rd position
        mutations = [
            {"position": i, "alt": "A"} 
            for i in range(3, len(long_sequence) + 1, 3)
        ]
        analyser.apply_mutations(mutations)
        # Check that mutations were applied
        mutated = str(analyser.mutated_sequence)
        for i in range(2, len(long_sequence), 3):  # 0-based indexing
            if i < len(mutated):
                assert mutated[i] == "A"


class TestWithFixtures:
    """Test cases using pytest fixtures."""
    
    def test_analyser_fixture(self, extended_analyser):
        """Test using the extended analyser fixture."""
        assert extended_analyser.dna_type == "mitochondrial"
        assert len(str(extended_analyser.input_sequence)) == 12
    
    def test_mutations_fixture(self, extended_analyser, multiple_mutations):
        """Test applying sample mutations using fixtures."""
        extended_analyser.apply_mutations(multiple_mutations)
        
        assert extended_analyser.mutated_sequence is not None
        assert extended_analyser.mutated_amino_acid_translation is not None
        assert str(extended_analyser.mutated_sequence) == "ATTCGCTCATAA"
