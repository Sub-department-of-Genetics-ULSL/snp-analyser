import textwrap
import unittest.mock as mock

import pytest
from Bio.Seq import Seq

from analyser_cli.fasta_reader import FastaReader


class TestFastaReader:
    """Test cases for FastaReader class."""

    @pytest.fixture
    def reader(self):
        return FastaReader()

    @pytest.fixture
    def example_fasta(self):
        return textwrap.dedent("""
    >sequence A
    ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatat
    tctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc
    >sequence B
    ggtaagtgctctagtacaaacacccccaatattgtgatataattaaaattatattcatat
    tctgttgccagattttacacttttaggctatattagagccatcttctttgaagcgttgtc
    tatgcatcgatcgacgactg
    """)
    
    def test_reading_file_ids(self, example_fasta):
        reader = FastaReader()
        with mock.patch('builtins.open', mock.mock_open(read_data=example_fasta)):
            reader.read_file("file.fasta")

            assert "sequence A" in reader.ids
            assert "sequence B" in reader.ids

    def test_reading_file_sequences(self, example_fasta):
        reader = FastaReader()
        sequence_a = Seq(
            "ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattc" \
            "atattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc"
        )

        with mock.patch('builtins.open', mock.mock_open(read_data=example_fasta)):
            reader.read_file("file.fasta")

            assert sequence_a == reader.get_sequence("sequence A")
