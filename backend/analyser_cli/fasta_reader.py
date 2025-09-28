from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser


class FastaReader:
    """Class for reading and containing genomic sequences from FASTA file."""
    
    def __init__(self):
        self._sequences = {}

    def read_file(self, filename: str):
        self._sequences = {}
        
        with open(filename) as handle:
            for values in SimpleFastaParser(handle):
                id = values[0]
                sequence = values[1]
                self._sequences[id] = Seq(sequence)

    @property
    def ids(self) -> list[str]:
        return [x for x in self._sequences.keys()]
    
    def get_sequence(self, id: str) -> Seq|None:
        return self._sequences.get(id)
