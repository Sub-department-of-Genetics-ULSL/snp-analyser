from pathlib import Path

from Bio import Entrez, SeqIO, SeqRecord

from .variables import (
    LATIN_TO_NCBI_ID_MAPPING,
    COMMON_TO_LATIN_MAPPING,
    DATA_FOLDER_NAME,
    ENTREZ_EMAIL,
)
from .cds_annotation import genetic_code_name, read_cds_annotation


class DataManager:
    def __init__(self):
        Entrez.email = ENTREZ_EMAIL

    @property
    def genomes(self) -> list[str]:
        return [x for x in LATIN_TO_NCBI_ID_MAPPING.keys()]
    
    @property
    def organisms(self) -> dict:
        return COMMON_TO_LATIN_MAPPING

    def download_genome(self, id: str):
        # Check if folder for data exists, if not - create it
        Path(f"{DATA_FOLDER_NAME}").mkdir(parents=True, exist_ok=True)

        with Entrez.efetch(
            db="nucleotide", rettype="gb", retmode="text", id=id
        ) as handle:
            with open(f"{DATA_FOLDER_NAME}/{id}.gb", "w") as file:
                file.write(handle.read())

    def get_genome(self, latin_name: str) -> SeqRecord.SeqRecord:
        try:
            id = LATIN_TO_NCBI_ID_MAPPING[latin_name]
        except KeyError:
            raise KeyError(f"{latin_name} not found in mapping dictionary!")
       
        file_path = f"{DATA_FOLDER_NAME}/{id}.gb"

        if not Path(file_path).is_file():
            self.download_genome(id)

        return SeqIO.read(file_path, "gb")
    
    def list_genes_for_organism(self, latin_name: str) -> list[str]:
        genome = self.get_genome(latin_name=latin_name)
        return [
            feature.qualifiers["gene"][0]
            for feature in genome.features
            if feature.type == "CDS" and feature.qualifiers.get("translation") is not None  # We only include features that are coding sequences (CDS)
        ]

    def get_gene_for_organism(self, latin_name: str, gene: str):
        genome = self.get_genome(latin_name=latin_name)
        feature = self._find_gene_feature(genome, gene)
        annotation = read_cds_annotation(feature)

        return {
            "sequence": str(feature.extract(genome.seq)),
            "startInGenome": int(feature.location.start) + 1,
            "endInGenome": int(feature.location.end),
            "translTable": annotation.transl_table,
            "translTableName": genetic_code_name(annotation.transl_table),
            "codonStart": annotation.codon_start,
            "translExcept": [
                {
                    "start": exception.start,
                    "end": exception.end,
                    "aminoAcid": exception.amino_acid,
                }
                for exception in annotation.transl_except
            ],
            "ncbiTranslation": annotation.translation,
        }

    @staticmethod
    def _find_gene_feature(genome: SeqRecord.SeqRecord, gene: str):
        """Return the feature describing a gene, preferring its CDS.

        A gene is usually annotated twice: once as a ``gene`` feature and once as a
        ``CDS``. Only the CDS carries /transl_table, /codon_start and /transl_except,
        so it has to win.
        """
        fallback = None

        for feature in genome.features:
            names = feature.qualifiers.get("gene")
            if not names or names[0] != gene:
                continue
            if feature.type == "CDS":
                return feature
            if fallback is None:
                fallback = feature

        if fallback is not None:
            return fallback

        raise ValueError(f"Gene {gene} not found in genome!")
