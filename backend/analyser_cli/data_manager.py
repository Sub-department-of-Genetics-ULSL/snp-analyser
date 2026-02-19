from pathlib import Path

from Bio import Entrez, SeqIO, SeqRecord

from .variables import LATIN_TO_NCBI_ID_MAPPING, COMMON_TO_LATIN_MAPPING, DATA_FOLDER_NAME, ENTREZ_EMAIL


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
            if feature.type == "gene" and feature.qualifiers.get("gene") is not None
        ]

    def get_gene_for_organism(self, latin_name: str, gene: str):
        genome = self.get_genome(latin_name=latin_name)
        for feature in genome.features:
            if feature.qualifiers.get("gene") is not None and feature.qualifiers["gene"][0] == gene:
                return {"sequence": str(feature.extract(genome.seq)), "startInGenome": int(feature.location.start)}
            
        raise ValueError(f"Gene {gene} not found in genome!")
