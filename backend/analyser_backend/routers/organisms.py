import sys
sys.path.append("../..")

from fastapi import APIRouter

from backend.analyser_cli import DataManager


router = APIRouter()
manager = DataManager()

@router.get("/organisms/", tags=["organisms"])
def get_available_organisms():
    return manager.organisms

@router.get("/organisms/{latin_name}/", tags=["organisms"])
def get_genes_for_organism(latin_name: str):
    return {"organism": latin_name, "genes": manager.list_genes_for_organism(latin_name=latin_name)}

@router.get("/organisms/{latin_name}/{gene}/", tags=["organisms"])
def get_gene_for_organism(latin_name: str, gene: str):
    return {"organism": latin_name, "gene": gene} | manager.get_gene_for_organism(latin_name=latin_name, gene=gene)
