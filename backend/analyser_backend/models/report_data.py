from pydantic import BaseModel


class ReportData(BaseModel):
    organism: str
    gene: str
    mutations: list[str]
    predict_mutated_structure: bool = False
