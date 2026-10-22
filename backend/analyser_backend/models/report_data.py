from pydantic import BaseModel


class ReportData(BaseModel):
    organism: str
    gene: str
    mutations: list[str]
    predict_mutated_structure: bool = False


class DescribeRequest(BaseModel):
    """HGVS descriptions to validate and normalise for one gene of an organism."""

    organism: str
    gene: str
    mutations: list[str]
