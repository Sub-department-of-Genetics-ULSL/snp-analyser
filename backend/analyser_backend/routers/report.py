import re
import sys
sys.path.append("../..")

from fastapi import APIRouter

from models.report_data import ReportData
from backend.analyser_cli import Reporter, Analyser, DataManager

router = APIRouter()
manager = DataManager()
reporter = Reporter()

@router.post("/report/", tags=["report"])
def generate_report(data: ReportData):
    gene_data = manager.get_gene_for_organism(data.organism, data.gene)
    gene_sequence = gene_data["sequence"]
    start_genome_pos = gene_data.get("startInGenome", 1)
    
    analyser = Analyser(gene_sequence)
    
    parsed_mutations = []
    
    sub_pattern = re.compile(r"(\d+)([ACGT])>([ACGT])")
    del_pattern = re.compile(r"(?:del\.(\d+)([ACGT])|(\d+)del)")
    ins_pattern = re.compile(r"ins\.(\d+)(?:_\d+)?([ACGT]+)")

    for mutation in data.mutations:
        m_sub = sub_pattern.search(mutation)
        m_del = del_pattern.search(mutation)
        m_ins = ins_pattern.search(mutation)

        if m_sub:
            abs_pos = int(m_sub.group(1))
            rel_pos = abs_pos - start_genome_pos + 1
            parsed_mutations.append({
                "type": "sub",
                "position": rel_pos,
                "ref": m_sub.group(2),
                "alt": m_sub.group(3)
            })
        elif m_del:
            pos_str = m_del.group(1) if m_del.group(1) else m_del.group(3)
            abs_pos = int(pos_str)
            rel_pos = abs_pos - start_genome_pos + 1
            parsed_mutations.append({
                "type": "del",
                "position": rel_pos,
                "ref": m_del.group(2) if m_del.group(2) else None,
                "alt": None
            })
        elif m_ins:
            abs_pos = int(m_ins.group(1))
            rel_pos = abs_pos - start_genome_pos + 1
            parsed_mutations.append({
                "type": "ins",
                "position": rel_pos,
                "ref": None,
                "alt": m_ins.group(2)
            })

    analyser.apply_mutations(parsed_mutations)

    report = reporter.generate_report(analyser, parsed_mutations)
    return report
