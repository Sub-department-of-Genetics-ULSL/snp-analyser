from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient

from backend.analyser_backend.models.report_data import ReportData
from backend.analyser_backend.routers import report
from backend.analyser_cli import Analyser
from backend.analyser_cli.data_manager import DataManager
from backend.analyser_cli.hgvs import ReferenceContext, describe_variants, reverse_complement

# Codons M K Q F G * on genome positions 11..28.
CDS = "ATGAAACAGTTCGGATAA"
FLANK_5, FLANK_3 = "GCGCGCGCGC", "CCGCCGCCGC"


def make_record(strand: int = 1) -> SeqRecord:
    gene = CDS if strand == 1 else reverse_complement(CDS)
    record = SeqRecord(
        Seq(FLANK_5 + gene + FLANK_3), id="NC_000999.1", annotations={"topology": "circular"}
    )
    record.features = [
        SeqFeature(
            FeatureLocation(10, 28, strand=strand),
            type="CDS",
            qualifiers={"gene": ["ND2"], "transl_table": ["2"], "protein_id": ["YP_000999.1"]},
        )
    ]
    return record


class FakeManager(DataManager):
    def __init__(self, record: SeqRecord):
        super().__init__()
        self.record = record

    def get_genome(self, latin_name: str):
        if latin_name != "test_organism":
            raise KeyError(latin_name)
        return self.record


@pytest.fixture
def client(monkeypatch):
    monkeypatch.setattr(report, "manager", FakeManager(make_record(strand=-1)))
    app = FastAPI()
    app.include_router(report.router, prefix="/api")
    return TestClient(app)


def test_get_gene_for_organism_returns_ncbi_1_based_start(monkeypatch):
    genome = SeqRecord(Seq("AAATGCGTAA"), id="test")
    genome.features = [
        SeqFeature(
            FeatureLocation(3, 8),
            type="gene",
            qualifiers={"gene": ["ND2"]},
        )
    ]

    manager = DataManager()
    monkeypatch.setattr(manager, "get_genome", lambda latin_name: genome)

    gene_data = manager.get_gene_for_organism("canis_lupus_familiaris", "ND2")

    assert gene_data["sequence"] == "TGCGT"
    assert gene_data["startInGenome"] == 4
    assert gene_data["endInGenome"] == 8
    assert gene_data["strand"] == 1


def test_get_gene_for_organism_reads_a_complement_location():
    gene_data = FakeManager(make_record(strand=-1)).get_gene_for_organism("test_organism", "ND2")

    # The sequence is the CDS in transcription order, c.1 is the highest genome position.
    assert gene_data["sequence"] == CDS
    assert gene_data["strand"] == -1
    assert (gene_data["startInGenome"], gene_data["endInGenome"]) == (11, 28)
    assert (gene_data["cdsStartInGenome"], gene_data["cdsEndInGenome"]) == (28, 11)
    assert gene_data["accession"] == "NC_000999.1"
    assert gene_data["proteinId"] == "YP_000999.1"
    assert gene_data["codingNumberingAvailable"] is True


def _mutated(context: ReferenceContext, *descriptions: str) -> str:
    result = describe_variants(context, list(descriptions))
    assert not result.errors, result.errors
    analyser = Analyser(context.cds_sequence, transl_table=2)
    analyser.apply_mutations(report._analyser_mutations(result.variants))
    return str(analyser.mutated_sequence)


@pytest.mark.parametrize(
    "descriptions, expected",
    [
        # A deletion-insertion replaces the bases in place.
        (["c.4_5delinsGT"], "ATGGTACAGTTCGGATAA"),
        (["c.1_2delinsCC"], "CCGAAACAGTTCGGATAA"),
        (["c.4_6delinsT"], "ATGTCAGTTCGGATAA"),
        (["c.7_9dup"], "ATGAAACAGCAGTTCGGATAA"),
        (["c.18_*1insG"], "ATGAAACAGTTCGGATAAG"),
        (["c.4del", "c.10T>C", "c.12_13insAAA"], "ATGAACAGCTCAAAGGATAA"),
    ],
)
def test_analyser_receives_the_variants_in_place(descriptions, expected):
    context = ReferenceContext.from_genbank(make_record(), make_record().features[0])
    assert _mutated(context, *descriptions) == expected


def test_report_components_describe_complement_strand_variants(monkeypatch):
    monkeypatch.setattr(report, "manager", FakeManager(make_record(strand=-1)))
    data = ReportData(organism="test_organism", gene="ND2", mutations=["m.25T>C"])

    report_text, report_mutations, analyser, _, gene_context = report._generate_report_components(data)

    assert str(analyser.mutated_sequence) == "ATGGAACAGTTCGGATAA"
    assert report_mutations[0]["hgvs_m"] == "NC_000999.1:m.25T>C"
    assert report_mutations[0]["hgvs_c"] == "NC_000999.1(ND2):c.4A>G"
    assert report_mutations[0]["hgvs_p"] == "YP_000999.1:p.(Lys2Glu)"
    assert "HGVS (coding DNA): NC_000999.1(ND2):c.4A>G" in report_text
    assert "HGVS (protein): YP_000999.1:p.(Lys2Glu)" in report_text
    assert gene_context["context"].strand == -1


def test_report_components_reject_invalid_descriptions(monkeypatch):
    monkeypatch.setattr(report, "manager", FakeManager(make_record()))
    data = ReportData(organism="test_organism", gene="ND2", mutations=["c.4C>G"])

    with pytest.raises(ValueError, match="is A, not C"):
        report._generate_report_components(data)


def test_frameshift_residue_is_the_first_changed_amino_acid():
    context = ReferenceContext.from_genbank(make_record(), make_record().features[0])
    result = describe_variants(context, ["c.4_5insT", "c.13G>T"])

    assert report._frameshift_residue(result.variants) == result.variants[0].protein.first_residue
    assert report._frameshift_residue(describe_variants(context, ["c.13G>T"]).variants) is None


def test_describe_endpoint_returns_normalised_descriptions(client):
    response = client.post(
        "/api/hgvs/describe",
        json={"organism": "test_organism", "gene": "ND2", "mutations": ["m.25T>C", "c.4a>g"]},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["hgvsVersion"] == "21.1"
    assert payload["variants"][0]["c"] == "NC_000999.1(ND2):c.4A>G"
    assert payload["variants"][0]["m"] == "NC_000999.1:m.25T>C"
    assert payload["variants"][0]["inputs"] == [0]
    assert payload["errors"][0]["inputs"] == [1]
    assert "case-sensitive" in payload["errors"][0]["message"]
    assert payload["reference"]["strand"] == -1


@pytest.mark.parametrize(
    "organism, gene", [("unknown_organism", "ND2"), ("test_organism", "ND9")]
)
def test_describe_endpoint_reports_unknown_genes(client, organism, gene):
    response = client.post(
        "/api/hgvs/describe", json={"organism": organism, "gene": gene, "mutations": ["c.4A>G"]}
    )
    assert response.status_code == 404


@pytest.mark.parametrize(
    "mutations, message",
    [
        ([], "At least one variant"),
        (["c.4C>G"], "is A, not C"),
        (["c.4_9insA"], "flanking"),
        (["c.4A>G", "c.4del"], "overlap"),
    ],
)
def test_report_job_is_not_queued_for_invalid_descriptions(client, mutations, message):
    response = client.post(
        "/api/report/jobs",
        json={"organism": "test_organism", "gene": "ND2", "mutations": mutations},
    )

    assert response.status_code == 400
    assert message in response.json()["detail"]
