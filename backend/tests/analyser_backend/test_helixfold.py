from backend.analyser_backend.services import helixfold as helixfold_service
from backend.analyser_backend.services.helixfold import (
    backbone_geometry,
    build_prediction_cache_key,
    normalize_mutations,
    protein_sequence_for_prediction,
    HelixFoldService,
    transform_pdb_coordinates,
)
from pathlib import Path


def test_normalize_mutations_is_order_independent():
    mutations = [
        {"type": "ins", "position": 5, "end_position": 5, "alt": "AA"},
        {"type": "sub", "position": 2, "ref": "A", "alt": "G"},
    ]

    reversed_mutations = list(reversed(mutations))

    assert normalize_mutations(mutations) == normalize_mutations(reversed_mutations)


def test_prediction_cache_key_changes_with_mutations():
    base_kwargs = {
        "organism": "human",
        "gene": "CYTB",
        "gene_sequence": "ATGAAA",
        "protein_sequence": "MK",
        "model_path": __import__("pathlib").Path("/models/helixfold.pdparams"),
        "script_path": __import__("pathlib").Path("/opt/HelixFold-single/helixfold_single_inference.py"),
    }

    first = build_prediction_cache_key(
        mutations=[{"type": "sub", "position": 2, "ref": "T", "alt": "C"}],
        **base_kwargs,
    )
    second = build_prediction_cache_key(
        mutations=[{"type": "sub", "position": 2, "ref": "T", "alt": "A"}],
        **base_kwargs,
    )

    assert first != second


def test_protein_sequence_for_prediction_trims_stop():
    assert protein_sequence_for_prediction("MKT*AA") == "MKT"


def test_transform_pdb_coordinates_updates_atom_lines():
    pdb_text = (
        "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C\n"
        "END\n"
    )

    transformed = transform_pdb_coordinates(
        pdb_text,
        rotation=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        translation=(1.0, 2.0, 3.0),
    )

    assert "  11.000  22.000  33.000" in transformed


def test_resolve_python_bin_prefers_helixfold_venv(tmp_path, monkeypatch):
    repo_dir = tmp_path / "HelixFold-single"
    repo_python = repo_dir / ".venv" / "bin" / "python"
    repo_python.parent.mkdir(parents=True)
    repo_python.write_text("", encoding="utf-8")

    monkeypatch.setenv("HELIXFOLD_SINGLE_REPO_DIR", str(repo_dir))
    monkeypatch.setenv("HELIXFOLD_SINGLE_PYTHON_BIN", "/usr/bin/python3")

    service = HelixFoldService()

    assert service._python_bin == str(repo_python.resolve())


def test_backbone_geometry_flags_unrelaxed_prediction():
    broken = (
        "ATOM      1  N   MET A   1      12.746 -16.314  40.391  1.00 51.62           N\n"
        "ATOM      2  CA  MET A   1      12.349 -16.103  39.002  1.00 51.62           C\n"
        "ATOM      3  C   MET A   1      11.187 -15.120  38.912  1.00 51.62           C\n"
        "ATOM      9  N   PRO A   2      11.419 -15.100  38.292  1.00 52.87           N\n"
        "ATOM     10  CA  PRO A   2      10.280 -14.184  38.200  1.00 52.87           C\n"
        "ATOM     11  C   PRO A   2       9.042 -14.840  37.592  1.00 52.87           C\n"
    )

    report = backbone_geometry(broken)

    assert report["measured"] is True
    assert report["residues"] == 2
    assert report["peptide_bond_mean"] < 1.0
    assert report["ca_distance_mean"] < 3.5


def test_backbone_geometry_accepts_physical_chain():
    physical = (
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 90.00           N\n"
        "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 90.00           C\n"
        "ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 90.00           C\n"
        "ATOM      4  N   ALA A   2       3.339   1.420   0.000  1.00 90.00           N\n"
        "ATOM      5  CA  ALA A   2       4.132   2.640   0.000  1.00 90.00           C\n"
        "ATOM      6  C   ALA A   2       5.590   2.640   0.000  1.00 90.00           C\n"
    )

    report = backbone_geometry(physical)

    assert 1.29 <= report["peptide_bond_mean"] <= 1.37
    assert 3.7 <= report["ca_distance_mean"] <= 3.9


def test_backbone_geometry_handles_single_residue():
    single = "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C\n"

    assert backbone_geometry(single)["measured"] is False


def test_relax_structure_skipped_when_disabled(tmp_path, monkeypatch):
    monkeypatch.setenv("HELIXFOLD_RELAX_ENABLED", "false")
    service = HelixFoldService(cache_root=tmp_path / "cache")

    info = service._relax_structure(tmp_path / "unrelaxed.pdb", tmp_path / "relaxed.pdb")

    assert info["relaxed"] is False
    assert "skipped" in info


def test_prediction_cache_key_depends_on_pipeline_version(monkeypatch):
    kwargs = {
        "organism": "human",
        "gene": "CYTB",
        "gene_sequence": "ATGAAA",
        "protein_sequence": "MK",
        "mutations": [{"type": "sub", "position": 2, "ref": "T", "alt": "C"}],
        "model_path": Path("/models/helixfold.pdparams"),
        "script_path": Path("/opt/HelixFold-single/helixfold_single_inference.py"),
    }

    before = build_prediction_cache_key(**kwargs)
    monkeypatch.setattr(helixfold_service, "PIPELINE_VERSION", 99)

    assert build_prediction_cache_key(**kwargs) != before
