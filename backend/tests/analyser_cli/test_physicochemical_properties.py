from types import SimpleNamespace

import pytest
from rpy2.rinterface_lib.embedded import RRuntimeError

from backend.analyser_cli.sub_reporters.physicochemical_properties import PhysicochemicalProperties


class _FakePeptides:
    """Stand-in for the R package, mirroring its minimum length requirements."""

    def __init__(self):
        self.seen = []
        self.seen_scale = None

    def _value(self, sequence: str):
        self.seen.append(sequence)
        return [1.2344 if sequence == "MKS" else 1.23449]

    def pI(self, sequence: str, pKscale: str = "EMBOSS"):
        self.seen_scale = pKscale
        return self._value(sequence)

    def mw(self, sequence: str):
        return self._value(sequence)

    def charge(self, sequence: str):
        return self._value(sequence)

    def hydrophobicity(self, sequence: str, scale: str):
        return self._value(sequence)

    def instaIndex(self, sequence: str):
        if len(sequence) < 3:
            raise RRuntimeError("dim(X) must have a positive length")
        return self._value(sequence)


@pytest.fixture
def fake_peptides(monkeypatch):
    peptides = _FakePeptides()
    monkeypatch.setattr(
        "backend.analyser_cli.sub_reporters.physicochemical_properties.importr",
        lambda package_name: peptides,
    )
    return peptides


def _analyser(original_translation: str, mutated_translation: str) -> SimpleNamespace:
    return SimpleNamespace(
        input_sequence="ATGAAAAGTTAA",
        mutated_sequence="ATGAGAAGTTAA",
        amino_acid_translation=original_translation,
        mutated_amino_acid_translation=mutated_translation,
    )


def test_generate_content_rounds_to_three_decimals_before_comparing(fake_peptides):
    reporter = PhysicochemicalProperties()

    content = reporter.generate_content(_analyser("MKS*", "MRS*"), [{"position": 5, "alt": "G"}])

    assert "Original Sequence: 1.234" in content
    assert "Mutated Sequence: 1.234" in content
    assert content.count("Change: -") == 5


def test_properties_are_computed_on_the_translation_up_to_the_first_stop(fake_peptides):
    reporter = PhysicochemicalProperties()

    reporter.generate_content(_analyser("MKS*W", "MRS*W"), [{"position": 5, "alt": "G"}])

    assert set(fake_peptides.seen) == {"MKS", "MRS"}
    assert fake_peptides.seen_scale == "Bjellqvist"


def test_mutation_stopping_at_the_first_codon_reports_no_peptide(fake_peptides):
    reporter = PhysicochemicalProperties()

    content = reporter.generate_content(_analyser("MKS*", "*KS*"), [{"position": 1, "alt": "T"}])

    assert content == "Mutated sequence stops at the first codon, no peptide is produced."


def test_no_mutations_reports_nothing_to_compare(fake_peptides):
    reporter = PhysicochemicalProperties()

    assert reporter.generate_content(_analyser("MKS*", "MKS*"), []) == "No mutations applied."


def test_property_undefined_for_a_short_peptide_is_reported_as_not_available(fake_peptides):
    reporter = PhysicochemicalProperties()

    content = reporter.generate_content(_analyser("MK*", "MR*"), [{"position": 5, "alt": "G"}])

    assert "Instability Index Original Sequence: n/a" in content
    assert "Instability Index Mutated Sequence: n/a" in content
    assert "Change: n/a" in content
    assert "Molecular Weight Original Sequence: 1.234" in content
