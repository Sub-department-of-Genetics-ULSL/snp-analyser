"""Every cached GenBank CDS has to translate exactly like its /translation qualifier.

The GenBank flat files are downloaded on demand and are not committed (``*.gb`` is
ignored), so the tests are skipped when no genome has been cached yet.
"""

from pathlib import Path

import pytest
from Bio import SeqIO

from analyser_cli.analyser import Analyser
from analyser_cli.cds_annotation import read_cds_annotation

GENOME_FOLDERS = [
    Path(__file__).resolve().parents[2] / "analyser_backend" / "genomic_data",
    Path(__file__).resolve().parents[2] / "genomic_data",
]


def _cached_genomes() -> list[Path]:
    for folder in GENOME_FOLDERS:
        genomes = sorted(folder.glob("*.gb"))
        if genomes:
            return genomes
    return []


def _coding_features(genome_path: Path):
    record = SeqIO.read(genome_path, "gb")
    for feature in record.features:
        if feature.type == "CDS" and feature.qualifiers.get("translation"):
            yield record, feature


@pytest.mark.parametrize(
    "genome_path", _cached_genomes(), ids=lambda path: path.stem
)
def test_translation_matches_the_genbank_qualifier(genome_path):
    checked = 0

    for record, feature in _coding_features(genome_path):
        annotation = read_cds_annotation(feature)
        analyser = Analyser(
            str(feature.extract(record.seq)),
            transl_table=annotation.transl_table,
            codon_start=annotation.codon_start,
            transl_except=annotation.transl_except,
        )

        gene = feature.qualifiers.get("gene", ["?"])[0]
        # GenBank never stores the terminating stop codon in /translation.
        assert str(analyser.amino_acid_translation).rstrip("*") == annotation.translation, (
            f"{genome_path.name} {gene} (transl_table={annotation.transl_table})"
        )
        checked += 1

    assert checked > 0


@pytest.mark.parametrize(
    "genome_path", _cached_genomes(), ids=lambda path: path.stem
)
def test_every_coding_sequence_declares_a_genetic_code(genome_path):
    for record, feature in _coding_features(genome_path):
        annotation = read_cds_annotation(feature)
        assert annotation.transl_table in (2, 5), (
            f"{genome_path.name} uses an unexpected genetic code {annotation.transl_table}"
        )


@pytest.mark.parametrize(
    "genome_path", _cached_genomes(), ids=lambda path: path.stem
)
def test_hgvs_protein_reference_matches_the_analyser(genome_path):
    """p. descriptions are numbered on the same protein the report translates."""
    from analyser_cli.hgvs_protein import translate_protein

    for record, feature in _coding_features(genome_path):
        annotation = read_cds_annotation(feature)
        sequence = str(feature.extract(record.seq))
        analyser = Analyser(
            sequence,
            transl_table=annotation.transl_table,
            codon_start=annotation.codon_start,
            transl_except=annotation.transl_except,
        )
        protein, terminated = translate_protein(
            sequence,
            analyser.codon_table,
            codon_start=annotation.codon_start,
            exceptions={
                (exception.start - analyser.coding_offset) // 3: exception.amino_acid
                for exception in annotation.transl_except
            },
        )

        gene = feature.qualifiers.get("gene", ["?"])[0]
        assert terminated, f"{genome_path.name} {gene}"
        assert protein == str(analyser.amino_acid_translation), f"{genome_path.name} {gene}"


@pytest.mark.parametrize(
    "genome_path", _cached_genomes(), ids=lambda path: path.stem
)
def test_annotated_partial_stop_codons_terminate_the_translation(genome_path):
    for record, feature in _coding_features(genome_path):
        annotation = read_cds_annotation(feature)
        terminal = [
            exception
            for exception in annotation.transl_except
            if exception.amino_acid == "*" and exception.end >= len(feature)
        ]
        if not terminal:
            continue

        analyser = Analyser(
            str(feature.extract(record.seq)),
            transl_table=annotation.transl_table,
            codon_start=annotation.codon_start,
            transl_except=annotation.transl_except,
        )

        gene = feature.qualifiers.get("gene", ["?"])[0]
        assert str(analyser.amino_acid_translation).endswith("*"), f"{genome_path.name} {gene}"
