import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from analyser_cli.analyser import Analyser, get_codon_table
from analyser_cli.cds_annotation import (
    DEFAULT_CODON_START,
    DEFAULT_TRANSL_TABLE,
    TranslationException,
    read_cds_annotation,
)
from analyser_cli.data_manager import DataManager
from analyser_cli.sub_reporters.mutation_effects import MutationEffectsReport


# Codon 2 (AAA) is changed to AGA, which terminates translation in vertebrate
# mitochondria (table 2) but encodes serine in invertebrate mitochondria (table 5).
SEQUENCE = "ATGAAACGATAA"
MUTATIONS = [{"type": "sub", "position": 5, "alt": "G"}]


class TestGetCodonTable:
    @pytest.mark.parametrize("table_id", [1, 2, 3, 4, 5, 6, 9, 11, 24, 33])
    def test_resolves_every_ncbi_id(self, table_id):
        assert get_codon_table(table_id).id == table_id

    def test_accepts_the_id_as_text(self):
        assert get_codon_table("5").id == 5

    def test_defaults_to_the_standard_code(self):
        assert get_codon_table().id == DEFAULT_TRANSL_TABLE

    @pytest.mark.parametrize("table_id", [0, 7, 8, 99, "vertebrate", None])
    def test_unknown_id_raises(self, table_id):
        with pytest.raises(ValueError):
            get_codon_table(table_id)


class TestAnalyserGeneticCode:
    def test_defaults_to_the_standard_code(self):
        analyser = Analyser(SEQUENCE)
        assert analyser.transl_table == 1
        assert analyser.codon_table.id == 1

    def test_aga_terminates_translation_in_vertebrate_mitochondria(self):
        analyser = Analyser(SEQUENCE, transl_table=2)
        analyser.apply_mutations(MUTATIONS)
        assert str(analyser.mutated_amino_acid_translation) == "M*R*"

    def test_aga_encodes_serine_in_invertebrate_mitochondria(self):
        analyser = Analyser(SEQUENCE, transl_table=5)
        analyser.apply_mutations(MUTATIONS)
        assert str(analyser.mutated_amino_acid_translation) == "MSR*"

    def test_unknown_table_raises(self):
        with pytest.raises(ValueError):
            Analyser(SEQUENCE, transl_table=99)

    def test_codon_table_name_is_reported(self):
        assert Analyser(SEQUENCE, transl_table=5).codon_table_name == "Invertebrate Mitochondrial"

    def test_genetic_code_description_carries_the_ncbi_id(self):
        assert Analyser(SEQUENCE, transl_table=2).genetic_code_description == (
            "Vertebrate Mitochondrial (NCBI transl_table=2)"
        )


class TestAlternativeStartCodons:
    @pytest.mark.parametrize("start_codon", ["ATT", "ATA", "GTG", "TTG"])
    def test_alternative_start_codon_translates_to_methionine(self, start_codon):
        analyser = Analyser(f"{start_codon}CGATAA", transl_table=5)
        assert str(analyser.amino_acid_translation).startswith("M")

    def test_non_start_codon_keeps_its_own_amino_acid(self):
        analyser = Analyser("CGACGATAA", transl_table=2)
        assert str(analyser.amino_acid_translation) == "RR*"

    def test_mutation_destroying_the_start_codon_is_visible(self):
        analyser = Analyser("ATGCGATAA", transl_table=2)
        analyser.apply_mutations([{"type": "sub", "position": 2, "alt": "C"}])

        assert str(analyser.amino_acid_translation) == "MR*"
        assert str(analyser.mutated_amino_acid_translation) == "TR*"


class TestCodonStart:
    def test_codon_start_shifts_the_reading_frame(self):
        # The leading GG belongs to a partial codon that is not translated.
        analyser = Analyser("GGATGCGATAA", transl_table=2, codon_start=3)
        assert analyser.coding_offset == 2
        assert str(analyser.amino_acid_translation) == "MR*"

    def test_codon_start_one_translates_from_the_first_base(self):
        analyser = Analyser("ATGCGATAA", transl_table=2, codon_start=DEFAULT_CODON_START)
        assert analyser.coding_offset == 0
        assert str(analyser.amino_acid_translation) == "MR*"


class TestTranslationExceptions:
    def test_partial_terminal_codon_becomes_a_stop(self):
        # "ATGCGAT" ends with a single T completed to TAA by polyadenylation.
        analyser = Analyser(
            "ATGCGAT",
            transl_table=2,
            transl_except=(TranslationException(6, 7, "*"),),
        )
        assert str(analyser.amino_acid_translation) == "MR*"

    def test_internal_exception_overrides_the_codon(self):
        analyser = Analyser(
            "ATGCGATAA",
            transl_table=2,
            transl_except=(TranslationException(3, 6, "U"),),
        )
        assert str(analyser.amino_acid_translation) == "MU*"

    def test_terminal_exception_follows_the_end_of_a_shortened_sequence(self):
        # "ATGCGACGAT" ends with a single T completed to TAA by polyadenylation.
        analyser = Analyser(
            "ATGCGACGAT",
            transl_table=2,
            transl_except=(TranslationException(9, 10, "*"),),
        )
        assert str(analyser.amino_acid_translation) == "MRR*"

        # An in-frame deletion keeps the lone T at the 3' end, so the stop codon is still
        # completed by polyadenylation - three bases further to the left.
        analyser.apply_mutations([{"type": "del", "position": position} for position in (4, 5, 6)])

        assert str(analyser.mutated_sequence) == "ATGCGAT"
        assert str(analyser.mutated_amino_acid_translation) == "MR*"

    def test_a_partial_stop_codon_is_lost_when_its_bases_change(self):
        # T-- completed to TAA is a stop codon, C-- completed to CAA encodes glutamine.
        analyser = Analyser(
            "ATGCGACGAT",
            transl_table=2,
            transl_except=(TranslationException(9, 10, "*"),),
        )
        analyser.apply_mutations([{"type": "sub", "position": 10, "alt": "C"}])

        assert str(analyser.mutated_amino_acid_translation) == "MRRQ"

    def test_a_shifted_partial_codon_is_completed_by_the_poly_a_tail(self):
        analyser = Analyser(
            "ATGCGACGATA",
            transl_table=2,
            transl_except=(TranslationException(9, 11, "*"),),
        )
        assert str(analyser.amino_acid_translation) == "MRR*"

        # Deleting one base leaves a lone A, completed to AAA (lysine), not to a stop codon.
        analyser.apply_mutations([{"type": "del", "position": 4}])

        assert str(analyser.mutated_sequence) == "ATGGACGATA"
        assert str(analyser.mutated_amino_acid_translation) == "MDDK"

    def test_terminal_exception_is_dropped_when_the_last_codon_becomes_complete(self):
        analyser = Analyser(
            "ATGCGACGAT",
            transl_table=2,
            transl_except=(TranslationException(9, 10, "*"),),
        )
        assert str(analyser.amino_acid_translation) == "MRR*"

        analyser.apply_mutations([{"type": "del", "position": 4}])

        assert str(analyser.mutated_sequence) == "ATGGACGAT"
        assert str(analyser.mutated_amino_acid_translation) == "MDD"

    def test_internal_exception_is_dropped_after_an_indel(self):
        analyser = Analyser(
            "ATGCGATAA",
            transl_table=2,
            transl_except=(TranslationException(3, 6, "U"),),
        )
        analyser.apply_mutations([{"type": "ins", "position": 1, "alt": "A"}])

        assert "U" not in str(analyser.mutated_amino_acid_translation)

    def test_the_serialisable_dict_form_is_accepted(self):
        analyser = Analyser(
            "ATGCGAT",
            transl_table=2,
            transl_except=[{"start": 6, "end": 7, "aminoAcid": "*"}],
        )
        assert analyser.transl_except == (TranslationException(6, 7, "*"),)
        assert str(analyser.amino_acid_translation) == "MR*"


class TestReadCdsAnnotation:
    def test_reads_every_qualifier(self):
        feature = SeqFeature(
            FeatureLocation(0, 9),
            type="CDS",
            qualifiers={
                "gene": ["ND1"],
                "transl_table": ["2"],
                "codon_start": ["1"],
                "translation": ["MR"],
            },
        )
        annotation = read_cds_annotation(feature)

        assert annotation.transl_table == 2
        assert annotation.codon_start == 1
        assert annotation.translation == "MR"
        assert annotation.transl_except == ()

    def test_missing_qualifiers_fall_back_to_the_standard_code(self):
        feature = SeqFeature(FeatureLocation(0, 9), type="gene", qualifiers={"gene": ["ND1"]})
        annotation = read_cds_annotation(feature)

        assert annotation.transl_table == DEFAULT_TRANSL_TABLE
        assert annotation.codon_start == DEFAULT_CODON_START
        assert annotation.translation is None

    def test_single_base_exception_is_mapped_to_a_gene_offset(self):
        feature = SeqFeature(
            FeatureLocation(10, 20),
            type="CDS",
            qualifiers={"gene": ["ND1"], "transl_except": ["(pos:20,aa:TERM)"]},
        )
        assert read_cds_annotation(feature).transl_except == (TranslationException(9, 10, "*"),)

    def test_two_base_exception_is_mapped_to_a_gene_range(self):
        feature = SeqFeature(
            FeatureLocation(10, 20),
            type="CDS",
            qualifiers={"gene": ["ND1"], "transl_except": ["(pos:19..20,aa:TERM)"]},
        )
        assert read_cds_annotation(feature).transl_except == (TranslationException(8, 10, "*"),)

    def test_exception_on_the_complement_strand_is_mapped_in_transcription_order(self):
        feature = SeqFeature(
            FeatureLocation(10, 20, strand=-1),
            type="CDS",
            qualifiers={"gene": ["ND1"], "transl_except": ["(pos:complement(11),aa:TERM)"]},
        )
        # Genome position 11 is the last base that is transcribed on the minus strand.
        assert read_cds_annotation(feature).transl_except == (TranslationException(9, 10, "*"),)

    def test_exception_outside_the_feature_is_ignored(self):
        feature = SeqFeature(
            FeatureLocation(10, 20),
            type="CDS",
            qualifiers={"gene": ["ND1"], "transl_except": ["(pos:5,aa:TERM)"]},
        )
        assert read_cds_annotation(feature).transl_except == ()


class TestDataManagerReadsTheCdsFeature:
    @staticmethod
    def _genome():
        genome = SeqRecord(Seq("ATGCGATAAGGG"), id="test")
        genome.features = [
            SeqFeature(FeatureLocation(0, 9), type="gene", qualifiers={"gene": ["ND1"]}),
            SeqFeature(
                FeatureLocation(0, 9),
                type="CDS",
                qualifiers={
                    "gene": ["ND1"],
                    "transl_table": ["5"],
                    "codon_start": ["1"],
                    "transl_except": ["(pos:9,aa:TERM)"],
                    "translation": ["MR"],
                },
            ),
        ]
        return genome

    def test_cds_qualifiers_win_over_the_gene_feature(self, monkeypatch):
        manager = DataManager()
        monkeypatch.setattr(manager, "get_genome", lambda latin_name: self._genome())

        gene_data = manager.get_gene_for_organism("homo_sapiens", "ND1")

        assert gene_data["translTable"] == 5
        assert gene_data["codonStart"] == 1
        assert gene_data["translExcept"] == [{"start": 8, "end": 9, "aminoAcid": "*"}]
        assert gene_data["ncbiTranslation"] == "MR"

    def test_features_without_a_cds_use_the_standard_code(self, monkeypatch):
        genome = SeqRecord(Seq("ATGCGATAA"), id="test")
        genome.features = [
            SeqFeature(FeatureLocation(0, 9), type="tRNA", qualifiers={"gene": ["trnF"]})
        ]
        manager = DataManager()
        monkeypatch.setattr(manager, "get_genome", lambda latin_name: genome)

        gene_data = manager.get_gene_for_organism("homo_sapiens", "trnF")

        assert gene_data["translTable"] == DEFAULT_TRANSL_TABLE
        assert gene_data["translExcept"] == []

    def test_unknown_gene_raises(self, monkeypatch):
        manager = DataManager()
        monkeypatch.setattr(manager, "get_genome", lambda latin_name: self._genome())

        with pytest.raises(ValueError):
            manager.get_gene_for_organism("homo_sapiens", "COX1")


class TestMutationEffectsUsesCodonTable:
    def test_effect_type_depends_on_the_genetic_code(self):
        vertebrate = Analyser(SEQUENCE, transl_table=2)
        vertebrate.apply_mutations(MUTATIONS)
        invertebrate = Analyser(SEQUENCE, transl_table=5)
        invertebrate.apply_mutations(MUTATIONS)

        report = MutationEffectsReport()
        vertebrate_content = report.generate_content(vertebrate, MUTATIONS)
        invertebrate_content = report.generate_content(invertebrate, MUTATIONS)
        assert "Effect Type: Nonsense (premature stop codon)" in vertebrate_content
        assert "HGVS (protein): p.(Lys2Ter)" in vertebrate_content
        assert "Effect Type: Missense" in invertebrate_content
        assert "HGVS (protein): p.(Lys2Ser)" in invertebrate_content

    def test_genetic_code_is_included_in_report(self):
        analyser = Analyser(SEQUENCE, transl_table=5)
        analyser.apply_mutations(MUTATIONS)
        content = MutationEffectsReport().generate_content(analyser, MUTATIONS)
        assert "Genetic Code: Invertebrate Mitochondrial (NCBI transl_table=5)" in content

    def test_codon_numbering_follows_codon_start(self):
        analyser = Analyser("GG" + SEQUENCE, transl_table=2, codon_start=3)
        mutations = [{"type": "sub", "position": 7, "alt": "G"}]
        analyser.apply_mutations(mutations)

        content = MutationEffectsReport().generate_content(analyser, mutations)

        assert "Codon Position: Codon 2, Position 2" in content

    def test_peptide_length_ends_at_a_premature_stop_codon(self):
        # CGA > AGA is a stop codon in the vertebrate mitochondrial code.
        analyser = Analyser("ATGCGATCGTAA", transl_table=2)
        mutations = [{"type": "sub", "position": 4, "alt": "A"}]
        analyser.apply_mutations(mutations)

        content = MutationEffectsReport().generate_content(analyser, mutations)

        assert "Mutated: M*S*" in content
        assert "Peptide Length: 3 -> 1 amino acid (up to the first stop codon)" in content
