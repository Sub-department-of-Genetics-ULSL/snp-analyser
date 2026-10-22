"""Predicted protein consequences, checked against the examples of the HGVS recommendations."""

import pytest

from analyser_cli.cds_annotation import TranslationException, get_codon_table
from analyser_cli.hgvs_protein import (
    EXTENSION,
    FRAMESHIFT,
    INFRAME_DELETION,
    INFRAME_DELINS,
    INFRAME_DUPLICATION,
    INFRAME_INSERTION,
    MISSENSE,
    NONSENSE,
    START,
    SYNONYMOUS,
    CdsEdit,
    apply_edits,
    predict_consequence,
    three_letter,
    translate_protein,
)

STANDARD = get_codon_table(1)
VERTEBRATE_MITOCHONDRIAL = get_codon_table(2)


def predict(cds, *edits, table=STANDARD, **options):
    return predict_consequence(cds, list(edits), table, **options)


def substitution(position, base):
    return CdsEdit(position - 1, position, base)


class TestTranslation:
    def test_alternative_start_codon_initiates_with_methionine(self):
        assert translate_protein("ATTAAATAA", VERTEBRATE_MITOCHONDRIAL) == ("MK*", True)

    def test_partial_stop_codon_is_completed_by_the_poly_a_tail(self):
        assert translate_protein("ATGAAAT", VERTEBRATE_MITOCHONDRIAL) == ("MK*", True)
        assert translate_protein("ATGAAATA", VERTEBRATE_MITOCHONDRIAL) == ("MK*", True)

    def test_a_frame_without_stop_codon_is_open(self):
        assert translate_protein("ATGAAAC", VERTEBRATE_MITOCHONDRIAL) == ("MKQ", False)

    def test_edits_are_applied_in_reference_coordinates(self):
        assert apply_edits("ATGCC", [CdsEdit(2, 4, "TT"), CdsEdit(0, 1, "G")]) == "GTTTC"
        assert apply_edits("ATGCC", [CdsEdit(2, 2, "A"), CdsEdit(2, 3, "")]) == "ATACC"


class TestSubstitutions:
    def test_missense(self):
        consequence = predict("ATGTGGTAA", substitution(6, "T"))  # TGG -> TGT, Trp2 -> Cys
        assert (consequence.edit, consequence.kind) == ("Trp2Cys", MISSENSE)
        assert consequence.description("NP_003997.1") == "NP_003997.1:p.(Trp2Cys)"

    def test_nonsense(self):
        consequence = predict("ATGTGGTAA", substitution(6, "A"))  # TGG -> TGA
        assert (consequence.edit, consequence.kind) == ("Trp2Ter", NONSENSE)

    def test_silent_change_names_the_unchanged_residue(self):
        # HGVS: p.(Leu54=), not p.Leu54Leu.
        consequence = predict("ATGCTGTAA", substitution(6, "A"))
        assert (consequence.edit, consequence.kind) == ("Leu2=", SYNONYMOUS)

    @pytest.mark.parametrize("position, base", [(1, "G"), (2, "C"), (3, "A")])
    def test_a_changed_start_codon_has_an_unknown_consequence(self, position, base):
        # HGVS: p.(Met1?) - p.Met1Thr is never a correct description.
        consequence = predict("ATGTGGTAA", substitution(position, base))
        assert (consequence.edit, consequence.kind) == ("Met1?", START)

    def test_an_alternative_start_codon_is_methionine_one(self):
        # ATT initiates with Met in the vertebrate mitochondrial code.
        consequence = predict("ATTTGGTAA", substitution(3, "C"), table=VERTEBRATE_MITOCHONDRIAL)
        assert consequence.edit == "Met1?"

    def test_an_incomplete_cds_has_no_start_codon_rule(self):
        consequence = predict("TCGTGGTAA", substitution(1, "G"), complete_start=False)
        assert consequence.edit == "Ser1Ala"

    def test_stop_codon_change_keeping_a_stop_is_silent(self):
        consequence = predict("ATGTGGTAA", substitution(9, "G"))  # TAA -> TAG
        assert consequence.edit == "Ter3="


class TestExtensions:
    def test_lost_stop_codon_reading_into_a_new_one(self):
        # HGVS: p.Ter110GlnextTer17 - the count is the number of amino acids added.
        cds = "ATGTGGTAA"
        consequence = predict(cds + "GCTTGA", substitution(7, "C"))
        assert (consequence.edit, consequence.kind) == ("Ter3GlnextTer2", EXTENSION)

    def test_lost_stop_codon_reading_into_the_poly_a_tail(self):
        consequence = predict("ATGTGGTAA", substitution(7, "C"))
        assert consequence.edit == "Ter3GlnextTer?"

    def test_poly_a_completed_stop_codon_can_be_lost(self):
        # T-- is completed to TAA; C-- to CAA (Gln) in the vertebrate mitochondrial code.
        consequence = predict(
            "ATGTGGT",
            substitution(7, "C"),
            table=VERTEBRATE_MITOCHONDRIAL,
            transl_except=[TranslationException(6, 7, "*")],
        )
        assert consequence.edit == "Ter3GlnextTer?"

    def test_insertion_right_before_the_stop_codon_is_an_extension(self):
        consequence = predict("ATGTGGTAA", CdsEdit(6, 6, "GCT"))
        assert consequence.edit == "Ter3AlaextTer1"


class TestFrameshifts:
    def test_hgvs_example_starts_at_the_first_changed_amino_acid(self):
        # HGVS: c.6_13dup in ATGGATGCATACGAGATGAGG is p.Glu5ValfsTer5.
        cds = "ATGGATGCATACGAGATGAGGTAA"
        consequence = predict(cds, CdsEdit(13, 13, cds[5:13]))
        assert (consequence.edit, consequence.kind) == ("Glu5ValfsTer5", FRAMESHIFT)
        assert consequence.first_residue == 5

    def test_a_deletion_creating_an_immediate_stop_is_nonsense(self):
        # HGVS: c.12del in ATGGATGCATACGTCACG gives p.(Tyr4*), not p.Tyr4TerfsTer1.
        consequence = predict("ATGGATGCATACGTCACGTAA", CdsEdit(11, 12, ""))
        assert (consequence.edit, consequence.kind) == ("Tyr4Ter", NONSENSE)

    def test_a_frame_without_new_stop_codon(self):
        consequence = predict("ATGAAAAAAAAATAA", CdsEdit(3, 4, ""))
        assert consequence.edit.endswith("fsTer?")

    def test_a_deletion_including_the_stop_codon_is_a_frameshift(self):
        # HGVS: deletions starting 5' of and including the stop codon are frameshifts.
        consequence = predict("ATGTGGGCTTAA", CdsEdit(6, 12, ""))
        assert consequence.kind == FRAMESHIFT
        assert consequence.edit.startswith("Ala3")

    def test_a_compensated_frameshift_is_a_deletion_insertion(self):
        # M K W H F *: deleting c.4 and inserting C after c.12 shifts only codons 2-4.
        cds = "ATGAAATGGCATTTTTAA"
        consequence = predict(cds, CdsEdit(3, 4, ""), CdsEdit(12, 12, "C"))
        assert consequence.kind == INFRAME_DELINS
        assert consequence.edit == "Lys2_His4delinsAsnGlyIle"

    def test_shifted_codons_encoding_the_same_residues_are_not_reported(self):
        cds = "ATGAAACCCGGGTTTTAA"
        consequence = predict(cds, CdsEdit(3, 4, ""), CdsEdit(9, 9, "C"))
        assert consequence.edit == "Lys2Asn"


class TestInFrameChanges:
    def test_deletion_in_a_stretch_takes_the_c_terminal_residue(self):
        # HGVS: MetLeuTrpTrpGlu -> MetLeuTrpGlu is p.Trp4del.
        consequence = predict("ATGCTGTGGTGGGAATAA", CdsEdit(6, 9, ""))
        assert (consequence.edit, consequence.kind) == ("Trp4del", INFRAME_DELETION)

    def test_duplication_in_a_stretch_takes_the_c_terminal_residue(self):
        # HGVS: MetLeuTrpTrpGlu -> MetLeuTrpTrpTrpGlu is p.Trp4dup.
        consequence = predict("ATGCTGTGGTGGGAATAA", CdsEdit(6, 6, "TGG"))
        assert (consequence.edit, consequence.kind) == ("Trp4dup", INFRAME_DUPLICATION)

    def test_insertion_with_a_stop_codon_is_an_insertion(self):
        # HGVS: c.9_10insGGGTAG gives p.(Met3_His4insGlyTer).
        consequence = predict("ATGAAAATGCATCAGTAA", CdsEdit(9, 9, "GGGTAG"))
        assert (consequence.edit, consequence.kind) == ("Met3_His4insGlyTer", INFRAME_INSERTION)

    def test_insertion(self):
        # HGVS: MetLysGlyHisGlnGlnCys -> MetLysGlyHisAlaGlnGlnCys is p.His4_Gln5insAla.
        cds = "ATGAAAGGCCATCAGCAGTGCTAA"
        consequence = predict(cds, CdsEdit(12, 12, "GCT"))
        assert consequence.edit == "His4_Gln5insAla"

    def test_delins_of_two_residues(self):
        consequence = predict("ATGAAATGGTAA", CdsEdit(3, 9, "TCTTGT"))
        assert (consequence.edit, consequence.kind) == ("Lys2_Trp3delinsSerCys", INFRAME_DELINS)

    def test_a_merged_codon_change_is_a_substitution(self):
        # HGVS: c.235_237delinsTAT is p.Lys79Tyr.
        consequence = predict("ATGAAGTAA", CdsEdit(3, 6, "TAT"))
        assert consequence.edit == "Lys2Tyr"

    def test_duplication_of_the_stop_codon_changes_nothing(self):
        consequence = predict("ATGTGGTAA", CdsEdit(9, 9, "TAA"))
        assert (consequence.edit, consequence.kind) == ("=", SYNONYMOUS)


def test_three_letter_code_uses_ter_for_stop():
    assert three_letter("MW*") == "MetTrpTer"
    assert three_letter("UOX") == "SecPylXaa"
