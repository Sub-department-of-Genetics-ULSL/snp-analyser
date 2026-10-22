import pytest

from analyser_cli.protein_alignment import GAP, align_translations, group_residues


class TestAlignTranslations:
    def test_equal_translations_have_no_differences(self):
        alignment = align_translations("MKWA*", "MKWA*")

        assert alignment.changed_residues == 0
        assert alignment.identical_residues == 5
        # A stop codon is not an amino acid.
        assert alignment.original_length == alignment.mutated_length == 4

    def test_substitution_is_marked_in_place(self):
        alignment = align_translations("MKWA*", "MKWG*")

        assert alignment.original == "MKWA*"
        assert alignment.mutated == "MKWG*"
        assert alignment.differences == (False, False, False, True, False)

    def test_deleted_residue_becomes_a_gap(self):
        alignment = align_translations("MKWA*", "MKA*")

        assert alignment.original == "MKWA*"
        assert alignment.mutated == "MK-A*"
        assert alignment.changed_residues == 1
        assert alignment.original_length == 4
        assert alignment.mutated_length == 3
        assert alignment.is_homologous is True

    def test_a_gap_inside_a_repeat_is_pushed_to_the_right(self):
        alignment = align_translations("MKQQQQAWT*", "MKQQQAWT*")

        assert alignment.mutated == "MKQQQ-AWT*"

    def test_an_inserted_repeat_residue_is_pushed_to_the_right(self):
        alignment = align_translations("MKQQQAWT*", "MKQQQQAWT*")

        assert alignment.original == "MKQQQ-AWT*"

    def test_inserted_residue_gaps_the_original(self):
        alignment = align_translations("MKWA*", "MKWQA*")

        assert alignment.original == "MKW-A*"
        assert alignment.mutated == "MKWQA*"
        assert alignment.original_length == 4
        assert alignment.mutated_length == 5

    def test_truncated_translation_is_padded(self):
        alignment = align_translations("MEFILSLIGS*", "MEFI*")

        assert len(alignment.original) == len(alignment.mutated)
        assert alignment.mutated.startswith("MEFI")
        assert alignment.mutated.count("-") == 6

    def test_peptide_length_ends_at_a_premature_stop_codon(self):
        # The whole translation stays in the alignment, but translation ends at the
        # first stop codon, so the peptide is only MEFI.
        alignment = align_translations("MEFILSLIGS*", "MEFI*LSKIGS*")

        assert alignment.mutated.replace("-", "") == "MEFI*LSKIGS*"
        assert alignment.original_length == 10
        assert alignment.mutated_length == 4

    def test_peptide_length_without_stop_codon_counts_every_residue(self):
        assert align_translations("MKW*", "MKWQRS").mutated_length == 6

    def test_missing_mutated_translation_is_all_gaps(self):
        alignment = align_translations("MKWA*", "")

        assert alignment.mutated == "-----"
        assert alignment.has_mutated is False

    def test_missing_original_translation_is_all_gaps(self):
        alignment = align_translations("", "MKWA*")

        assert alignment.original == "-----"
        assert alignment.has_mutated is True

    def test_two_empty_translations(self):
        alignment = align_translations("", "")

        assert alignment.original == ""
        assert alignment.mutated == ""
        assert alignment.has_mutated is False


class TestNonHomologousTails:
    def test_a_frameshift_tail_is_not_aligned(self):
        alignment = align_translations(
            "MKTAYIAKQRQISFVKSHFSRQ*", "MKTAYIAKQRQCLRSSVILFPW*", frameshift_residue=12
        )

        assert alignment.original == "MKTAYIAKQRQISFVKSHFSRQ*"
        assert alignment.mutated == "MKTAYIAKQRQCLRSSVILFPW*"
        assert alignment.divergence_start == 12
        assert alignment.divergence_reason == "frameshift"
        assert alignment.is_homologous is False

    def test_residues_before_the_frameshift_are_still_aligned(self):
        # The substitution at residue 4 happens before the frameshift, so it has to stay a
        # substitution instead of being swallowed by the divergent tail. The codon of
        # residue 11 is the first one read in another frame, but it still encodes Ile,
        # so the frameshift starts at residue 12, the first amino acid it changes (HGVS).
        alignment = align_translations(
            "MKWACDEFGHIKLMNP*", "MKWQCDEFGHIRSTVL*", frameshift_residue=11
        )

        assert alignment.divergence_start == 12
        assert alignment.differences[:11] == (
            False, False, False, True, False, False, False, False, False, False, False,
        )
        assert alignment.identical_residues == 10

    def test_an_in_frame_indel_before_the_frameshift_moves_the_boundary(self):
        # The mutated head is one residue shorter, so the tails must not be split by index.
        alignment = align_translations(
            "MKWACDEFGHIKLMNPQRST*", "MKWCDEFGHIRSTVLPWQ*", frameshift_residue=11
        )

        assert alignment.original == "MKWACDEFGHIKLMNPQRST*"
        assert alignment.mutated == "MKW-CDEFGHIRSTVLPWQ*-"
        assert alignment.divergence_start == 12

    def test_a_frameshift_never_opens_a_gap_inside_the_tail(self):
        alignment = align_translations("MKWAGGGCDEF*", "MKWQTTT*", frameshift_residue=4)

        assert alignment.original == "MKWAGGGCDEF*"
        assert alignment.mutated == "MKWQTTT*----"
        assert GAP not in alignment.mutated.rstrip(GAP)

    def test_a_frameshift_tail_is_reported_as_changed_even_where_it_matches(self):
        alignment = align_translations("MKWACDE*", "MKWQCDE*", frameshift_residue=4)

        # The shared 'CDE*' is read in another frame, so it is a coincidence not a match.
        assert alignment.differences == (False, False, False, True, True, True, True, True)
        assert alignment.changed_residues == 5

    def test_a_frameshift_past_the_stop_codon_is_ignored(self):
        alignment = align_translations("MKWA*", "MKWG*", frameshift_residue=99)

        assert alignment.is_homologous is True
        assert alignment.changed_residues == 1

    def test_a_frameshift_on_the_first_residue_leaves_no_homologous_head(self):
        alignment = align_translations("MKWA*", "QLPT*", frameshift_residue=1)

        assert alignment.divergence_start == 1
        assert alignment.changed_residues == 5

    def test_a_premature_stop_codon_truncates_instead_of_deleting(self):
        alignment = align_translations("MKTAYIAKQRQISFVK*", "MKTAYIAK*")

        assert alignment.original == "MKTAYIAKQRQISFVK*"
        assert alignment.mutated == "MKTAYIAK*--------"
        assert alignment.divergence_start == 9
        assert alignment.divergence_reason == "truncation"

    def test_an_in_frame_deletion_is_still_aligned(self):
        alignment = align_translations("MKTAYIAKQRQISFVK*", "MKTAYIAKQRISFVK*")

        assert alignment.is_homologous is True
        assert alignment.changed_residues == 1

    def test_a_frameshift_flag_on_an_unchanged_translation_changes_nothing(self):
        alignment = align_translations("MKWA*", "MKWA*", frameshift_residue=4)

        assert alignment.is_homologous is True
        assert alignment.changed_residues == 0


class TestAlignmentLines:
    def test_lines_are_wrapped_and_numbered(self):
        alignment = align_translations("M" * 130, "M" * 130)
        lines = alignment.lines(line_length=60, group_size=5)

        assert [line.original_start for line in lines] == [1, 61, 121]
        assert [line.original_end for line in lines] == [60, 120, 130]
        assert [len(line.original) for line in lines] == [60, 60, 10]

    def test_line_length_is_rounded_down_to_whole_groups(self):
        alignment = align_translations("M" * 130, "M" * 130)
        lines = alignment.lines(line_length=58, group_size=5)

        assert [len(line.original) for line in lines] == [55, 55, 20]

    def test_positions_skip_gaps(self):
        alignment = align_translations("MKWA*", "MKA*")
        line = alignment.lines(line_length=60, group_size=5)[0]

        assert (line.original_start, line.original_end) == (1, 5)
        assert (line.mutated_start, line.mutated_end) == (1, 4)

    def test_a_line_made_only_of_gaps_has_no_positions(self):
        alignment = align_translations("M" * 130, "")
        lines = alignment.lines(line_length=60, group_size=5)

        assert lines[1].mutated_start is None
        assert lines[1].mutated_end is None

    def test_the_divergence_is_reported_on_the_line_that_holds_it(self):
        alignment = align_translations("M" * 130, "M" * 70 + "K" * 60, frameshift_residue=71)
        lines = alignment.lines(line_length=60, group_size=5)

        assert alignment.divergence_start == 71
        assert [line.divergence for line in lines] == [None, 10, None]

    def test_a_homologous_alignment_has_no_divergence(self):
        alignment = align_translations("M" * 130, "M" * 130)
        lines = alignment.lines(line_length=60, group_size=5)

        assert [line.divergence for line in lines] == [None, None, None]

    def test_invalid_line_length_raises(self):
        with pytest.raises(ValueError):
            align_translations("MKWA*", "MKWA*").lines(line_length=0)


class TestGroupResidues:
    def test_residues_are_grouped_by_five(self):
        assert group_residues("MKWACDEFGHIK", 5) == "MKWAC DEFGH IK"

    def test_group_size_below_one_keeps_the_string(self):
        assert group_residues("MKWAC", 0) == "MKWAC"
