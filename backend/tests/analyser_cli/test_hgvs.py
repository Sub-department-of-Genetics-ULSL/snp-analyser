"""HGVS (version 21.1) descriptions of variants in one coding sequence of a genome."""

from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from analyser_cli.hgvs import (
    HGVS_VERSION,
    HgvsError,
    ReferenceContext,
    describe_variants,
    reverse_complement,
)

# Codons M K Q F G * - an A run (c.4_6), a CAG unit (c.7_9) and TTC (c.10_12).
CDS = "ATGAAACAGTTCGGATAA"
UPSTREAM = "GCGCGCGCGC"
DOWNSTREAM = "CCGCCGCCGC"
GENOME_FOLDERS = [
    Path(__file__).resolve().parents[2] / "analyser_backend" / "genomic_data",
    Path(__file__).resolve().parents[2] / "genomic_data",
]


def make_context(
    genome: str,
    location,
    *,
    gene: str = "ND2",
    symbol: str | None = "MT-ND2",
    circular: bool = True,
    transl_table: int = 1,
) -> ReferenceContext:
    record = SeqRecord(
        Seq(genome),
        id="NC_000999.1",
        annotations={"topology": "circular" if circular else "linear"},
    )
    cds = SeqFeature(
        location,
        type="CDS",
        qualifiers={
            "gene": [gene],
            "transl_table": [str(transl_table)],
            "codon_start": ["1"],
            "protein_id": ["YP_000999.1"],
        },
    )
    record.features = [cds]
    if symbol:
        record.features.append(
            SeqFeature(
                location,
                type="gene",
                qualifiers={"gene": [gene], "nomenclature": [f"Official Symbol: {symbol} | Name: x"]},
            )
        )
    return ReferenceContext.from_genbank(record, cds)


@pytest.fixture
def plus():
    """The CDS on the plus strand, genome positions 11..28."""
    return make_context(UPSTREAM + CDS + DOWNSTREAM, FeatureLocation(10, 28, strand=1))


@pytest.fixture
def minus():
    """The same CDS read from complement(11..28): c.1 is genome position 28."""
    return make_context(
        UPSTREAM + reverse_complement(CDS) + DOWNSTREAM, FeatureLocation(10, 28, strand=-1)
    )


def describe(context, *descriptions):
    result = describe_variants(context, list(descriptions))
    assert not result.errors, [error.message for error in result.errors]
    return result


def single(context, description):
    result = describe(context, description)
    assert len(result.variants) == 1
    return result.variants[0]


def error_of(context, description) -> str:
    result = describe_variants(context, [description])
    assert not result.variants
    assert len(result.errors) == 1
    return result.errors[0].message


class TestReferenceContext:
    def test_plus_strand_coordinates(self, plus):
        assert plus.cds_sequence == CDS
        assert plus.strand == 1
        assert (plus.first_genome_position, plus.last_genome_position) == (11, 28)
        assert plus.symbol == "MT-ND2"
        assert plus.c_reference == "NC_000999.1(MT-ND2):c."
        assert plus.m_reference == "NC_000999.1:m."

    def test_complement_location_reverses_the_genome_numbering(self, minus):
        # complement(11..28): c.1 is the A of the start codon at genome position 28.
        assert minus.cds_sequence == CDS
        assert minus.strand == -1
        assert minus.first_genome_position == 28
        assert minus.last_genome_position == 11
        assert [minus.genome_position(index) for index in (0, 1, 2, 17)] == [28, 27, 26, 11]

    def test_a_joined_location_is_rejected(self):
        location = CompoundLocation([FeatureLocation(10, 16, 1), FeatureLocation(19, 28, 1)])
        with pytest.raises(HgvsError, match="not one contiguous stretch"):
            make_context(UPSTREAM + CDS + DOWNSTREAM, location)

    def test_a_cds_across_the_origin_of_a_circular_genome_is_contiguous(self):
        genome = CDS[6:] + "GCGCGCGCGC" + CDS[:6]
        location = CompoundLocation(
            [FeatureLocation(len(genome) - 6, len(genome), 1), FeatureLocation(0, 12, 1)]
        )
        context = make_context(genome, location)

        assert context.cds_sequence == CDS
        variant = single(context, "c.7C>T")
        assert variant.m == "NC_000999.1:m.1C>T"

    def test_symbol_falls_back_to_the_genbank_gene_name(self):
        context = make_context(UPSTREAM + CDS + DOWNSTREAM, FeatureLocation(10, 28, 1), symbol=None)
        assert context.symbol == "ND2"

    def test_a_symbol_with_hgvs_characters_is_not_used(self):
        context = make_context(
            UPSTREAM + CDS + DOWNSTREAM, FeatureLocation(10, 28, 1), symbol="mt:ND2"
        )
        assert context.symbol == "ND2"
        assert "mt:ND2" in context.aliases


class TestComplementStrand:
    def test_c_and_m_describe_the_same_change_with_complementary_bases(self, minus):
        variant = single(minus, "c.4A>G")
        # c.4 is genome position 25, where the plus strand has T.
        assert variant.m == "NC_000999.1:m.25T>C"
        assert variant.c == "NC_000999.1(MT-ND2):c.4A>G"
        assert variant.p == "YP_000999.1:p.(Lys2Glu)"

    def test_m_input_is_mapped_through_the_reversed_numbering(self, minus):
        assert single(minus, "m.25T>C").c == "NC_000999.1(MT-ND2):c.4A>G"

    def test_m_input_uses_plus_strand_bases(self, minus):
        assert "is T, not A" in error_of(minus, "m.25A>G")

    def test_the_3_prime_rule_follows_each_reference_sequence(self, minus):
        # Coding AAA (c.4_6) is TTT at m.23_25: the c. description names the 3' end of the
        # gene (c.6 = m.23), the m. description the 3' end of the genome (m.25 = c.4).
        variant = single(minus, "c.4del")
        assert variant.c == "NC_000999.1(MT-ND2):c.6del"
        assert variant.m == "NC_000999.1:m.25del"

    def test_a_duplication_is_found_on_both_strands(self, minus):
        variant = single(minus, "c.6_7insA")
        assert variant.c == "NC_000999.1(MT-ND2):c.6dup"
        assert variant.m == "NC_000999.1:m.25dup"

    def test_repeat_units_follow_the_orientation(self, minus):
        variant = single(minus, "c.9_10insCAGCAG")
        assert variant.c == "NC_000999.1(MT-ND2):c.7_9CAG[3]"
        assert variant.m == "NC_000999.1:m.20_22CTG[3]"


class TestParsing:
    @pytest.mark.parametrize(
        "description",
        [
            "m.14A>G",
            "c.4A>G",
            "NC_000999.1:m.14A>G",
            "NC_000999.1(MT-ND2):m.14A>G",
            "NC_000999.1(MT-ND2):c.4A>G",
            "NC_000999.1(ND2):c.4A>G",
            "  c.4A>G  ",
        ],
    )
    def test_accepted_forms(self, plus, description):
        variant = single(plus, description)
        assert variant.m == "NC_000999.1:m.14A>G"
        assert variant.c == "NC_000999.1(MT-ND2):c.4A>G"

    @pytest.mark.parametrize(
        "description, message",
        [
            ("C.4A>G", "lower case"),
            ("M.14A>G", "lower case"),
            ("c.4a>g", "case-sensitive"),
            ("c.4DEL", "case-sensitive"),
            ("c.4_5INSA", "case-sensitive"),
            ("c.04A>G", "leading zeros"),
            ("c.0A>G", "leading zeros"),
            ("c.4 A>G", "spaces"),
            ("c.4_9insA", "flanking (adjacent) positions"),
            ("c.5_4insA", "flanking (adjacent) positions"),
            ("c.4insA", "two flanking positions"),
            ("c.6_4del", "two different positions"),
            ("c.4_4del", "two different positions"),
            ("c.4_6del3", "first and last position"),
            ("c.4del3", "first and last position"),
            ("c.4=", "unchanged"),
            ("c.(4_6)del", "uncertain"),
            ("c.4_5insN", "ambiguous"),
            ("c.4+1G>A", "intronic"),
            ("c.4A>A", "does not change"),
            ("c.4A>G;c.5A>G", "not a supported"),
            ("g.14A>G", "Unsupported coordinate prefix"),
            ("n.4A>G", "Unsupported coordinate prefix"),
            ("p.4A>G", "Unsupported coordinate prefix"),
            ("o.14A>G", "Unsupported coordinate prefix"),
            ("4A>G", "Unsupported mutation notation"),
            ("del.5", "Unsupported mutation notation"),
            ("NC_000999:m.14A>G", "version number"),
            ("NC_000111.1:m.14A>G", "selected genome is NC_000999.1"),
            ("NC_000999.1(MT-ND5):c.4A>G", "selected gene is MT-ND2"),
            ("NC_000999.1:c.4A>G", "has to name the gene"),
            ("m.[14A>G];[15A>G]", "Only one allele"),
            ("m.14A>G(;)15A>G", "unknown phase"),
            ("c.19A>G", "outside of the MT-ND2 coding sequence"),
            ("m.500A>G", "outside of the 38 nt reference"),
        ],
    )
    def test_invalid_descriptions_are_explained(self, plus, description, message):
        assert message in error_of(plus, description)

    def test_reference_nucleotide_is_checked(self, plus):
        assert "is A, not C" in error_of(plus, "c.4C>G")

    def test_a_listed_deleted_sequence_has_to_match(self, plus):
        assert "is AAA, not AAT" in error_of(plus, "c.4_6delAAT")
        assert single(plus, "c.4_6delAAA").c == "NC_000999.1(MT-ND2):c.4_6del"

    def test_a_listed_duplicated_sequence_has_to_match(self, plus):
        assert "not CAA" in error_of(plus, "c.7_9dupCAA")
        assert single(plus, "c.7_9dupCAG").c == "NC_000999.1(MT-ND2):c.7_9dup"

    def test_an_allele_of_variants_in_cis_is_split(self, plus):
        result = describe(plus, "c.[4A>G;10T>C]")
        assert [variant.c for variant in result.variants] == [
            "NC_000999.1(MT-ND2):c.4A>G",
            "NC_000999.1(MT-ND2):c.10T>C",
        ]
        assert result.allele_c == "NC_000999.1(MT-ND2):c.[4A>G;10T>C]"
        assert result.allele_m == "NC_000999.1:m.[14A>G;20T>C]"

    def test_c_numbering_needs_a_complete_cds(self):
        location = FeatureLocation(BeforePosition(10), 28, strand=1)
        context = make_context(UPSTREAM + CDS + DOWNSTREAM, location)

        assert not context.c_available
        assert "incomplete" in error_of(context, "c.4A>G")
        variant = single(context, "m.14A>G")
        assert variant.c is None
        assert describe(context, "m.14A>G").allele_c is None

    def test_valid_and_invalid_inputs_are_reported_separately(self, plus):
        result = describe_variants(plus, ["c.4A>G", "c.4a>g", "c.10T>C"])
        assert [variant.inputs for variant in result.variants] == [(0,), (2,)]
        assert [error.inputs for error in result.errors] == [(1,)]


class TestNormalisation:
    def test_deletion_in_a_run_is_shifted_3_prime(self, plus):
        variant = single(plus, "c.4del")
        assert variant.c == "NC_000999.1(MT-ND2):c.6del"
        assert variant.m == "NC_000999.1:m.16del"

    def test_insertion_copying_the_preceding_base_is_a_duplication(self, plus):
        variant = single(plus, "c.4_5insA")
        assert variant.kind == "duplication"
        assert variant.c == "NC_000999.1(MT-ND2):c.6dup"

    def test_dup_is_shifted_3_prime(self, plus):
        assert single(plus, "c.4dup").c == "NC_000999.1(MT-ND2):c.6dup"

    def test_two_extra_tandem_copies_are_a_repeat(self, plus):
        variant = single(plus, "c.9_10insCAGCAG")
        assert variant.kind == "repeat"
        assert variant.c == "NC_000999.1(MT-ND2):c.7_9CAG[3]"
        assert variant.m == "NC_000999.1:m.17_19CAG[3]"

    def test_repeats_of_units_other_than_3n_are_insertions_in_c(self):
        # HGVS: NM_024312.4:c.1741_1742insTATATATA, not c.1738TA[6].
        cds = "ATGGCTAGCGGCTAA"  # M A S G * with a single TA at c.6_7
        context = make_context(UPSTREAM + cds + DOWNSTREAM, FeatureLocation(10, 25, 1))
        variant = single(context, "c.7_8insTATA")

        assert variant.m == "NC_000999.1:m.16_17TA[3]"
        assert variant.c == "NC_000999.1(MT-ND2):c.7_8insTATA"
        assert "units of 3" in error_of(context, "c.6_7TA[3]")

    def test_one_extra_copy_of_a_whole_repeat_is_a_duplication(self):
        # HGVS: NM_024312.4:c.2692_2693dup, not c.2686A[10].
        cds = "ATGGCTATAGGCTAA"  # TA TA at c.6_9
        context = make_context(UPSTREAM + cds + DOWNSTREAM, FeatureLocation(10, 25, 1))
        variant = single(context, "c.9_10insTATA")

        assert variant.c == "NC_000999.1(MT-ND2):c.6_9dup"
        assert variant.m == "NC_000999.1:m.16_19dup"

    def test_repeat_input_is_understood(self, plus):
        assert single(plus, "c.7_9CAG[3]").c == "NC_000999.1(MT-ND2):c.7_9CAG[3]"
        assert single(plus, "m.17_19CAG[2]").c == "NC_000999.1(MT-ND2):c.7_9dup"
        assert "reference number" in error_of(plus, "c.7_9CAG[1]")

    def test_a_replaced_reverse_complement_is_an_inversion(self, plus):
        variant = single(plus, "c.10_12delinsGAA")
        assert variant.kind == "inversion"
        assert variant.c == "NC_000999.1(MT-ND2):c.10_12inv"
        assert single(plus, "c.10_12inv").m == "NC_000999.1:m.20_22inv"

    def test_one_replaced_nucleotide_is_a_substitution(self, plus):
        assert single(plus, "c.4delinsG").c == "NC_000999.1(MT-ND2):c.4A>G"

    def test_a_delins_is_trimmed_to_the_changed_bases(self, plus):
        assert single(plus, "c.4_6delinsAGA").c == "NC_000999.1(MT-ND2):c.5A>G"

    @pytest.mark.parametrize("description", ["c.4_6delinsAAA", "c.4delinsA"])
    def test_a_description_without_a_change_is_rejected(self, plus, description):
        assert "does not change" in error_of(plus, description)

    def test_a_shift_crosses_the_origin_of_a_circular_genome(self):
        # The stop codon TAA ends at the last genome position and the genome starts with AA.
        cds = "ATG" + "GCT" * 5 + "TAA"
        genome = "AAC" + "G" * 16 + cds
        context = make_context(genome, FeatureLocation(19, 40, 1))

        variant = single(context, "c.21del")
        assert variant.m == "NC_000999.1:m.2del"
        assert variant.c == "NC_000999.1(MT-ND2):c.*2del"
        # The stop codon is still TAA once the poly(A) tail completes it.
        assert variant.p == "YP_000999.1:p.(Ter7=)"
        assert single(context, "m.40_1insT").m == "NC_000999.1:m.40_1insT"

    def test_a_description_shifted_out_of_the_gene_is_still_accepted(self, plus):
        # c.18 is the last A of TAA; the m. form a 3' shifted description could take is
        # still inside the CDS here, so check a variant described from the m. side.
        variant = single(plus, "m.28del")
        assert variant.c == "NC_000999.1(MT-ND2):c.18del"

    def test_insertion_after_the_last_base_is_part_of_the_gene(self, plus):
        assert single(plus, "c.18_*1insT").c == "NC_000999.1(MT-ND2):c.18_*1insT"
        # c.*1 and c.*2 are C, so an inserted C is a duplication shifted out of the CDS.
        assert single(plus, "c.18_*1insC").c == "NC_000999.1(MT-ND2):c.*2dup"

    def test_insertion_in_front_of_c1_is_outside_of_the_cds(self, plus):
        assert "outside of the MT-ND2 coding sequence" in error_of(plus, "c.-1_1insT")

    def test_insertion_in_front_of_c1_equivalent_to_a_dup_of_c1_is_inside(self, plus):
        # c.1 is A, so an A inserted in front of it duplicates c.1 (3' rule).
        variant = single(plus, "c.-1_1insA")
        assert variant.c == "NC_000999.1(MT-ND2):c.1dup"
        assert variant.p == "YP_000999.1:p.(Met1?)"


class TestMerging:
    def test_adjacent_substitutions_are_one_delins(self, plus):
        result = describe(plus, "c.4A>G", "c.5A>T")
        assert [variant.c for variant in result.variants] == ["NC_000999.1(MT-ND2):c.4_5delinsGT"]
        assert result.variants[0].inputs == (0, 1)

    def test_a_substitution_next_to_an_insertion_is_one_delins(self, plus):
        result = describe(plus, "c.10T>G", "c.10_11insCC")
        assert [variant.c for variant in result.variants] == ["NC_000999.1(MT-ND2):c.10delinsGCC"]

    def test_two_changes_one_base_apart_in_one_codon_are_one_delins(self, plus):
        # HGVS: c.235_237delinsTAT (p.Lys79Tyr), not c.[235A>T;237G>T].
        result = describe(plus, "c.4A>T", "c.6A>T")
        assert [variant.c for variant in result.variants] == ["NC_000999.1(MT-ND2):c.4_6delinsTAT"]
        assert result.variants[0].p == "YP_000999.1:p.(Lys2Tyr)"

    def test_two_changes_one_base_apart_in_two_codons_stay_apart(self, plus):
        result = describe(plus, "c.6A>T", "c.8A>T")
        assert [variant.c for variant in result.variants] == [
            "NC_000999.1(MT-ND2):c.6A>T",
            "NC_000999.1(MT-ND2):c.8A>T",
        ]

    def test_deletions_in_one_run_are_one_deletion(self, plus):
        result = describe(plus, "c.4del", "c.6del")
        assert [variant.c for variant in result.variants] == ["NC_000999.1(MT-ND2):c.5_6del"]

    def test_overlapping_variants_are_rejected(self, plus):
        result = describe_variants(plus, ["c.4A>G", "c.4_5del"])
        assert not result.variants
        assert "overlap" in result.errors[0].message
        assert result.errors[0].inputs == (0, 1)

    def test_the_same_variant_twice_is_rejected(self, plus):
        result = describe_variants(plus, ["c.4A>G", "m.14A>G"])
        assert "more than once" in result.errors[0].message

    def test_variants_cancelling_each_other_are_rejected(self, plus):
        result = describe_variants(plus, ["c.4del", "c.6_7insA"])
        assert "cancel each other out" in result.errors[0].message

    def test_analysis_placement_stays_inside_the_cds(self, plus):
        variant = single(plus, "c.4del")
        assert (variant.cds_start, variant.cds_end, variant.alt) == (3, 4, "")
        insertion = single(plus, "c.4_5insA")
        assert (insertion.cds_start, insertion.cds_end, insertion.alt) == (4, 4, "A")


def applied(context, result) -> str:
    """The CDS with every described row applied."""
    from analyser_cli.hgvs_protein import CdsEdit, apply_edits

    return apply_edits(
        context.cds_sequence,
        [CdsEdit(variant.cds_start, variant.cds_end, variant.alt) for variant in result.variants],
    )


class TestAlleles:
    """All inputs are one allele: every change is described in the context of the others."""

    def test_changes_that_only_could_touch_stay_apart(self):
        # HGVS keeps g.[3C>G;7dup] for AGCGTTTAGC > AGGGTTTTAGC (SVD-WG010 was rejected).
        cds = "ATGAGCGTTTAGCTAA"
        context = make_context(UPSTREAM + cds + DOWNSTREAM, FeatureLocation(10, 26, 1))
        result = describe(context, "c.6C>G", "c.7_8insT")

        assert result.allele_c == "NC_000999.1(MT-ND2):c.[6C>G;10dup]"

    def test_a_shift_stops_at_the_next_change(self, plus):
        # c.4del could shift to c.6, but c.6 is changed too, so the two become one delins.
        result = describe(plus, "c.4del", "c.6A>G")

        assert [variant.c for variant in result.variants] == ["NC_000999.1(MT-ND2):c.5_6delinsG"]
        assert applied(plus, result) == "ATGAGCAGTTCGGATAA"

    def test_changes_that_touch_in_some_placement_are_one_change(self, plus):
        # Reference A A A (c.4_6) becomes G A: consecutive nucleotides changed, one
        # delins however the input placed the deletion (as Mutalyzer normalises it).
        for inputs in (["c.4A>G", "c.5del"], ["c.4A>G", "c.6del"], ["c.4_5delinsG"]):
            result = describe(plus, *inputs)

            assert applied(plus, result) == "ATGGACAGTTCGGATAA"
            assert result.allele_c == "NC_000999.1(MT-ND2):c.4_5delinsG"

    def test_consecutive_deletions_are_one_deletion(self, plus, minus):
        for context in (plus, minus):
            result = describe(context, "c.9del", "c.10del")
            assert result.allele_c == "NC_000999.1(MT-ND2):c.9_10del"
            assert len(result.variants) == 1

    def test_merging_does_not_depend_on_the_strand(self, plus, minus):
        assert describe(plus, "c.3G>C", "c.6del").allele_c == "NC_000999.1(MT-ND2):c.3_4delinsC"
        assert describe(minus, "c.3G>C", "c.6del").allele_c == "NC_000999.1(MT-ND2):c.3_4delinsC"

    def test_insertions_meeting_at_one_position_keep_their_order(self, plus):
        # c.9_10insG shifts nowhere; c.8_9insCA shifts onto it (CAG + C/A run).
        result = describe(plus, "c.8_9insT", "c.9_10insC")

        assert applied(plus, result) == "ATGAAACATGCTTCGGATAA"

    def test_a_delins_with_unchanged_bases_is_split(self, plus):
        # c.4_9 AAACAG > GAACAT changes two codons apart: two substitutions.
        result = describe(plus, "c.4_9delinsGAACAT")

        assert [variant.c for variant in result.variants] == [
            "NC_000999.1(MT-ND2):c.4A>G",
            "NC_000999.1(MT-ND2):c.9G>T",
        ]

    def test_a_delins_one_unchanged_base_apart_in_one_codon_stays_one(self, plus):
        assert single(plus, "c.4_6delinsTAT").c == "NC_000999.1(MT-ND2):c.4_6delinsTAT"

    def test_complement_strand_changes_merged_in_c_are_one_variant_in_m(self, minus):
        # Coding A A A (c.4_6) is TTT at m.23_25. In c. the deletion shifted 3' touches
        # c.6A>G (c.5_6delinsG), so m. describes the same single change: TT > C.
        result = describe(minus, "c.4del", "c.6A>G")

        assert result.allele_c == "NC_000999.1(MT-ND2):c.5_6delinsG"
        assert result.allele_m == "NC_000999.1:m.23_24delinsC"
        assert len(result.variants) == 1
        assert applied(minus, result) == "ATGAGCAGTTCGGATAA"

    def test_the_described_allele_reproduces_the_input(self, plus):
        inputs = ["c.4A>G", "c.6del", "c.8_9insTT", "c.13G>T"]
        result = describe(plus, *inputs)
        again = describe(plus, result.allele_c)

        assert applied(plus, again) == applied(plus, result) == "ATGGACATTGTTCTGATAA"
        assert again.allele_m == result.allele_m


class TestProteinAllele:
    def test_consecutive_amino_acid_changes_are_one_delins(self, plus):
        # HGVS: p.Arg76_Cys77delinsSerTrp, not p.[Arg76Ser;Cys77Trp].
        result = describe(plus, "c.4A>G", "c.8A>T")
        assert result.allele_p == "YP_000999.1:p.(Lys2_Gln3delinsGluLeu)"

    def test_separate_changes_are_listed(self, plus):
        result = describe(plus, "c.4A>G", "c.11T>G")
        assert result.allele_p == "YP_000999.1:p.[(Lys2Glu;Phe4Cys)]"

    def test_nothing_after_a_premature_stop_codon_is_described(self, plus):
        # c.5A>T makes AAA > ATA; c.4A>T then TAA ends the protein at codon 2.
        result = describe(plus, "c.4A>T", "c.13G>T")
        assert result.allele_p == "YP_000999.1:p.(Lys2Ter)"

    def test_an_inserted_stop_codon_ends_the_allele(self, plus):
        # An inserted TAG right before Gln3 is an immediate stop: nonsense, p.(Gln3Ter).
        result = describe(plus, "c.6_7insTAG", "c.13G>T")
        assert result.allele_p == "YP_000999.1:p.(Gln3Ter)"

    def test_an_inserted_sequence_ending_with_a_stop_codon_ends_the_allele(self, plus):
        result = describe(plus, "c.6_7insGGGTAG", "c.13G>T")
        assert result.allele_p == "YP_000999.1:p.(Lys2_Gln3insGlyTer)"


class TestSizeLimits:
    @pytest.mark.parametrize(
        "description, message",
        [
            ("m.20T[2000000000]", "longer than the whole reference"),
            ("c.4_5ins" + "A" * 100, "longer than the whole reference"),
            ("c.1" + "0" * 5000 + "A>G", "outside of the reference"),
        ],
    )
    def test_oversized_descriptions_are_rejected_without_building_them(self, plus, description, message):
        assert message in error_of(plus, description)


class TestAgainstBiocommonsParser:
    """Every single-variant description must be valid for an independent HGVS parser."""

    def test_generated_descriptions_parse(self, plus, minus):
        parser = pytest.importorskip("hgvs.parser").Parser()
        descriptions = []
        for context in (plus, minus):
            for text in ["c.4A>G", "c.4del", "c.4_5insT", "c.6dup", "c.4_5delinsGT",
                         "c.10_12inv", "c.18_*1insT", "c.18_*1insC", "c.4_6delinsTAT",
                         "c.13_15del", "c.2T>C", "c.16T>C", "c.13_14del"]:
                for variant in describe(context, text).variants:
                    descriptions.extend([variant.m, variant.c, variant.p])

        for description in descriptions:
            parser.parse_hgvs_variant(description)


def _cached_genomes():
    for folder in GENOME_FOLDERS:
        genomes = sorted(folder.glob("*.gb"))
        if genomes:
            return genomes
    return []


@pytest.mark.parametrize("genome_path", _cached_genomes(), ids=lambda path: path.stem)
def test_every_cached_cds_maps_onto_its_genome(genome_path):
    record = SeqIO.read(genome_path, "gb")
    for feature in record.features:
        if feature.type != "CDS":
            continue
        context = ReferenceContext.from_genbank(record, feature)
        assert context.cds_sequence == str(feature.extract(record.seq)).upper()
        expected_first = feature.location.end if feature.location.strand == -1 else feature.location.start + 1
        assert context.first_genome_position == int(expected_first)


@pytest.fixture
def human():
    path = next((folder / "NC_012920.1.gb" for folder in GENOME_FOLDERS if (folder / "NC_012920.1.gb").exists()), None)
    if path is None:
        pytest.skip("NC_012920.1 is not cached")
    record = SeqIO.read(path, "gb")

    def context(gene):
        feature = next(
            item for item in record.features
            if item.type == "CDS" and item.qualifiers.get("gene") == [gene]
        )
        return ReferenceContext.from_genbank(record, feature)

    return context


@pytest.mark.parametrize(
    "gene, description, m, c, p",
    [
        # LHON, MT-ND6 on the complement strand.
        ("ND6", "m.14484T>C", "NC_012920.1:m.14484T>C", "NC_012920.1(MT-ND6):c.190A>G",
         "YP_003024037.1:p.(Met64Val)"),
        ("ND6", "NC_012920.1(MT-ND6):c.190A>G", "NC_012920.1:m.14484T>C",
         "NC_012920.1(MT-ND6):c.190A>G", "YP_003024037.1:p.(Met64Val)"),
        # NARP, MT-ATP6.
        ("ATP6", "m.8993T>G", "NC_012920.1:m.8993T>G", "NC_012920.1(MT-ATP6):c.467T>G",
         "YP_003024031.1:p.(Leu156Arg)"),
        # LHON, MT-ND4 and MT-ND1.
        ("ND4", "m.11778G>A", "NC_012920.1:m.11778G>A", "NC_012920.1(MT-ND4):c.1019G>A",
         "YP_003024035.1:p.(Arg340His)"),
        ("ND1", "m.3460G>A", "NC_012920.1:m.3460G>A", "NC_012920.1(MT-ND1):c.154G>A",
         "YP_003024026.1:p.(Ala52Thr)"),
        # A T-- stop codon completed by the poly(A) tail is lost.
        ("COX3", "c.784T>C", "NC_012920.1:m.9990T>C", "NC_012920.1(MT-CO3):c.784T>C",
         "YP_003024032.1:p.(Ter262GlnextTer?)"),
        # The last base of MT-ATP8 shifts 3' out of the gene.
        ("ATP8", "c.207del", "NC_012920.1:m.8573del", "NC_012920.1(MT-ATP8):c.*1del",
         "YP_003024030.1:p.(Ter69=)"),
    ],
)
def test_known_human_variants(human, gene, description, m, c, p):
    variant = single(human(gene), description)
    assert (variant.m, variant.c, variant.p) == (m, c, p)


@pytest.mark.parametrize(
    "gene, inputs, allele_m",
    [
        # A shift stopped by the next change: m.3467 must not be read as reference A.
        ("ND1", ["m.3464T>C", "m.3467A>G", "m.3468del"], "NC_012920.1:m.[3464T>C;3467_3468delinsG]"),
        ("ATP8", ["c.40del", "c.42_43insG", "c.44del"], "NC_012920.1:m.[8405del;8408C>G]"),
        # Checked against Mutalyzer: changes that touch in some placement are one change,
        # merged 5' to 3'; changes always one nucleotide apart stay apart (SVD-WG010 rejected).
        ("ATP6", ["m.8541del", "m.8543del"], "NC_012920.1:m.8541_8542del"),
        ("ATP6", ["m.8541G>A", "m.8543del"], "NC_012920.1:m.8541_8542delinsA"),
        ("ATP6", ["m.8577A>C", "m.8579del", "m.8581G>A"], "NC_012920.1:m.[8577del;8581G>A]"),
        ("ATP6", ["m.8540T>A", "m.8541_8542insT"], "NC_012920.1:m.[8540T>A;8543dup]"),
        ("ND6", ["c.22del", "c.24G>A"], "NC_012920.1:m.14650_14651delinsT"),
    ],
)
def test_review_examples_reproduce_the_input(human, gene, inputs, allele_m):
    context = human(gene)
    result = describe(context, *inputs)
    again = describe(context, result.allele_m)

    assert result.allele_m == allele_m
    assert applied(context, again) == applied(context, result)


@pytest.mark.parametrize(
    "gene, inputs, allele_p",
    [
        ("ATP6", ["c.6_7insGGGTAA", "m.8993T>G"], "YP_003024031.1:p.(Asn2_Glu3insGlyTer)"),
        ("ND1", ["c.30del", "c.600_601insA", "c.700A>G"], "YP_003024026.1:p.(Ile10MetfsTer5)"),
    ],
)
def test_review_examples_of_the_protein_allele(human, gene, inputs, allele_p):
    assert describe(human(gene), *inputs).allele_p == allele_p


def test_hgvs_version_is_declared():
    assert HGVS_VERSION == "21.1"
