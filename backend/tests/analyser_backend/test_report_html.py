from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from backend.analyser_backend.routers import report
from backend.analyser_cli.hgvs import (
    DescriptionResult,
    ReferenceContext,
    VariantDescription,
    describe_variants,
    reverse_complement,
)
from backend.analyser_cli.hgvs_protein import FRAMESHIFT, MISSENSE, ProteinConsequence

# Codons M K Q F G * of an ND2 gene on genome positions 3914..3931.
CDS = "ATGAAACAGTTCGGATAA"
GENETIC_CODE = "Vertebrate Mitochondrial (NCBI transl_table=2)"


def make_context(strand: int = 1) -> ReferenceContext:
    gene = CDS if strand == 1 else reverse_complement(CDS)
    record = SeqRecord(
        Seq("GC" * 1956 + "G" + gene + "CCGCCGCCGC"),
        id="NC_002008.4",
        annotations={"topology": "circular"},
    )
    location = FeatureLocation(3913, 3913 + len(CDS), strand=strand)
    record.features = [
        SeqFeature(
            location,
            type="CDS",
            qualifiers={"gene": ["ND2"], "transl_table": ["2"], "protein_id": ["NP_008472.1"]},
        ),
        SeqFeature(
            location,
            type="gene",
            qualifiers={"gene": ["ND2"], "nomenclature": ["Official Symbol: MT-ND2 | Name: x"]},
        ),
    ]
    return ReferenceContext.from_genbank(record, record.features[0])


def variant(protein: ProteinConsequence | None = None) -> VariantDescription:
    return VariantDescription(
        inputs=(0,),
        texts=("m.3917A>G",),
        kind="substitution",
        m_edits=("3917A>G",),
        c_edits=("4A>G",),
        m="NC_002008.4:m.3917A>G",
        c="NC_002008.4(MT-ND2):c.4A>G",
        cds_start=3,
        cds_end=4,
        ref="A",
        alt="G",
        protein=protein,
        p=protein.description("NP_008472.1") if protein else None,
    )


def render(original, mutated, result=None, context=None):
    context = context or make_context()
    result = result or DescriptionResult(variants=(variant(),), errors=())
    return report._build_report_html(
        "dog",
        "ND2",
        "Detailed report",
        context,
        result,
        original,
        mutated_amino_acid_translation=mutated,
        genetic_code=GENETIC_CODE,
        original_pdb_path=None,
        mutated_pdb_path=None,
        mutated_prediction=None,
    )


def test_build_report_html_includes_reference_metadata_and_hgvs_descriptions():
    context = make_context()
    result = describe_variants(context, ["m.3914A>G", "c.4_5insTT", "c.7_8del"])
    html_report = render("M" * 200, "M" * 100 + "K" + "M" * 99, result, context)

    assert "<strong>Genome Range:</strong> 3914..3931" in html_report
    assert "<strong>Reference:</strong> NC_002008.4" in html_report
    assert "<strong>Protein:</strong> NP_008472.1" in html_report
    assert "<strong>Gene:</strong> ND2 (MT-ND2)" in html_report
    assert "<strong>HGVS Nomenclature:</strong> version 21.1" in html_report
    assert f"<strong>Genetic Code:</strong> {GENETIC_CODE}" in html_report
    assert (
        "<th>Genome (m.)</th><th>Coding DNA (c.)</th><th>Protein (p., predicted)</th>"
        in html_report
    )
    assert (
        "<tr><td>NC_002008.4:m.3914A&gt;G</td><td>NC_002008.4(MT-ND2):c.1A&gt;G</td>"
        "<td>NP_008472.1:p.(Met1?)</td><td>substitution</td>"
        "<td>Start codon changed (consequence unknown)</td></tr>"
    ) in html_report
    assert "<td>NC_002008.4:m.3917_3918insTT</td>" in html_report
    assert "<td>NC_002008.4:m.3920_3921del</td>" in html_report
    # All variants together, described as one allele.
    assert "<code>NC_002008.4:m.[3914A&gt;G;3917_3918insTT;3920_3921del]</code>" in html_report
    assert "<code>NC_002008.4(MT-ND2):c.[1A&gt;G;4_5insTT;7_8del]</code>" in html_report
    assert "<code>NP_008472.1:p.(Met1?)</code>" in html_report
    assert "Descriptions follow HGVS Nomenclature version 21.1" in html_report
    # Variant types (DNA) and predicted consequences (protein) are counted separately.
    assert "<li><span>substitution</span><strong>1</strong></li>" in html_report
    assert "<li><span>Frameshift</span><strong>2</strong></li>" in html_report


def test_build_report_html_shows_a_complement_location():
    context = make_context(strand=-1)
    result = describe_variants(context, ["c.4A>G"])
    html_report = render("MKQFG*", "MEQFG*", result, context)

    assert "<strong>Genome Range:</strong> complement(3914..3931)" in html_report
    assert "<td>NC_002008.4:m.3928T&gt;C</td><td>NC_002008.4(MT-ND2):c.4A&gt;G</td>" in html_report


def test_build_report_html_compares_the_two_translations():
    html_report = render("MKWA", "MKWG")

    assert '<span class="aln-label">Original</span>' in html_report
    assert '<span class="aln-label">Mutated</span>' in html_report
    # Only the fourth residue changed.
    assert 'MKW<span class="aln-diff">A</span>' in html_report
    assert 'MKW<span class="aln-diff">G</span>' in html_report
    assert "<li><span>Changed residues</span><strong>1</strong></li>" in html_report
    assert "<li><span>Unchanged residues</span><strong>3</strong></li>" in html_report


def test_build_report_html_marks_a_deleted_residue_as_a_gap():
    html_report = render("MKWA", "MKA")

    assert 'MK<span class="aln-gap">-</span>A' in html_report
    assert "<li><span>Original length</span><strong>4 aa</strong></li>" in html_report
    assert "<li><span>Mutated length</span><strong>3 aa</strong></li>" in html_report


def test_build_report_html_groups_residues_by_five_and_wraps_lines():
    html_report = render("M" * 130, "K" + "M" * 129)

    # 130 residues wrap into lines of 60, so three blocks numbered 1, 61 and 121.
    assert html_report.count('<div class="aln-block">') == 3
    assert '<span class="aln-pos">61</span>' in html_report
    assert '<span class="aln-pos">121</span>' in html_report
    assert '<span class="aln-pos aln-pos-end">120</span>' in html_report
    # A space is inserted every five residues.
    assert 'aln-diff">K</span>MMMM MMMMM' in html_report


def test_build_report_html_does_not_align_a_frameshift_tail():
    frameshift = ProteinConsequence("Cys5GlyfsTer6", FRAMESHIFT, 5)
    result = DescriptionResult(
        variants=(variant(ProteinConsequence("Lys2Gln", MISSENSE, 2)), variant(frameshift)),
        errors=(),
    )
    html_report = render("MKWACDEF*", "MQWAGHIKL*", result)

    assert "The frameshift changes the protein from residue 5 on" in html_report
    assert '<div class="aln-line aln-divergence">' in html_report
    # The substitution happens before the frameshift, so it stays a plain substitution.
    assert 'M<span class="aln-diff">K</span>WA' in html_report
    assert 'M<span class="aln-diff">Q</span>WA' in html_report
    # The tail is read in another frame, so its residues are not reported as conserved.
    assert "<li><span>Changed residues</span><strong>7</strong></li>" in html_report
    assert "<li><span>Unchanged residues</span><strong>3</strong></li>" in html_report


def test_build_report_html_shows_a_premature_stop_codon_as_a_truncation():
    html_report = render("MKWACDE*", "MKWA*")

    assert (
        "A premature stop codon at codon 5 ends the mutated peptide after 4 amino acids, "
        "so the remaining 3 residues of the original are lost" in html_report
    )
    # Lengths count amino acids up to the first stop codon (MKWACDE / MKWA).
    assert "<li><span>Original length</span><strong>7 aa</strong></li>" in html_report
    assert "<li><span>Mutated length</span><strong>4 aa</strong></li>" in html_report
    # The lost tail is padded at the end instead of being spread through the alignment.
    assert (
        'MKWA<span class="aln-diff">*</span> <span class="aln-gap">-</span>'
        '<span class="aln-gap">-</span><span class="aln-gap">-</span>'
    ) in html_report


def test_build_report_html_keeps_an_in_frame_deletion_aligned():
    html_report = render("MKWACDE*", "MKWCDE*")

    assert "The frameshift changes" not in html_report
    assert '<div class="aln-line aln-divergence">' not in html_report
    assert "<li><span>Changed residues</span><strong>1</strong></li>" in html_report


def test_build_report_html_explains_missing_c_numbering():
    context = make_context()
    object.__setattr__(context, "cds_complete", False)
    result = describe_variants(context, ["m.3917A>G"])
    html_report = render("MKQFG*", "MEQFG*", result, context)

    assert "<td>NC_002008.4:m.3917A&gt;G</td><td>n/a</td>" in html_report
    assert "not available: the ND2 CDS is annotated as incomplete" in html_report


def test_parse_report_html_reads_mutations_from_extended_table():
    html_text = """
    <html>
      <body>
        <section class="card">
          <h2>Metadata</h2>
          <div><strong>Organism:</strong> dog</div>
          <div><strong>Gene:</strong> ND2</div>
          <div><strong>Genome Range:</strong> 3914..4955</div>
          <div><strong>Generated:</strong> 2026-08-10T12:00:00+00:00</div>
        </section>
        <table>
          <tbody>
            <tr><td>m.3914A>G</td><td>1</td><td>3914</td><td>SUB</td><td>Substitution</td></tr>
          </tbody>
        </table>
      </body>
    </html>
    """

    job = report._parse_report_html("job-1", html_text, "http://127.0.0.1:8000")

    assert job is not None
    assert job["mutations"] == ["m.3914A>G"]


def test_parse_report_html_reads_the_hgvs_descriptions_back():
    context = make_context()
    result = describe_variants(context, ["m.3914A>G", "c.7_8del"])
    html_report = render("MKQFG*", "MEQFG*", result, context)

    job = report._parse_report_html("job-2", html_report, "http://127.0.0.1:8000")

    assert job["mutations"] == ["NC_002008.4:m.3914A>G", "NC_002008.4:m.3920_3921del"]
