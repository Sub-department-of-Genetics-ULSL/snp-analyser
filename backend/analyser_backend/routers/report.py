import html
import json
import re
import sqlite3
import threading
import uuid
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

from fastapi import APIRouter, BackgroundTasks, HTTPException, Request

from backend.analyser_backend.models.report_data import DescribeRequest, ReportData
from backend.analyser_cli import Reporter, Analyser, DataManager
from backend.analyser_cli.hgvs import (
    HGVS_VERSION,
    DescriptionResult,
    ReferenceContext,
    VariantDescription,
    describe_variants,
)
from backend.analyser_cli.hgvs_protein import FRAMESHIFT
from backend.analyser_cli.protein_alignment import (
    GAP as ALIGNMENT_GAP,
    align_translations,
    group_residues,
)
from backend.analyser_backend.services.helixfold import (
    HelixFoldPrediction,
    HelixFoldService,
    protein_sequence_for_prediction,
)

router = APIRouter()
manager = DataManager()
reporter = Reporter()
API_PREFIX = "/api"
REPORTS_DIR = Path(__file__).resolve().parent.parent / "generated_reports"
PDB_FILES_DIR = Path(__file__).resolve().parent.parent / "pdb_files"
TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "templates"
REPORT_TEMPLATE_PATH = TEMPLATES_DIR / "report.html"
REPORTS_DIR.mkdir(parents=True, exist_ok=True)
PDB_FILES_DIR.mkdir(parents=True, exist_ok=True)
REPORT_JOBS = {}
REPORT_JOBS_LOCK = threading.Lock()
HELIXFOLD_SERVICE: HelixFoldService | None = None
HELIXFOLD_LOAD_ERROR: str | None = None
try:
    HELIXFOLD_SERVICE = HelixFoldService()
except Exception as exc:
    HELIXFOLD_LOAD_ERROR = f"{type(exc).__name__}: {exc}"

# The protein alignment is wrapped into lines of this many residues, grouped by five,
# so the report never needs a horizontal scrollbar.
ALIGNMENT_LINE_LENGTH = 60
ALIGNMENT_GROUP_SIZE = 5


def get_helixfold_startup_status() -> tuple[bool, str]:
    if HELIXFOLD_LOAD_ERROR is not None:
        return False, f"Failed to initialize HelixFold service: {HELIXFOLD_LOAD_ERROR}"

    if HELIXFOLD_SERVICE is None:
        return False, "HelixFold service is unavailable: unknown initialization error."

    if HELIXFOLD_SERVICE.available:
        return True, "HelixFold service is available."

    reason = HELIXFOLD_SERVICE.unavailable_reason or "Unknown HelixFold-single configuration error."
    return False, reason


def get_reference_context(organism: str, gene: str) -> ReferenceContext:
    """HGVS coordinates of a gene; ``LookupError`` when the organism or gene is unknown."""
    try:
        return manager.get_reference_context(organism, gene)
    except KeyError as exc:
        raise LookupError(f"Unknown organism '{organism}'") from exc
    except ValueError as exc:
        if "not found" in str(exc):
            raise LookupError(str(exc)) from exc
        raise


def describe_mutations(
    organism: str, gene: str, mutations: list[str]
) -> tuple[ReferenceContext, DescriptionResult]:
    """Validate and normalise the HGVS descriptions submitted for one gene."""
    context = get_reference_context(organism, gene)
    return context, describe_variants(context, mutations)


def _analyser_mutations(variants: tuple[VariantDescription, ...]) -> list[dict]:
    """Turn normalised variants into the single-base edits the Analyser applies."""
    mutations = []
    for variant in variants:
        start, end = variant.cds_start, variant.cds_end
        if end - start == 1 and len(variant.alt) == 1:
            mutations.append(
                {"type": "sub", "position": start + 1, "ref": variant.ref, "alt": variant.alt}
            )
            continue

        if variant.alt:
            # The Analyser applies edits right to left and keeps the given order for one
            # position, so the new bases go in after the last replaced base before that
            # base is deleted - exactly where the replaced stretch was.
            mutations.append({"type": "ins", "position": max(end, start), "alt": variant.alt})
        mutations.extend(
            {"type": "del", "position": position} for position in range(end, start, -1)
        )
    return mutations


def _report_mutation(variant: VariantDescription) -> dict:
    """The variant as the report, the text sub-reports and the HelixFold cache see it."""
    insertion = variant.cds_start == variant.cds_end
    return {
        "type": variant.kind,
        "position": variant.cds_start if insertion else variant.cds_start + 1,
        "end_position": variant.cds_end,
        "ref": variant.ref,
        "alt": variant.alt,
        "hgvs_m": variant.m,
        "hgvs_c": variant.c,
        "hgvs_p": variant.p,
        "consequence": variant.consequence,
        "protein_kind": variant.protein.kind if variant.protein else None,
        "first_residue": variant.protein.first_residue if variant.protein else None,
        "length_change": variant.length_change,
        "frameshift": variant.frameshift,
    }


def _frameshift_residue(variants: tuple[VariantDescription, ...]) -> int | None:
    """The first amino acid a frameshift changes, which is where HGVS starts its ``fs``."""
    for variant in variants:
        if variant.protein is not None and variant.protein.kind == FRAMESHIFT:
            return variant.protein.first_residue
    return None


def _genome_range(context: ReferenceContext) -> str:
    first, last = context.first_genome_position, context.last_genome_position
    if context.strand == -1:
        return f"complement({last}..{first})"
    return f"{first}..{last}"


def _build_static_url(path: Path) -> str:
    relative_path = path.relative_to(PDB_FILES_DIR).as_posix()
    return f"{API_PREFIX}/pdb-files/{relative_path}"


def _build_original_pdb_path(organism: str, gene: str) -> Path:
    return PDB_FILES_DIR / organism / f"{gene}.pdb"


@lru_cache(maxsize=1)
def _load_report_template() -> str:
    return REPORT_TEMPLATE_PATH.read_text(encoding="utf-8")


def _render_report_template(context: dict[str, str]) -> str:
    template = _load_report_template()
    for key, value in context.items():
        template = template.replace(f"__{key}__", value)
    return template


def _build_alignment_residues(residues: str, differences: tuple[bool, ...]) -> str:
    """Render one alignment line, marking the residues that changed."""
    rendered = []

    for index, residue in enumerate(residues):
        if index and index % ALIGNMENT_GROUP_SIZE == 0:
            rendered.append(" ")

        escaped = html.escape(residue)
        if index < len(differences) and differences[index]:
            css_class = "aln-gap" if residue == ALIGNMENT_GAP else "aln-diff"
            rendered.append(f'<span class="{css_class}">{escaped}</span>')
        else:
            rendered.append(escaped)

    return "".join(rendered)


def _build_alignment_line(label: str, residues: str, differences, start, end) -> str:
    position_start = str(start) if start is not None else ""
    position_end = str(end) if end is not None else ""

    return (
        '<div class="aln-line">'
        f'<span class="aln-label">{html.escape(label)}</span>'
        f'<span class="aln-pos">{position_start}</span>'
        f'<span class="aln-seq">{_build_alignment_residues(residues, differences)}</span>'
        f'<span class="aln-pos aln-pos-end">{position_end}</span>'
        "</div>"
    )


def _build_translation_alignment_html(alignment) -> str:
    """Render the original and mutated translation as wrapped, aligned blocks."""
    if not alignment.original:
        return '<p class="muted">The gene has no protein translation.</p>'

    if not alignment.has_mutated:
        return '<p class="muted">No mutations were applied, so there is nothing to compare.</p>'

    blocks = [_build_divergence_notice(alignment)]
    for line in alignment.lines(ALIGNMENT_LINE_LENGTH, ALIGNMENT_GROUP_SIZE):
        marker = group_residues(
            "".join("|" if different else " " for different in line.differences),
            ALIGNMENT_GROUP_SIZE,
        )
        blocks.append(
            '<div class="aln-block">'
            + _build_alignment_line(
                "Original", line.original, line.differences, line.original_start, line.original_end
            )
            + _build_alignment_line(
                "Mutated", line.mutated, line.differences, line.mutated_start, line.mutated_end
            )
            + '<div class="aln-line aln-marker">'
            '<span class="aln-label"></span><span class="aln-pos"></span>'
            f'<span class="aln-seq">{marker}</span>'
            '<span class="aln-pos aln-pos-end"></span>'
            "</div>"
            + _build_divergence_marker(line)
            + "</div>"
        )

    return f'<div class="alignment">{"".join(blocks)}</div>'


def _build_divergence_marker(line) -> str:
    """Point at the residue where the mutated translation stops being homologous."""
    if line.divergence is None:
        return ""

    arrow = "".join("^" if index == line.divergence else " " for index in range(len(line.original)))

    return (
        '<div class="aln-line aln-divergence">'
        '<span class="aln-label"></span><span class="aln-pos"></span>'
        f'<span class="aln-seq">{group_residues(arrow, ALIGNMENT_GROUP_SIZE)}</span>'
        '<span class="aln-pos aln-pos-end"></span>'
        "</div>"
    )


def _build_divergence_notice(alignment) -> str:
    """Warn that the residues past a frameshift or a premature stop cannot be compared."""
    if alignment.is_homologous:
        return ""

    if alignment.divergence_reason == "truncation":
        kept = alignment.mutated_length
        lost = alignment.original_length - kept
        explanation = (
            f"A premature stop codon at codon {kept + 1} ends the mutated peptide after "
            f"{kept} amino acid{'' if kept == 1 else 's'}, so the remaining {lost} "
            f"residue{' of the original is' if lost == 1 else 's of the original are'} "
            f"lost rather than changed."
        )
    else:
        explanation = (
            f"The frameshift changes the protein from residue {alignment.divergence_start} on "
            f"(the first amino acid changed, where HGVS starts its 'fs' description). "
            f"Everything after it is read in another frame, so the two tails are different "
            f"peptides and are shown side by side instead of aligned residue by residue."
        )

    return f'<p class="aln-notice">{html.escape(explanation)}</p>'


def _count_items(values: list[str]) -> str:
    counts: dict[str, int] = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return "".join(
        f"<li><span>{html.escape(key)}</span><strong>{count}</strong></li>"
        for key, count in sorted(counts.items())
    )


def _build_hgvs_summary(context: ReferenceContext, result: DescriptionResult) -> str:
    """The variants of the report together as one allele, with the references used."""
    coding = result.allele_c
    if not context.c_available:
        coding = f"not available: {context.c_unavailable_reason}"
    rows = [
        ("Genome (m.)", result.allele_m),
        ("Coding DNA (c.)", coding),
        ("Protein (p., predicted)", result.allele_p),
    ]
    items = "".join(
        f"<li><span>{html.escape(label)}</span><code>{html.escape(value or 'n/a')}</code></li>"
        for label, value in rows
    )
    references = [f"{context.accession} (genome, m.)"]
    if context.c_available:
        references.append(f"{context.symbol} CDS annotated in {context.accession} (c.)")
    if context.protein_id:
        references.append(f"{context.protein_id} (protein, p.)")
    note = (
        f"All variants are in cis (one mitochondrial DNA molecule). Reference sequences: "
        f"{', '.join(references)}. Descriptions follow HGVS Nomenclature version {HGVS_VERSION}; "
        f"protein changes are predicted, hence the parentheses."
    )
    return f'<ul class="hgvs-list">{items}</ul><p class="muted">{html.escape(note)}</p>'


def _build_report_html(
    organism: str,
    gene: str,
    report_text: str,
    context: ReferenceContext,
    result: DescriptionResult,
    amino_acid_translation: str,
    *,
    mutated_amino_acid_translation: str,
    genetic_code: str,
    original_pdb_path: Path | None,
    mutated_pdb_path: Path | None,
    mutated_prediction: HelixFoldPrediction | None,
):
    variants = result.variants
    original_pdb_url = _build_static_url(original_pdb_path) if original_pdb_path else None
    mutated_pdb_url = _build_static_url(mutated_pdb_path) if mutated_pdb_path else None
    prediction_available = mutated_pdb_url is not None
    alignment = align_translations(
        amino_acid_translation,
        mutated_amino_acid_translation,
        frameshift_residue=_frameshift_residue(variants),
    )

    variant_rows = [
        f"<tr><td>{html.escape(variant.m)}</td>"
        f"<td>{html.escape(variant.c or 'n/a')}</td>"
        f"<td>{html.escape(variant.p or 'n/a')}</td>"
        f"<td>{html.escape(variant.kind)}</td>"
        f"<td>{html.escape(variant.consequence)}</td></tr>"
        for variant in variants
    ]
    variant_table_rows = "\n".join(variant_rows) or "<tr><td colspan='5'>No variants</td></tr>"

    prediction_summary = ""
    if mutated_prediction is not None:
        prediction_summary = f"""
        <section class=\"card\">
          <h2>HelixFold cache</h2>
          <ul class=\"stat-list\">
            <li><span>Available</span><strong>yes</strong></li>
            <li><span>Cache key</span><strong>{html.escape(mutated_prediction.cache_key[:16])}…</strong></li>
            <li><span>Protein length</span><strong>{len(mutated_prediction.protein_sequence)}</strong></li>
          </ul>
        </section>
        """
    else:
        prediction_summary = """
        <section class=\"card\">
          <h2>HelixFold cache</h2>
          <p class=\"muted\">Mutated structure prediction was not requested for this report.</p>
        </section>
        """

    viewer_payload = {
        "original": {
            "url": original_pdb_url,
            "available": original_pdb_url is not None,
            "label": "Original protein",
        },
        "mutated": {
            "url": mutated_pdb_url,
            "available": prediction_available,
            "label": "Mutated protein",
        },
        "metadata": mutated_prediction.metadata if mutated_prediction is not None else None,
    }

    viewer_initial_note = (
        "Original and mutated structures are available."
        if prediction_available and original_pdb_url
        else "Only the original structure is available for this report."
    )

    return _render_report_template(
        {
            "TITLE_ORGANISM": html.escape(organism),
            "TITLE_GENE": html.escape(gene),
            "ORGANISM": html.escape(organism),
            "GENE": html.escape(gene if gene == context.symbol else f"{gene} ({context.symbol})"),
            "REFERENCE_SEQUENCE": html.escape(context.accession),
            "PROTEIN_REFERENCE": html.escape(context.protein_id or "n/a"),
            "GENOME_RANGE": html.escape(_genome_range(context)),
            "GENETIC_CODE": html.escape(genetic_code),
            "HGVS_VERSION": html.escape(HGVS_VERSION),
            "GENERATED_AT": html.escape(datetime.now(timezone.utc).isoformat()),
            "TOTAL_VARIANTS": str(len(variants)),
            "ORIGINAL_PDB_STATUS": "available" if original_pdb_url else "missing",
            "MUTATED_PDB_STATUS": "available" if prediction_available else "missing",
            "VARIANT_TYPE_ITEMS": _count_items([variant.kind for variant in variants]),
            "CONSEQUENCE_ITEMS": _count_items([variant.consequence for variant in variants]),
            "HGVS_SUMMARY": _build_hgvs_summary(context, result),
            "PREDICTION_SUMMARY": prediction_summary,
            "AMINO_ACID_TRANSLATION": html.escape(amino_acid_translation),
            "MUTATED_AMINO_ACID_TRANSLATION": html.escape(
                mutated_amino_acid_translation or "not available"
            ),
            "TRANSLATION_ALIGNMENT": _build_translation_alignment_html(alignment),
            "ORIGINAL_PROTEIN_LENGTH": str(alignment.original_length),
            "MUTATED_PROTEIN_LENGTH": str(alignment.mutated_length),
            "CHANGED_RESIDUES": str(alignment.changed_residues),
            "IDENTICAL_RESIDUES": str(alignment.identical_residues),
            "VARIANT_TABLE_ROWS": variant_table_rows,
            "ORIGINAL_BUTTON_DISABLED": "disabled" if not original_pdb_url else "",
            "MUTATED_BUTTON_DISABLED": "disabled" if not prediction_available else "",
            "OVERLAY_BUTTON_DISABLED": (
                "disabled" if not (prediction_available and original_pdb_url) else ""
            ),
            "VIEWER_INITIAL_NOTE": viewer_initial_note,
            "REPORT_TEXT": html.escape(report_text),
            "VIEWER_PAYLOAD": json.dumps(viewer_payload),
        }
    )


def _create_job_payload(job_id: str, data: ReportData, report_base_url: str):
    now = datetime.now(timezone.utc).isoformat()
    return {
        "id": job_id,
        "status": "pending",
        "organism": data.organism,
        "gene": data.gene,
        "mutations": data.mutations,
        "predict_mutated_structure": data.predict_mutated_structure,
        "created_at": now,
        "updated_at": now,
        "report_url": None,
        "report_base_url": report_base_url,
        "error": None,
    }


def _report_html_path(job_id: str) -> Path:
    return REPORTS_DIR / f"{job_id}.html"


def _report_metadata_path(job_id: str) -> Path:
    return REPORTS_DIR / f"{job_id}.json"


def _report_jobs_db_path() -> Path:
    return REPORTS_DIR.parent / "report_jobs.sqlite3"


def _open_report_jobs_db() -> sqlite3.Connection:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(_report_jobs_db_path())
    connection.row_factory = sqlite3.Row
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS report_jobs (
            job_id TEXT PRIMARY KEY,
            payload_json TEXT NOT NULL,
            report_file TEXT,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
        """
    )
    return connection


def _serialize_persisted_job(job: dict) -> dict:
    return {
        "id": job["id"],
        "status": job["status"],
        "organism": job["organism"],
        "gene": job["gene"],
        "mutations": job["mutations"],
        "predict_mutated_structure": job.get("predict_mutated_structure", False),
        "created_at": job["created_at"],
        "updated_at": job["updated_at"],
        "report_url": job.get("report_url"),
        "report_base_url": job.get("report_base_url"),
        "error": job.get("error"),
    }


def _serialize_job(job: dict):
    return {
        "id": job["id"],
        "status": job["status"],
        "organism": job["organism"],
        "gene": job["gene"],
        "mutations": job["mutations"],
        "predict_mutated_structure": job.get("predict_mutated_structure", False),
        "created_at": job["created_at"],
        "updated_at": job["updated_at"],
        "report_url": job["report_url"],
        "report_base_url": job.get("report_base_url"),
        "error": job["error"],
    }


def _persist_report_job_to_db(job: dict, report_file_name: str | None = None) -> None:
    persisted_job = _serialize_persisted_job(job) | {
        "status": job["status"],
        "report_file": report_file_name,
    }
    now = datetime.now(timezone.utc).isoformat()
    with _open_report_jobs_db() as connection:
        connection.execute(
            """
            INSERT INTO report_jobs (job_id, payload_json, report_file, created_at, updated_at)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(job_id) DO UPDATE SET
                payload_json = excluded.payload_json,
                report_file = excluded.report_file,
                created_at = excluded.created_at,
                updated_at = excluded.updated_at
            """,
            (
                job["id"],
                json.dumps(persisted_job, sort_keys=True),
                report_file_name,
                job["created_at"],
                now,
            ),
        )


def _load_report_job_from_db(job_id: str, report_base_url: str) -> dict | None:
    if not _report_jobs_db_path().exists():
        return None

    with _open_report_jobs_db() as connection:
        row = connection.execute(
            "SELECT payload_json FROM report_jobs WHERE job_id = ?",
            (job_id,),
        ).fetchone()

    if row is None:
        return None

    data = json.loads(row["payload_json"])
    if data.get("report_base_url") != report_base_url:
        data["report_base_url"] = report_base_url
        if data.get("status") == "completed":
            data["report_url"] = f"{report_base_url}{API_PREFIX}/generated-reports/{job_id}.html"
    return data


def _load_report_jobs_from_db(report_base_url: str) -> dict[str, dict]:
    if not _report_jobs_db_path().exists():
        return {}

    with _open_report_jobs_db() as connection:
        rows = connection.execute(
            "SELECT job_id, payload_json FROM report_jobs ORDER BY updated_at DESC"
        ).fetchall()

    jobs: dict[str, dict] = {}
    for row in rows:
        data = json.loads(row["payload_json"])
        data["report_base_url"] = report_base_url
        if data.get("status") == "completed":
            data["report_url"] = f"{report_base_url}{API_PREFIX}/generated-reports/{row['job_id']}.html"
        jobs[row["job_id"]] = data
    return jobs


def _parse_report_html(job_id: str, html_text: str, report_base_url: str) -> dict | None:
    organism_match = re.search(r"<strong>Organism:</strong>\s*([^<]+)</div>", html_text)
    gene_match = re.search(r"<strong>Gene:</strong>\s*([^<]+)</div>", html_text)
    generated_match = re.search(r"<strong>Generated:</strong>\s*([^<]+)</div>", html_text)
    if organism_match is None and gene_match is None:
        return None

    mutation_rows = re.findall(r"<tr>\s*<td>(.*?)</td>.*?</tr>", html_text, flags=re.DOTALL)
    mutations = [html.unescape(row).strip() for row in mutation_rows if row.strip()]

    predicted = bool(re.search(r"<span>Mutated PDB</span>\s*<strong>available</strong>", html_text))

    timestamp = generated_match.group(1).strip() if generated_match else datetime.now(timezone.utc).isoformat()
    job = {
        "id": job_id,
        "status": "completed",
        "organism": html.unescape(organism_match.group(1).strip()) if organism_match else None,
        "gene": html.unescape(gene_match.group(1).strip()) if gene_match else None,
        "mutations": mutations,
        "predict_mutated_structure": predicted,
        "created_at": timestamp,
        "updated_at": timestamp,
        "report_url": f"{report_base_url}{API_PREFIX}/generated-reports/{job_id}.html",
        "report_base_url": report_base_url,
        "error": None,
    }
    return job


def _seed_report_job_from_file(job_id: str, report_base_url: str) -> dict | None:
    metadata_path = _report_metadata_path(job_id)
    if metadata_path.exists():
        data = json.loads(metadata_path.read_text(encoding="utf-8"))
        data["report_base_url"] = report_base_url
        if _report_html_path(job_id).exists():
            data["report_url"] = f"{report_base_url}{API_PREFIX}/generated-reports/{job_id}.html"
        return data

    html_path = _report_html_path(job_id)
    if not html_path.exists():
        return None

    return _parse_report_html(job_id, html_path.read_text(encoding="utf-8"), report_base_url)
def _persist_completed_job(job_id: str, job: dict, report_file_path: Path) -> None:
    metadata_path = _report_metadata_path(job_id)
    persisted_job = _serialize_persisted_job(job) | {
        "status": "completed",
        "report_file": report_file_path.name,
    }
    metadata_path.write_text(json.dumps(persisted_job, indent=2, sort_keys=True), encoding="utf-8")
    _persist_report_job_to_db(persisted_job, report_file_path.name)


def _job_from_metadata(job_id: str, metadata_path: Path, report_base_url: str) -> dict | None:
    if not metadata_path.exists():
        return None

    data = json.loads(metadata_path.read_text(encoding="utf-8"))
    html_path = _report_html_path(job_id)
    if html_path.exists():
        data["report_url"] = f"{report_base_url}{API_PREFIX}/generated-reports/{job_id}.html"
    _persist_report_job_to_db(data, data.get("report_file"))
    return data


def _job_from_report_file(job_id: str, report_base_url: str) -> dict | None:
    html_path = _report_html_path(job_id)
    if not html_path.exists():
        return None

    parsed_job = _parse_report_html(job_id, html_path.read_text(encoding="utf-8"), report_base_url)
    if parsed_job is not None:
        _persist_report_job_to_db(parsed_job, html_path.name)
    return parsed_job


def _collect_report_jobs(report_base_url: str):
    jobs = _load_report_jobs_from_db(report_base_url)

    for report_file in REPORTS_DIR.glob("*.html"):
        job_id = report_file.stem
        if job_id in jobs:
            continue
        legacy_job = _job_from_report_file(job_id, report_base_url)
        if legacy_job is not None:
            jobs[job_id] = legacy_job

    for metadata_file in REPORTS_DIR.glob("*.json"):
        job_id = metadata_file.stem
        if job_id in jobs:
            continue
        metadata_job = _job_from_metadata(job_id, metadata_file, report_base_url)
        if metadata_job is not None:
            jobs[job_id] = metadata_job

    with REPORT_JOBS_LOCK:
        jobs.update(REPORT_JOBS)
        for job in REPORT_JOBS.values():
            if job.get("status") in {"pending", "running", "completed", "failed"}:
                _persist_report_job_to_db(job, _report_html_path(job["id"]).name if _report_html_path(job["id"]).exists() else None)

    return sorted(jobs.values(), key=lambda item: item["created_at"], reverse=True)


def _predict_mutated_structure(
    *,
    organism: str,
    gene: str,
    gene_sequence: str,
    report_mutations: list[dict],
    analyser: Analyser,
) -> HelixFoldPrediction | None:
    if HELIXFOLD_SERVICE is None:
        details = HELIXFOLD_LOAD_ERROR or "unknown initialization error"
        raise RuntimeError(
            "Mutated protein prediction was requested, but HelixFold-single failed to initialize: "
            f"{details}"
        )

    if not HELIXFOLD_SERVICE.available:
        details = HELIXFOLD_SERVICE.unavailable_reason or "unknown configuration error"
        raise RuntimeError(
            "Mutated protein prediction was requested, but HelixFold-single is not configured: "
            f"{details}"
        )

    mutated_translation = protein_sequence_for_prediction(str(analyser.mutated_amino_acid_translation))
    if not mutated_translation:
        raise RuntimeError(
            "The mutated sequence cannot be predicted because it does not produce a valid protein sequence."
        )

    return HELIXFOLD_SERVICE.predict_mutated_structure(
        organism=organism,
        gene=gene,
        gene_sequence=gene_sequence,
        mutations=report_mutations,
        protein_sequence=mutated_translation,
    )


def _generate_report_components(data: ReportData):
    context, result = describe_mutations(data.organism, data.gene, data.mutations)
    result.raise_for_errors()
    if not result.variants:
        raise ValueError("At least one variant is required")

    gene_sequence = context.cds_sequence
    analyser = Analyser(
        gene_sequence,
        transl_table=context.transl_table,
        codon_start=context.codon_start,
        transl_except=context.transl_except,
    )
    analyser.apply_mutations(_analyser_mutations(result.variants))
    report_mutations = [_report_mutation(variant) for variant in result.variants]
    report_text = reporter.generate_report(analyser, report_mutations)
    mutated_prediction = None

    if data.predict_mutated_structure:
        mutated_prediction = _predict_mutated_structure(
            organism=data.organism,
            gene=data.gene,
            gene_sequence=gene_sequence,
            report_mutations=report_mutations,
            analyser=analyser,
        )

    gene_context = {"context": context, "result": result}
    return report_text, report_mutations, analyser, mutated_prediction, gene_context


def _process_report_job(job_id: str):
    with REPORT_JOBS_LOCK:
        job = REPORT_JOBS.get(job_id)
        if job is None:
            return
        job["status"] = "running"
        job["updated_at"] = datetime.now(timezone.utc).isoformat()
        _persist_report_job_to_db(job)
        payload = ReportData(
            organism=job["organism"],
            gene=job["gene"],
            mutations=job["mutations"],
            predict_mutated_structure=job.get("predict_mutated_structure", False),
        )

    try:
        report_text, _, analyser, mutated_prediction, gene_context = _generate_report_components(payload)
        original_pdb_path = _build_original_pdb_path(payload.organism, payload.gene)
        mutated_pdb_path = mutated_prediction.mutated_pdb_path if mutated_prediction else None
        html_report = _build_report_html(
            payload.organism,
            payload.gene,
            report_text,
            gene_context["context"],
            gene_context["result"],
            str(analyser.amino_acid_translation),
            mutated_amino_acid_translation=str(analyser.mutated_amino_acid_translation or ""),
            genetic_code=analyser.genetic_code_description,
            original_pdb_path=original_pdb_path if original_pdb_path.exists() else None,
            mutated_pdb_path=mutated_pdb_path if mutated_pdb_path and mutated_pdb_path.exists() else None,
            mutated_prediction=mutated_prediction,
        )
        report_file_path = REPORTS_DIR / f"{job_id}.html"
        report_file_path.write_text(html_report, encoding="utf-8")

        with REPORT_JOBS_LOCK:
            completed_job = REPORT_JOBS.get(job_id)
            if completed_job is not None:
                completed_job["status"] = "completed"
                report_base_url = completed_job.get("report_base_url", "").rstrip("/")
                completed_job["report_url"] = f"{report_base_url}{API_PREFIX}/generated-reports/{job_id}.html"
                completed_job["updated_at"] = datetime.now(timezone.utc).isoformat()
                _persist_completed_job(job_id, completed_job, report_file_path)
    except Exception as exc:
        with REPORT_JOBS_LOCK:
            failed_job = REPORT_JOBS.get(job_id)
            if failed_job is not None:
                failed_job["status"] = "failed"
                failed_job["error"] = str(exc)
                failed_job["updated_at"] = datetime.now(timezone.utc).isoformat()


def _reference_payload(context: ReferenceContext) -> dict:
    return {
        "accession": context.accession,
        "gene": context.gene,
        "symbol": context.symbol,
        "proteinId": context.protein_id,
        "strand": context.strand,
        "codingNumberingAvailable": context.c_available,
        "cdsStartInGenome": context.first_genome_position,
        "cdsEndInGenome": context.last_genome_position,
    }


@router.post("/hgvs/describe", tags=["hgvs"])
def describe_hgvs(data: DescribeRequest):
    """Validate HGVS descriptions of one gene and return them normalised (m., c. and p.)."""
    try:
        context, result = describe_mutations(data.organism, data.gene, data.mutations)
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    return result.to_dict() | {"reference": _reference_payload(context)}


@router.post("/report/jobs", tags=["report"])
def create_report_job(data: ReportData, background_tasks: BackgroundTasks, request: Request):
    if not data.mutations:
        raise HTTPException(status_code=400, detail="At least one variant is required")

    try:
        _, result = describe_mutations(data.organism, data.gene, data.mutations)
        result.raise_for_errors()
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc

    job_id = uuid.uuid4().hex
    job_payload = _create_job_payload(job_id, data, str(request.base_url).rstrip("/"))
    with REPORT_JOBS_LOCK:
        REPORT_JOBS[job_id] = job_payload
    _persist_report_job_to_db(job_payload)

    background_tasks.add_task(_process_report_job, job_id)
    return _serialize_job(job_payload)


@router.get("/report/jobs", tags=["report"])
def list_report_jobs(request: Request):
    jobs = _collect_report_jobs(str(request.base_url).rstrip("/"))
    return {"jobs": jobs}


@router.get("/report/jobs/{job_id}", tags=["report"])
def get_report_job(job_id: str, request: Request):
    with REPORT_JOBS_LOCK:
        job = REPORT_JOBS.get(job_id)
    if job is None:
        db_job = _load_report_job_from_db(job_id, str(request.base_url).rstrip("/"))
        if db_job is not None:
            return db_job
        metadata_job: dict | None = _job_from_metadata(
            job_id,
            _report_metadata_path(job_id),
            str(request.base_url).rstrip("/"),
        )
        if metadata_job is None:
            metadata_job = _job_from_report_file(job_id, str(request.base_url).rstrip("/"))
        if metadata_job is None:
            raise HTTPException(status_code=404, detail="Report job not found")
        return metadata_job
    return _serialize_job(job)


@router.post("/report/", tags=["report"])
def generate_report(data: ReportData):
    try:
        report_text, _, _, _, _ = _generate_report_components(data)
    except LookupError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return report_text
