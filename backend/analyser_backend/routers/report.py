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

from backend.analyser_backend.models.report_data import ReportData
from backend.analyser_cli import Reporter, Analyser, DataManager
from backend.analyser_cli.cds_annotation import DEFAULT_CODON_START, DEFAULT_TRANSL_TABLE
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

GENE_PREFIX = "c"
GENOME_PREFIX = "m"
SUPPORTED_PREFIXES = (GENE_PREFIX, GENOME_PREFIX)
_PREFIX = r"(?P<prefix>[cm])\."
SUBSTITUTION_PATTERN = re.compile(
    rf"(?i){_PREFIX}(?P<pos>\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])"
)
INSERTION_PATTERN = re.compile(
    rf"(?i){_PREFIX}(?P<left>\d+)_(?P<right>\d+)ins(?P<seq>[ACGT]+)"
)
DELETION_PATTERN = re.compile(
    rf"(?i){_PREFIX}(?P<start>\d+)(?:_(?P<end>\d+))?del(?:[ACGT]+)?"
)
DELINS_PATTERN = re.compile(
    rf"(?i){_PREFIX}(?P<start>\d+)(?:_(?P<end>\d+))?delins(?P<seq>[ACGT]+)"
)
DUPLICATION_PATTERN = re.compile(
    rf"(?i){_PREFIX}(?P<start>\d+)(?:_(?P<end>\d+))?dup(?:[ACGT]+)?"
)
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


def _unsupported_prefix_message(prefix: str | None) -> str:
    supported = " or ".join(f"'{value}.'" for value in SUPPORTED_PREFIXES)
    return (
        f"Unsupported coordinate prefix '{prefix}.'. "
        f"Only {supported} coordinates are supported "
        f"('{GENE_PREFIX}.' is gene-relative, '{GENOME_PREFIX}.' is genome-absolute)."
    )


def _unsupported_notation_message(mutation: str) -> str:
    return (
        f"Unsupported mutation notation: '{mutation}'. "
        f"Every mutation requires the '{GENE_PREFIX}.' (gene-relative) or "
        f"'{GENOME_PREFIX}.' (genome-absolute) prefix, e.g. 'm.5367C>T' or 'c.10_12del'."
    )


def _to_relative_position(position: int, prefix: str, start_genome_pos: int) -> int:
    """Convert an HGVS position to a 1-based position relative to the gene start."""
    if prefix == GENE_PREFIX:
        return position
    if prefix == GENOME_PREFIX:
        return position - start_genome_pos + 1
    raise ValueError(_unsupported_prefix_message(prefix))


def _ensure_within_gene(
    mutation: str, gene_sequence: str, start_genome_pos: int, *positions: int
) -> None:
    gene_length = len(gene_sequence)
    for position in positions:
        if position < 1 or position > gene_length:
            raise ValueError(
                f"Mutation '{mutation}' points to gene position {position}, which is outside "
                f"of the gene ({GENE_PREFIX}.1..{gene_length}, "
                f"{GENOME_PREFIX}.{start_genome_pos}..{start_genome_pos + gene_length - 1})."
            )


def _build_sub_mutation(pos: int, ref: str, alt: str, prefix: str, start_genome_pos: int):
    rel_pos = _to_relative_position(pos, prefix, start_genome_pos)
    return {
        "type": "sub",
        "position": rel_pos,
        "ref": ref,
        "alt": alt
    }


def _build_del_mutations(start: int, end: int, prefix: str, start_genome_pos: int):
    rel_start = _to_relative_position(start, prefix, start_genome_pos)
    rel_end = _to_relative_position(end, prefix, start_genome_pos)
    if rel_end < rel_start:
        rel_start, rel_end = rel_end, rel_start

    return [
        {
            "type": "del",
            "position": rel_pos,
            "ref": None,
            "alt": None
        }
        for rel_pos in range(rel_start, rel_end + 1)
    ]


def _build_del_mutation(start: int, end: int, prefix: str, start_genome_pos: int):
    rel_start = _to_relative_position(start, prefix, start_genome_pos)
    rel_end = _to_relative_position(end, prefix, start_genome_pos)
    if rel_end < rel_start:
        rel_start, rel_end = rel_end, rel_start

    return {
        "type": "del",
        "position": rel_start,
        "end_position": rel_end,
        "ref": None,
        "alt": None
    }


def _annotate_deletion(mutation: dict, deleted_length: int):
    length_change = -deleted_length
    mutation["deleted_length"] = deleted_length
    mutation["length_change"] = length_change
    mutation["frameshift"] = (deleted_length % 3) != 0
    return mutation


def _build_ins_mutation(pos: int, sequence: str, prefix: str, start_genome_pos: int):
    rel_pos = _to_relative_position(pos, prefix, start_genome_pos)
    return {
        "type": "ins",
        "position": rel_pos,
        "ref": None,
        "alt": sequence,
        "inserted_length": len(sequence),
        "length_change": len(sequence)
    }


def _match_mutation(mutation: str):
    """Return (kind, match) for a supported HGVS notation or raise ValueError."""
    candidate = mutation.strip()
    for kind, pattern in (
        ("sub", SUBSTITUTION_PATTERN),
        ("delins", DELINS_PATTERN),
        ("ins", INSERTION_PATTERN),
        ("del", DELETION_PATTERN),
        ("dup", DUPLICATION_PATTERN),
    ):
        matched = pattern.fullmatch(candidate)
        if matched:
            return kind, matched

    unsupported_prefix = re.fullmatch(r"(?i)\s*([a-z])\.\S*", candidate)
    if unsupported_prefix:
        prefix = unsupported_prefix.group(1).lower()
        if prefix not in SUPPORTED_PREFIXES:
            raise ValueError(_unsupported_prefix_message(prefix))
    raise ValueError(_unsupported_notation_message(candidate))


def validate_mutation_notation(mutation: str) -> None:
    """Raise ValueError when the notation is not a supported c./m. mutation."""
    _match_mutation(mutation)


def _parse_mutation(mutation: str, start_genome_pos: int, gene_sequence: str):
    report_mutations = []
    analyser_mutations = []
    mutation = mutation.strip()
    kind, matched = _match_mutation(mutation)
    prefix = matched.group("prefix").lower()

    if kind == "sub":
        pos = int(matched.group("pos"))
        _ensure_within_gene(
            mutation,
            gene_sequence,
            start_genome_pos,
            _to_relative_position(pos, prefix, start_genome_pos),
        )
        parsed_sub_mutation = _build_sub_mutation(
            pos,
            matched.group("ref").upper(),
            matched.group("alt").upper(),
            prefix,
            start_genome_pos,
        )
        report_mutations.append(parsed_sub_mutation)
        analyser_mutations.append(parsed_sub_mutation)
        return report_mutations, analyser_mutations

    if kind == "delins":
        start = int(matched.group("start"))
        end = int(matched.group("end")) if matched.group("end") else start
        insertion_anchor = max(start, end)
        sequence = matched.group("seq").upper()
        rel_start = _to_relative_position(start, prefix, start_genome_pos)
        rel_end = _to_relative_position(end, prefix, start_genome_pos)
        if rel_end < rel_start:
            rel_start, rel_end = rel_end, rel_start
        _ensure_within_gene(mutation, gene_sequence, start_genome_pos, rel_start, rel_end)

        del_mutation = _build_del_mutation(start, end, prefix, start_genome_pos)
        deleted_length = rel_end - rel_start + 1
        inserted_length = len(sequence)
        length_change = inserted_length - deleted_length
        report_mutations.append(
            del_mutation | {
                "type": "delins",
                "alt": sequence,
                "deleted_length": deleted_length,
                "inserted_length": inserted_length,
                "length_change": length_change,
                "frameshift": (length_change % 3) != 0
            }
        )
        analyser_mutations.extend(_build_del_mutations(start, end, prefix, start_genome_pos))
        analyser_mutations.append(
            _build_ins_mutation(insertion_anchor, sequence, prefix, start_genome_pos)
        )
        return report_mutations, analyser_mutations

    if kind == "ins":
        left = int(matched.group("left"))
        _ensure_within_gene(
            mutation,
            gene_sequence,
            start_genome_pos,
            _to_relative_position(left, prefix, start_genome_pos),
        )
        ins_mutation = _build_ins_mutation(
            left,
            matched.group("seq").upper(),
            prefix,
            start_genome_pos,
        )
        report_mutations.append(ins_mutation)
        analyser_mutations.append(ins_mutation)
        return report_mutations, analyser_mutations

    if kind == "del":
        start = int(matched.group("start"))
        end = int(matched.group("end")) if matched.group("end") else start
        _ensure_within_gene(
            mutation,
            gene_sequence,
            start_genome_pos,
            _to_relative_position(start, prefix, start_genome_pos),
            _to_relative_position(end, prefix, start_genome_pos),
        )
        deleted_length = abs(end - start) + 1
        report_mutations.append(
            _annotate_deletion(_build_del_mutation(start, end, prefix, start_genome_pos), deleted_length)
        )
        analyser_mutations.extend(_build_del_mutations(start, end, prefix, start_genome_pos))
        return report_mutations, analyser_mutations

    start = int(matched.group("start"))
    end = int(matched.group("end")) if matched.group("end") else start
    insertion_anchor = max(start, end)

    rel_start = _to_relative_position(start, prefix, start_genome_pos)
    rel_end = _to_relative_position(end, prefix, start_genome_pos)
    if rel_end < rel_start:
        rel_start, rel_end = rel_end, rel_start
    _ensure_within_gene(mutation, gene_sequence, start_genome_pos, rel_start, rel_end)

    duplicated_sequence = gene_sequence[rel_start - 1:rel_end]
    duplicated_length = rel_end - rel_start + 1
    report_mutations.append({
        "type": "dup",
        "position": rel_start,
        "end_position": rel_end,
        "ref": None,
        "alt": duplicated_sequence,
        "duplicated_length": duplicated_length,
        "length_change": duplicated_length,
        "frameshift": (duplicated_length % 3) != 0
    })
    analyser_mutations.append(
        _build_ins_mutation(insertion_anchor, duplicated_sequence, prefix, start_genome_pos)
    )
    return report_mutations, analyser_mutations


def _build_report_notation(mutation: dict, start_genome_pos: int):
    """Build the genome-absolute (m.) HGVS notation for a gene-relative mutation."""
    mut_type = mutation.get("type", "sub")
    rel_pos = _get_required_position(mutation)
    rel_end_pos = int(mutation.get("end_position", rel_pos))
    pos = start_genome_pos + rel_pos - 1
    end_pos = start_genome_pos + rel_end_pos - 1
    locus = f"{pos}" if pos == end_pos else f"{pos}_{end_pos}"
    if mut_type == "sub":
        return f"{GENOME_PREFIX}.{pos}{mutation.get('ref', '?')}>{mutation.get('alt', '?')}"
    if mut_type == "del":
        return f"{GENOME_PREFIX}.{locus}del"
    if mut_type == "ins":
        return f"{GENOME_PREFIX}.{pos}_{pos + 1}ins{mutation.get('alt', '')}"
    if mut_type == "delins":
        return f"{GENOME_PREFIX}.{locus}delins{mutation.get('alt', '')}"
    if mut_type == "dup":
        return f"{GENOME_PREFIX}.{locus}dup"
    return f"{GENOME_PREFIX}.{pos}?"


def _get_required_position(mutation: dict) -> int:
    pos = mutation.get("position")
    if pos is None:
        raise ValueError("Mutation position is required")
    return int(pos)


def _format_position_range(start: int, end: int | None = None) -> str:
    final_end = start if end is None else end
    return str(start) if start == final_end else f"{start}_{final_end}"


def _build_genome_position_display(mutation: dict, start_genome_pos: int) -> str:
    pos = _get_required_position(mutation)
    end_pos = int(mutation.get("end_position", pos))
    genome_start = start_genome_pos + pos - 1
    genome_end = start_genome_pos + end_pos - 1
    mut_type = mutation.get("type", "sub")

    if mut_type == "ins":
        return f"{genome_start}_{genome_start + 1}"

    return _format_position_range(genome_start, genome_end)


def _build_gene_position_display(mutation: dict) -> str:
    pos = _get_required_position(mutation)
    if mutation.get("type", "sub") == "ins":
        return f"{pos}_{pos + 1}"

    return _format_position_range(pos, mutation.get("end_position"))


def _build_report_consequence(mutation: dict):
    mut_type = mutation.get("type", "sub")
    if mut_type == "sub":
        return "Substitution"
    if mutation.get("frameshift") is True:
        return "Frameshift mutation"
    if mut_type in {"del", "ins", "delins", "dup"}:
        return "In-frame indel"
    return mut_type.upper()


def _build_static_url(path: Path) -> str:
    relative_path = path.relative_to(PDB_FILES_DIR).as_posix()
    return f"/pdb-files/{relative_path}"


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

    blocks = []
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
            "</div></div>"
        )

    return f'<div class="alignment">{"".join(blocks)}</div>'


def _build_report_html(
    organism: str,
    gene: str,
    report_text: str,
    report_mutations: list[dict],
    amino_acid_translation: str,
    *,
    mutated_amino_acid_translation: str,
    start_genome_pos: int,
    end_genome_pos: int,
    genetic_code: str,
    original_pdb_path: Path | None,
    mutated_pdb_path: Path | None,
    mutated_prediction: HelixFoldPrediction | None,
):
    total_mutations = len(report_mutations)
    frameshift_mutations = sum(1 for mutation in report_mutations if mutation.get("frameshift") is True)
    inframe_mutations = total_mutations - frameshift_mutations
    mutation_type_counts: dict[str, int] = {}
    for mutation in report_mutations:
        mutation_type = mutation.get("type", "unknown").upper()
        mutation_type_counts[mutation_type] = mutation_type_counts.get(mutation_type, 0) + 1

    original_pdb_url = _build_static_url(original_pdb_path) if original_pdb_path else None
    mutated_pdb_url = _build_static_url(mutated_pdb_path) if mutated_pdb_path else None
    prediction_available = mutated_pdb_url is not None
    alignment = align_translations(amino_acid_translation, mutated_amino_acid_translation)

    mutation_rows = []
    for mutation in report_mutations:
        gene_position = _build_gene_position_display(mutation)
        genome_position = _build_genome_position_display(mutation, start_genome_pos)
        mutation_rows.append(
            f"<tr><td>{html.escape(_build_report_notation(mutation, start_genome_pos))}</td>"
            f"<td>{html.escape(gene_position)}</td>"
            f"<td>{html.escape(genome_position)}</td>"
            f"<td>{html.escape(mutation.get('type', '').upper())}</td>"
            f"<td>{html.escape(_build_report_consequence(mutation))}</td></tr>"
        )
    mutation_table_rows = "\n".join(mutation_rows) or "<tr><td colspan='5'>No mutations</td></tr>"

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

    mutation_type_items = "".join(
        f"<li><span>{html.escape(key)}</span><strong>{count}</strong></li>"
        for key, count in sorted(mutation_type_counts.items())
    )
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
            "GENE": html.escape(gene),
            "GENOME_RANGE": html.escape(f"{start_genome_pos}..{end_genome_pos}"),
            "GENETIC_CODE": html.escape(genetic_code),
            "GENERATED_AT": html.escape(datetime.now(timezone.utc).isoformat()),
            "TOTAL_MUTATIONS": str(total_mutations),
            "FRAMESHIFT_MUTATIONS": str(frameshift_mutations),
            "INFRAME_MUTATIONS": str(inframe_mutations),
            "ORIGINAL_PDB_STATUS": "available" if original_pdb_url else "missing",
            "MUTATED_PDB_STATUS": "available" if prediction_available else "missing",
            "MUTATION_TYPE_ITEMS": mutation_type_items,
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
            "MUTATION_TABLE_ROWS": mutation_table_rows,
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
            data["report_url"] = f"{report_base_url}/generated-reports/{job_id}.html"
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
            data["report_url"] = f"{report_base_url}/generated-reports/{row['job_id']}.html"
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
        "report_url": f"{report_base_url}/generated-reports/{job_id}.html",
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
            data["report_url"] = f"{report_base_url}/generated-reports/{job_id}.html"
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
        data["report_url"] = f"{report_base_url}/generated-reports/{job_id}.html"
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
    gene_data = manager.get_gene_for_organism(data.organism, data.gene)
    gene_sequence = gene_data["sequence"]
    start_genome_pos = gene_data.get("startInGenome", 1)
    end_genome_pos = gene_data.get("endInGenome", start_genome_pos + len(gene_sequence) - 1)

    analyser = Analyser(
        gene_sequence,
        transl_table=gene_data.get("translTable", DEFAULT_TRANSL_TABLE),
        codon_start=gene_data.get("codonStart", DEFAULT_CODON_START),
        transl_except=gene_data.get("translExcept"),
    )

    report_mutations = []
    parsed_mutations = []
    for mutation in data.mutations:
        report_items, analyser_items = _parse_mutation(mutation, start_genome_pos, gene_sequence)
        report_mutations.extend(report_items)
        parsed_mutations.extend(analyser_items)

    analyser.apply_mutations(parsed_mutations)
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


    gene_context = {
        "start_genome_pos": start_genome_pos,
        "end_genome_pos": end_genome_pos,
    }
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
        report_text, report_mutations, analyser, mutated_prediction, gene_context = _generate_report_components(payload)
        original_pdb_path = _build_original_pdb_path(payload.organism, payload.gene)
        mutated_pdb_path = mutated_prediction.mutated_pdb_path if mutated_prediction else None
        html_report = _build_report_html(
            payload.organism,
            payload.gene,
            report_text,
            report_mutations,
            str(analyser.amino_acid_translation),
            mutated_amino_acid_translation=str(analyser.mutated_amino_acid_translation or ""),
            start_genome_pos=gene_context["start_genome_pos"],
            end_genome_pos=gene_context["end_genome_pos"],
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
                completed_job["report_url"] = f"{report_base_url}/generated-reports/{job_id}.html"
                completed_job["updated_at"] = datetime.now(timezone.utc).isoformat()
                _persist_completed_job(job_id, completed_job, report_file_path)
    except Exception as exc:
        with REPORT_JOBS_LOCK:
            failed_job = REPORT_JOBS.get(job_id)
            if failed_job is not None:
                failed_job["status"] = "failed"
                failed_job["error"] = str(exc)
                failed_job["updated_at"] = datetime.now(timezone.utc).isoformat()


@router.post("/report/jobs", tags=["report"])
def create_report_job(data: ReportData, background_tasks: BackgroundTasks, request: Request):
    if not data.mutations:
        raise HTTPException(status_code=400, detail="At least one mutation is required")

    for mutation in data.mutations:
        try:
            validate_mutation_notation(mutation)
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
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return report_text
