from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
from dataclasses import dataclass
from errno import ENOEXEC
from pathlib import Path
from threading import Lock
from typing import Any

from backend.analyser_backend.runtime_env import load_backend_env


load_backend_env()


@dataclass(slots=True, frozen=True)
class HelixFoldPrediction:
    cache_key: str
    cache_dir: Path
    mutated_pdb_path: Path
    protein_sequence: str
    metadata: dict[str, Any]


def normalize_mutations(mutations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []

    for mutation in mutations:
        normalized.append(
            {
                "type": str(mutation.get("type", "sub")),
                "position": int(mutation["position"]),
                "end_position": int(mutation.get("end_position", mutation["position"])),
                "ref": mutation.get("ref"),
                "alt": mutation.get("alt"),
            }
        )

    normalized.sort(
        key=lambda item: (
            item["position"],
            item["end_position"],
            item["type"],
            item.get("ref") or "",
            item.get("alt") or "",
        )
    )
    return normalized


def build_prediction_cache_key(
    *,
    organism: str,
    gene: str,
    gene_sequence: str,
    protein_sequence: str,
    mutations: list[dict[str, Any]],
    model_path: Path,
    script_path: Path,
) -> str:
    payload = {
        "organism": organism,
        "gene": gene,
        "gene_sequence": gene_sequence,
        "protein_sequence": protein_sequence,
        "mutations": normalize_mutations(mutations),
        "model_path": str(model_path),
        "script_path": str(script_path),
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    return digest


def protein_sequence_for_prediction(amino_acid_translation: str) -> str:
    return amino_acid_translation.split("*", 1)[0].strip()


def transform_pdb_coordinates(pdb_text: str, rotation: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]], translation: tuple[float, float, float]) -> str:
    transformed_lines: list[str] = []

    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            transformed_lines.append(line)
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        new_x = rotation[0][0] * x + rotation[0][1] * y + rotation[0][2] * z + translation[0]
        new_y = rotation[1][0] * x + rotation[1][1] * y + rotation[1][2] * z + translation[1]
        new_z = rotation[2][0] * x + rotation[2][1] * y + rotation[2][2] * z + translation[2]

        transformed_lines.append(
            f"{line[:30]}{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}{line[54:]}"
        )

    return "\n".join(transformed_lines) + ("\n" if pdb_text.endswith("\n") else "")


class HelixFoldService:
    def __init__(self, cache_root: Path | None = None) -> None:
        base_dir = cache_root or (Path(__file__).resolve().parents[1] / "pdb_files" / "cache" / "helixfold")
        self._cache_root = base_dir
        self._cache_root.mkdir(parents=True, exist_ok=True)
        self._lock = Lock()

        configured_enabled = os.getenv("HELIXFOLD_SINGLE_ENABLED", "").strip().lower()
        self._enabled = configured_enabled in {"1", "true", "yes", "on"}
        self._repo_dir = self._resolve_repo_dir()
        self._model_path = self._resolve_model_path()
        self._python_bin = self._resolve_python_bin()
        self._script_relpath = os.getenv("HELIXFOLD_SINGLE_SCRIPT_RELPATH", "helixfold_single_inference.py")
        self._timeout_seconds = int(os.getenv("HELIXFOLD_SINGLE_TIMEOUT_SECONDS", "7200"))

    @property
    def available(self) -> bool:
        self._refresh_configuration()
        return (
            self._enabled
            and self._repo_dir is not None
            and self._model_path is not None
            and self._python_bin is not None
            and self._resolve_script_path() is not None
        )

    @property
    def unavailable_reason(self) -> str | None:
        self._refresh_configuration()
        if self.available:
            return None
        return self._build_unavailable_reason()

    def predict_mutated_structure(
        self,
        *,
        organism: str,
        gene: str,
        gene_sequence: str,
        mutations: list[dict[str, Any]],
        protein_sequence: str,
    ) -> HelixFoldPrediction:
        self._refresh_configuration()
        if not self.available:
            raise RuntimeError(
                "HelixFold-single is disabled or not configured. "
                "Set HELIXFOLD_SINGLE_ENABLED=true and provide the repo/model paths."
            )

        assert self._repo_dir is not None
        assert self._model_path is not None
        assert self._python_bin is not None
        script_path = self._resolve_script_path()
        if script_path is None:
            raise RuntimeError(
                "HelixFold-single script not found. "
                "Check HELIXFOLD_SINGLE_REPO_DIR and HELIXFOLD_SINGLE_SCRIPT_RELPATH."
            )

        cache_key = build_prediction_cache_key(
            organism=organism,
            gene=gene,
            gene_sequence=gene_sequence,
            protein_sequence=protein_sequence,
            mutations=mutations,
            model_path=self._model_path,
            script_path=script_path,
        )

        cache_dir = self._cache_root / cache_key[:2] / cache_key
        cache_dir.mkdir(parents=True, exist_ok=True)
        metadata_path = cache_dir / "metadata.json"
        mutated_pdb_path = cache_dir / "mutated.pdb"

        with self._lock:
            if mutated_pdb_path.exists() and metadata_path.exists():
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                return HelixFoldPrediction(
                    cache_key=cache_key,
                    cache_dir=cache_dir,
                    mutated_pdb_path=mutated_pdb_path,
                    protein_sequence=protein_sequence,
                    metadata=metadata,
                )

            fasta_path = cache_dir / "input.fasta"
            fasta_path.write_text(f">{organism}_{gene}\n{protein_sequence}\n", encoding="utf-8")
            run_dir = cache_dir / "_run"
            run_dir.mkdir(parents=True, exist_ok=True)

            command = [
                self._python_bin,
                str(script_path),
                f"--init_model={self._model_path}",
                f"--fasta_file={fasta_path}",
                f"--output_dir={run_dir}",
            ]

            try:
                completed = subprocess.run(
                    command,
                    cwd=self._repo_dir,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=max(self._timeout_seconds, 60),
                )
            except OSError as exc:
                if exc.errno == ENOEXEC:
                    raise RuntimeError(
                        "HelixFold script points to a non-executable file. "
                        "Ensure HELIXFOLD_SINGLE_PYTHON_BIN references a working Python interpreter."
                    ) from exc
                raise RuntimeError(f"Failed to execute HelixFold-single: {exc}") from exc

            if completed.returncode != 0:
                stderr = completed.stderr.strip()
                stdout = completed.stdout.strip()
                details = stderr or stdout or "Unknown HelixFold-single error."
                raise RuntimeError(f"HelixFold-single failed: {self._with_runtime_hints(details)}")

            output_candidate = run_dir / "unrelaxed.pdb"
            if not output_candidate.exists():
                pdb_candidates = sorted(
                    run_dir.rglob("*.pdb"),
                    key=lambda path: path.stat().st_mtime,
                    reverse=True,
                )
                if not pdb_candidates:
                    raise RuntimeError(
                        "HelixFold-single completed but no PDB output was found in the output directory."
                    )
                output_candidate = pdb_candidates[0]

            shutil.copy2(output_candidate, mutated_pdb_path)
            metadata = {
                "engine": "helixfold_single",
                "repo_dir": str(self._repo_dir),
                "script": str(script_path),
                "model_path": str(self._model_path),
                "python_bin": self._python_bin,
                "timeout_seconds": self._timeout_seconds,
                "selected_output": output_candidate.name,
                "organism": organism,
                "gene": gene,
                "cache_key": cache_key,
                "protein_sequence_length": len(protein_sequence),
            }
            metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8")

            return HelixFoldPrediction(
                cache_key=cache_key,
                cache_dir=cache_dir,
                mutated_pdb_path=mutated_pdb_path,
                protein_sequence=protein_sequence,
                metadata=metadata,
            )

    def _resolve_repo_dir(self) -> Path | None:
        configured = os.getenv("HELIXFOLD_SINGLE_REPO_DIR", "").strip()
        if not configured:
            return None

        repo_dir = Path(configured).expanduser().resolve()
        if not repo_dir.exists() or not repo_dir.is_dir():
            return None

        return repo_dir

    def _resolve_model_path(self) -> Path | None:
        configured = os.getenv("HELIXFOLD_SINGLE_MODEL_PATH", "").strip()
        if not configured:
            return None

        model_path = Path(configured).expanduser().resolve()
        if not model_path.exists() or not model_path.is_file():
            return None

        return model_path

    def _resolve_script_path(self) -> Path | None:
        if self._repo_dir is None:
            return None

        script_path = (self._repo_dir / self._script_relpath).resolve()
        if not script_path.exists() or not script_path.is_file():
            return None

        return script_path

    def _resolve_python_bin(self) -> str | None:
        repo_python = self._repo_dir / ".venv" / "bin" / "python" if self._repo_dir else None
        if repo_python and repo_python.exists():
            return str(repo_python)

        configured = os.getenv("HELIXFOLD_SINGLE_PYTHON_BIN", "").strip()
        if configured:
            candidate = Path(configured).expanduser().resolve()
            if self._repo_dir is not None:
                expected_root = (self._repo_dir / ".venv").resolve()
                try:
                    candidate.relative_to(expected_root)
                except ValueError:
                    return None
            if candidate.exists():
                return str(candidate)

        return None

    def _with_runtime_hints(self, details: str) -> str:
        lowered = details.lower()
        hints: list[str] = []

        if "import paddle" in lowered or "module named 'paddle'" in lowered:
            hints.append(
                "Use a dedicated HelixFold environment and install paddlepaddle there, then set "
                "HELIXFOLD_SINGLE_PYTHON_BIN to that interpreter."
            )

        if "libnetcdf" in lowered or "pyharp" in lowered:
            hints.append(
                "This usually means an incompatible paddle stack in the HelixFold environment."
            )

        if "fused_gate_attention" in lowered:
            hints.append(
                "Your Paddle build does not expose the fused_gate_attention kernel expected by HelixFold."
            )

        if "dropout_nd" in lowered:
            hints.append(
                "Your Paddle build is missing the legacy dropout_nd op expected by HelixFold."
            )

        if not hints:
            return details

        return f"{details}\n\nHints:\n- " + "\n- ".join(hints)

    def _refresh_configuration(self) -> None:
        self._enabled = os.getenv("HELIXFOLD_SINGLE_ENABLED", "").strip().lower() in {"1", "true", "yes", "on"}
        self._repo_dir = self._resolve_repo_dir()
        self._model_path = self._resolve_model_path()
        self._python_bin = self._resolve_python_bin()
        self._script_relpath = os.getenv("HELIXFOLD_SINGLE_SCRIPT_RELPATH", "helixfold_single_inference.py")
        self._timeout_seconds = int(os.getenv("HELIXFOLD_SINGLE_TIMEOUT_SECONDS", "7200"))

    def _build_unavailable_reason(self) -> str:
        enabled_raw = os.getenv("HELIXFOLD_SINGLE_ENABLED", "").strip()
        if not self._enabled:
            return (
                "HELIXFOLD_SINGLE_ENABLED is not enabled "
                f"(current value: {enabled_raw!r})."
            )

        configured_repo_dir = os.getenv("HELIXFOLD_SINGLE_REPO_DIR", "").strip()
        if self._repo_dir is None:
            if not configured_repo_dir:
                return "HELIXFOLD_SINGLE_REPO_DIR is not set."
            return (
                "HELIXFOLD_SINGLE_REPO_DIR does not point to an existing directory: "
                f"{configured_repo_dir}"
            )

        configured_model_path = os.getenv("HELIXFOLD_SINGLE_MODEL_PATH", "").strip()
        if self._model_path is None:
            if not configured_model_path:
                return "HELIXFOLD_SINGLE_MODEL_PATH is not set."
            return (
                "HELIXFOLD_SINGLE_MODEL_PATH does not point to an existing file: "
                f"{configured_model_path}"
            )

        configured_python_bin = os.getenv("HELIXFOLD_SINGLE_PYTHON_BIN", "").strip()
        if self._python_bin is None:
            repo_python = self._repo_dir / ".venv" / "bin" / "python"
            if repo_python.exists():
                return (
                    "Unable to use Python from HelixFold repo virtual environment: "
                    f"{repo_python}"
                )
            if configured_python_bin:
                return (
                    "HELIXFOLD_SINGLE_PYTHON_BIN is invalid. It must exist and point "
                    "inside HELIXFOLD_SINGLE_REPO_DIR/.venv: "
                    f"{configured_python_bin}"
                )
            return (
                "Python interpreter for HelixFold is unavailable. Expected: "
                f"{repo_python}, or set HELIXFOLD_SINGLE_PYTHON_BIN "
                "to an interpreter inside that .venv."
            )

        script_path = self._resolve_script_path()
        if script_path is None:
            return (
                "HelixFold inference script not found at: "
                f"{(self._repo_dir / self._script_relpath).resolve()}"
            )

        return "Unknown HelixFold-single configuration error."
