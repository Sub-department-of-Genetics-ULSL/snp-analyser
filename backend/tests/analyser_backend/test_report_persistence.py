import json
from pathlib import Path

from backend.analyser_backend.routers import report


def test_collect_report_jobs_includes_legacy_html_reports(tmp_path, monkeypatch):
    reports_dir = tmp_path / "generated_reports"
    reports_dir.mkdir()
    html_path = reports_dir / "legacy-job.html"
    html_path.write_text(
        """
        <html>
          <body>
            <section class="card">
              <h2>Metadata</h2>
              <div><strong>Organism:</strong> human</div>
              <div><strong>Gene:</strong> CYTB</div>
              <div><strong>Generated:</strong> 2026-07-07T14:00:00+00:00</div>
            </section>
            <table>
              <tbody>
                <tr><td>m.14747G>A</td><td>SUB</td><td>Substitution</td></tr>
              </tbody>
            </table>
            <li><span>Mutated PDB</span><strong>available</strong></li>
          </body>
        </html>
        """,
        encoding="utf-8",
    )

    monkeypatch.setattr(report, "REPORTS_DIR", reports_dir)
    monkeypatch.setattr(report, "REPORT_JOBS", {})

    jobs = report._collect_report_jobs("http://127.0.0.1:8000")

    assert len(jobs) == 1
    assert jobs[0]["id"] == "legacy-job"
    assert jobs[0]["status"] == "completed"
    assert jobs[0]["organism"] == "human"
    assert jobs[0]["gene"] == "CYTB"
    assert jobs[0]["mutations"] == ["m.14747G>A"]
    assert jobs[0]["report_url"] == "http://127.0.0.1:8000/api/generated-reports/legacy-job.html"
    assert report._load_report_job_from_db("legacy-job", "http://127.0.0.1:8000") is not None


def test_collect_report_jobs_prefers_metadata_sidecar(tmp_path, monkeypatch):
    reports_dir = tmp_path / "generated_reports"
    reports_dir.mkdir()
    html_path = reports_dir / "job-123.html"
    html_path.write_text("<html>report</html>", encoding="utf-8")
    metadata_path = reports_dir / "job-123.json"
    metadata_path.write_text(
        json.dumps(
            {
                "id": "job-123",
                "status": "completed",
                "organism": "human",
                "gene": "CYTB",
                "mutations": ["14747G>A"],
                "predict_mutated_structure": True,
                "created_at": "2026-07-07T14:00:00+00:00",
                "updated_at": "2026-07-07T14:01:00+00:00",
                "report_url": None,
                "report_base_url": "http://old-host",
                "error": None,
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(report, "REPORTS_DIR", reports_dir)
    monkeypatch.setattr(report, "REPORT_JOBS", {})

    jobs = report._collect_report_jobs("http://127.0.0.1:8000")

    assert len(jobs) == 1
    assert jobs[0]["id"] == "job-123"
    assert jobs[0]["organism"] == "human"
    assert jobs[0]["report_url"] == "http://127.0.0.1:8000/api/generated-reports/job-123.html"


def test_collect_report_jobs_loads_sqlite_backed_history(tmp_path, monkeypatch):
    reports_dir = tmp_path / "generated_reports"
    reports_dir.mkdir()

    monkeypatch.setattr(report, "REPORTS_DIR", reports_dir)
    monkeypatch.setattr(report, "REPORT_JOBS", {})

    report._persist_report_job_to_db(
        {
            "id": "db-job",
            "status": "completed",
            "organism": "human",
            "gene": "CYTB",
            "mutations": ["m.14747G>A"],
            "predict_mutated_structure": True,
            "created_at": "2026-07-07T14:00:00+00:00",
            "updated_at": "2026-07-07T14:01:00+00:00",
            "report_url": "http://old-host/generated-reports/db-job.html",
            "report_base_url": "http://old-host",
            "error": None,
        },
        "db-job.html",
    )

    jobs = report._collect_report_jobs("http://127.0.0.1:8000")

    assert len(jobs) == 1
    assert jobs[0]["id"] == "db-job"
    assert jobs[0]["organism"] == "human"
    assert jobs[0]["gene"] == "CYTB"
    assert jobs[0]["mutations"] == ["m.14747G>A"]
    assert jobs[0]["report_url"] == "http://127.0.0.1:8000/api/generated-reports/db-job.html"
