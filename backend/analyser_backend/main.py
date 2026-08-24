from pathlib import Path
import logging

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.middleware.cors import CORSMiddleware
from routers import organisms_router, report_router
from routers.report import get_helixfold_startup_status


app = FastAPI()
app.include_router(organisms_router)
app.include_router(report_router)
logger = logging.getLogger(__name__)

origins = [
    "http://localhost.tiangolo.com",
    "https://localhost.tiangolo.com",
    "http://localhost",
    "http://localhost:8000",
    "http://127.0.0.1:8000"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app_root = Path(__file__).resolve().parent
generated_reports_dir = app_root / "generated_reports"
pdb_files_dir = app_root / "pdb_files"
generated_reports_dir.mkdir(parents=True, exist_ok=True)
pdb_files_dir.mkdir(parents=True, exist_ok=True)

app.mount(
    "/generated-reports",
    StaticFiles(directory=str(generated_reports_dir)),
    name="generated-reports",
)
app.mount(
    "/pdb-files",
    StaticFiles(directory=str(pdb_files_dir)),
    name="pdb-files",
)


@app.on_event("startup")
def _log_helixfold_startup_status() -> None:
    available, reason = get_helixfold_startup_status()
    if not available:
        logger.warning("HelixFold unavailable at startup: %s", reason)
