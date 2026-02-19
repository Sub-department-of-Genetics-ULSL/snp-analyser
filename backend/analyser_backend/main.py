from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware
from routers import organisms_router, report_router


app = FastAPI()
app.include_router(organisms_router)
app.include_router(report_router)

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