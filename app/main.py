"""
MST-Pipeline — Application entry point.

Wires FastAPI (API endpoints) + Plotly Dash (UI) into a single ASGI app.
Start with: uvicorn app.main:app --reload --host 0.0.0.0 --port 8050
"""
from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.middleware.wsgi import WSGIMiddleware

from app.api.upload import router as upload_router
from app.api.pipeline import router as pipeline_router
from app.config import DATA_DIR
from app.dashboard.app import app as dash_app
from app.dashboard.layout import create_layout
from app.db.database import init_db

# Eagerly import pages so their @callback decorators register with Dash
import app.dashboard.pages.unified_page  # noqa: F401


# -- FastAPI application -------------------------------------------------------

api = FastAPI(
    title="MST-Pipeline API",
    docs_url="/api/docs",
    openapi_url="/api/openapi.json",
)

# Register API routers BEFORE mounting Dash (order matters)
api.include_router(upload_router)
api.include_router(pipeline_router)

# Serve dataset output files
api.mount("/static-data", StaticFiles(directory=str(DATA_DIR)), name="static-data")


@api.on_event("startup")
def on_startup():
    """Initialize database and data directories on server start."""
    init_db()


# -- Mount Dash inside FastAPI -------------------------------------------------

# Set the Dash layout (done here to avoid circular imports)
dash_app.layout = create_layout()

# Mount the Dash WSGI app as a catch-all under "/"
# All requests not matching /api/* routes fall through to Dash
api.mount("/", WSGIMiddleware(dash_app.server))

# This is what uvicorn imports: app.main:app
app = api
