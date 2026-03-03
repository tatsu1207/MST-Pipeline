"""
MST-Pipeline — Pipeline control API endpoints.
"""
from fastapi import APIRouter, Depends
from pydantic import BaseModel
from sqlalchemy.orm import Session

from app.config import to_absolute
from app.db.database import get_db
from app.db.models import Dataset, Sample
from app.pipeline.runner import (
    cancel_pipeline,
    get_pathogen_status,
    get_pipeline_status,
    get_sourcetracker_status,
    launch_dada2_pipeline,
    launch_pathogen_run,
    launch_sourcetracker_run,
)

router = APIRouter(prefix="/api", tags=["pipeline"])


# -- DADA2 --------------------------------------------------------------------

class Dada2LaunchRequest(BaseModel):
    upload_id: int
    name: str = "DADA2 Run"
    threads: int | None = None
    error_model: str = "default"


@router.post("/dada2/launch")
def launch_dada2(req: Dada2LaunchRequest, db: Session = Depends(get_db)):
    """Launch a DADA2 pipeline. Truncation lengths are auto-detected."""
    dataset_id = launch_dada2_pipeline(
        upload_id=req.upload_id,
        name=req.name,
        threads=req.threads,
        error_model=req.error_model,
    )
    return {"dataset_id": dataset_id, "status": "launched"}


@router.post("/dada2/cancel/{dataset_id}")
def cancel_dada2(dataset_id: int):
    """Cancel a running DADA2 pipeline."""
    cancelled = cancel_pipeline(dataset_id)
    return {"cancelled": cancelled}


@router.get("/dada2/status/{dataset_id}")
def dada2_status(dataset_id: int):
    """Get DADA2 pipeline status."""
    return get_pipeline_status(dataset_id)


@router.get("/datasets")
def list_datasets(db: Session = Depends(get_db)):
    """List all completed datasets."""
    datasets = db.query(Dataset).order_by(Dataset.created_at.desc()).all()
    return [
        {
            "id": d.id,
            "name": d.name,
            "status": d.status,
            "sequencing_type": d.sequencing_type,
            "variable_region": d.variable_region,
            "sample_count": d.sample_count,
            "asv_count": d.asv_count,
            "asv_table_path": d.asv_table_path,
            "created_at": d.created_at.isoformat() if d.created_at else None,
        }
        for d in datasets
    ]


# -- SourceTracker -------------------------------------------------------------

class SourcetrackerLaunchRequest(BaseModel):
    dataset_id: int
    sample_names: list[str] | None = None
    selected_groups: list[str]
    feature_mode: str = "asv"
    src_depth: int = 1000
    snk_depth: int = 1000
    restarts: int = 10
    burnin: int = 100
    draws: int = 1
    name: str = "SourceTracker Run"


@router.post("/sourcetracker/launch")
def launch_st(req: SourcetrackerLaunchRequest, db: Session = Depends(get_db)):
    """Launch a SourceTracker run."""
    query = db.query(Sample).filter(Sample.dataset_id == req.dataset_id, Sample.asv_table_path.isnot(None))
    if req.sample_names:
        query = query.filter(Sample.sample_name.in_(req.sample_names))
    samples = query.all()
    paths = [str(to_absolute(s.asv_table_path)) for s in samples]
    if not paths:
        return {"error": "No per-sample ASV tables found for the given dataset/samples"}

    run_id = launch_sourcetracker_run(
        sample_table_paths=paths,
        selected_groups=req.selected_groups,
        feature_mode=req.feature_mode,
        src_depth=req.src_depth,
        snk_depth=req.snk_depth,
        restarts=req.restarts,
        burnin=req.burnin,
        draws=req.draws,
        dataset_id=req.dataset_id,
        name=req.name,
    )
    return {"run_id": run_id, "status": "launched"}


@router.get("/sourcetracker/status/{run_id}")
def st_status(run_id: int):
    """Get SourceTracker run status."""
    return get_sourcetracker_status(run_id)


@router.get("/source-groups")
def source_groups():
    """Get available source groups from MST.design."""
    from app.pipeline.sourcetracker import load_design, GROUP_COLORS
    design = load_design()
    groups = sorted(set(design.values()))
    return {
        "groups": groups,
        "colors": {g: GROUP_COLORS.get(g, "#888888") for g in groups},
        "sample_counts": {
            g: sum(1 for v in design.values() if v == g) for g in groups
        },
    }


# -- Pathogen ------------------------------------------------------------------

class PathogenLaunchRequest(BaseModel):
    dataset_id: int
    sample_names: list[str] | None = None
    silva_identity: float = 0.97
    name: str = "Pathogen Run"


@router.post("/pathogen/launch")
def launch_path(req: PathogenLaunchRequest, db: Session = Depends(get_db)):
    """Launch a pathogen detection run."""
    query = db.query(Sample).filter(Sample.dataset_id == req.dataset_id, Sample.asv_table_path.isnot(None))
    if req.sample_names:
        query = query.filter(Sample.sample_name.in_(req.sample_names))
    samples = query.all()
    paths = [str(to_absolute(s.asv_table_path)) for s in samples]
    if not paths:
        return {"error": "No per-sample ASV tables found for the given dataset/samples"}

    run_id = launch_pathogen_run(
        sample_table_paths=paths,
        silva_identity=req.silva_identity,
        dataset_id=req.dataset_id,
        name=req.name,
    )
    return {"run_id": run_id, "status": "launched"}


@router.get("/pathogen/status/{run_id}")
def path_status(run_id: int):
    """Get pathogen detection run status."""
    return get_pathogen_status(run_id)
