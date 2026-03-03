"""
MST-Pipeline — File upload API endpoints.
"""
import os
import shutil
from pathlib import Path

from fastapi import APIRouter, Depends, File, UploadFile
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from app.config import UPLOAD_DIR, to_absolute, to_relative
from app.db.database import get_db
from app.db.models import FastqFile, Upload
from app.pipeline.detect import detect_sequencing_type, extract_sample_name

router = APIRouter(prefix="/api", tags=["upload"])


@router.post("/upload")
async def upload_files(files: list[UploadFile] = File(...), db: Session = Depends(get_db)):
    """Upload FASTQ files, detect sequencing type, create DB records."""
    if not files:
        return {"error": "No files provided"}

    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

    # Save files to flat UPLOAD_DIR
    filenames = []
    total_size = 0.0
    for f in files:
        dest = UPLOAD_DIR / f.filename
        with open(dest, "wb") as out:
            content = await f.read()
            out.write(content)
            total_size += len(content) / (1024 * 1024)
        filenames.append(f.filename)

    # Detect sequencing type
    detection = detect_sequencing_type(filenames)

    upload = Upload(
        upload_dir=str(UPLOAD_DIR),
        sequencing_type=detection["type"],
        total_files=len(filenames),
        total_size_mb=round(total_size, 2),
        status="uploaded",
    )
    db.add(upload)
    db.flush()

    # Create FastqFile records
    for sample_name, file_map in detection["samples"].items():
        for direction, filename in [("R1", file_map.get("R1")), ("R2", file_map.get("R2"))]:
            if not filename:
                continue
            abs_path = UPLOAD_DIR / filename
            file_size = os.path.getsize(abs_path) / (1024 * 1024) if abs_path.exists() else 0

            read_dir = direction
            if detection["type"] == "single-end":
                read_dir = "single"

            fastq = FastqFile(
                upload_id=upload.id,
                sample_name=sample_name,
                filename=filename,
                file_path=to_relative(abs_path),
                read_direction=read_dir,
                file_size_mb=round(file_size, 2),
            )
            db.add(fastq)

    db.commit()

    return {
        "upload_id": upload.id,
        "sequencing_type": detection["type"],
        "samples": len(detection["samples"]),
        "total_files": len(filenames),
        "errors": detection["errors"],
    }


@router.get("/uploads/{upload_id}/files/{filename}")
def download_file(upload_id: int, filename: str, db: Session = Depends(get_db)):
    """Download an individual FASTQ file from an upload."""
    fastq = (
        db.query(FastqFile)
        .filter(FastqFile.upload_id == upload_id, FastqFile.filename == filename)
        .first()
    )
    if not fastq:
        return {"error": "File not found"}
    file_path = to_absolute(fastq.file_path)
    if not file_path.exists():
        # Fallback: file may be in flat UPLOAD_DIR
        file_path = UPLOAD_DIR / filename
    if not file_path.exists():
        return {"error": "File not found on disk"}
    return FileResponse(str(file_path), filename=filename)


@router.get("/uploads")
def list_uploads(db: Session = Depends(get_db)):
    """List all uploads."""
    uploads = db.query(Upload).order_by(Upload.created_at.desc()).all()
    return [
        {
            "id": u.id,
            "sequencing_type": u.sequencing_type,
            "variable_region": u.variable_region,
            "total_files": u.total_files,
            "total_size_mb": u.total_size_mb,
            "status": u.status,
            "created_at": u.created_at.isoformat() if u.created_at else None,
        }
        for u in uploads
    ]


@router.get("/uploads/{upload_id}")
def get_upload(upload_id: int, db: Session = Depends(get_db)):
    """Get upload details with file list."""
    upload = db.query(Upload).filter(Upload.id == upload_id).first()
    if not upload:
        return {"error": "Upload not found"}

    files = db.query(FastqFile).filter(FastqFile.upload_id == upload_id).all()
    return {
        "id": upload.id,
        "upload_dir": upload.upload_dir,
        "sequencing_type": upload.sequencing_type,
        "variable_region": upload.variable_region,
        "total_files": upload.total_files,
        "total_size_mb": upload.total_size_mb,
        "status": upload.status,
        "created_at": upload.created_at.isoformat() if upload.created_at else None,
        "files": [
            {
                "id": f.id,
                "sample_name": f.sample_name,
                "filename": f.filename,
                "read_direction": f.read_direction,
                "file_size_mb": f.file_size_mb,
                "read_count": f.read_count,
                "avg_read_length": f.avg_read_length,
            }
            for f in files
        ],
    }
