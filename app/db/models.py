"""
MST-Pipeline — SQLAlchemy 2.0 ORM Models
"""
from datetime import datetime

from sqlalchemy import (
    Boolean,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


class Base(DeclarativeBase):
    pass


# -- Uploads -------------------------------------------------------------------

class Upload(Base):
    __tablename__ = "uploads"

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_dir: Mapped[str] = mapped_column(String, nullable=False)
    sequencing_type: Mapped[str | None] = mapped_column(String)
    variable_region: Mapped[str | None] = mapped_column(String)
    primers_detected: Mapped[bool | None] = mapped_column(Boolean)
    study: Mapped[str | None] = mapped_column(String)
    total_files: Mapped[int | None] = mapped_column(Integer)
    total_size_mb: Mapped[float | None] = mapped_column(Float)
    status: Mapped[str] = mapped_column(String, default="uploaded")
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    fastq_files: Mapped[list["FastqFile"]] = relationship(
        back_populates="upload", cascade="all, delete-orphan"
    )


# -- Upload Metadata -----------------------------------------------------------

class UploadMetadata(Base):
    __tablename__ = "upload_metadata"
    __table_args__ = (
        UniqueConstraint("upload_id", "sample_name", "key", name="uq_upload_sample_key"),
        Index("idx_meta_upload", "upload_id"),
    )

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_id: Mapped[int] = mapped_column(
        ForeignKey("uploads.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    key: Mapped[str] = mapped_column(String, nullable=False)
    value: Mapped[str | None] = mapped_column(String)

    upload: Mapped["Upload"] = relationship()


# -- FASTQ Files ---------------------------------------------------------------

class FastqFile(Base):
    __tablename__ = "fastq_files"
    __table_args__ = (Index("idx_fastq_upload", "upload_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_id: Mapped[int] = mapped_column(
        ForeignKey("uploads.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    filename: Mapped[str] = mapped_column(String, nullable=False)
    file_path: Mapped[str] = mapped_column(String, nullable=False)
    read_direction: Mapped[str | None] = mapped_column(String)  # R1/R2/single
    file_size_mb: Mapped[float | None] = mapped_column(Float)
    read_count: Mapped[int | None] = mapped_column(Integer)
    avg_read_length: Mapped[int | None] = mapped_column(Integer)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    upload: Mapped["Upload"] = relationship(back_populates="fastq_files")


# -- Datasets (DADA2 run output) -----------------------------------------------

class Dataset(Base):
    __tablename__ = "datasets"

    id: Mapped[int] = mapped_column(primary_key=True)
    upload_id: Mapped[int | None] = mapped_column(ForeignKey("uploads.id"))
    name: Mapped[str] = mapped_column(String, nullable=False)
    status: Mapped[str] = mapped_column(String, default="pending")
    sequencing_type: Mapped[str | None] = mapped_column(String)
    variable_region: Mapped[str | None] = mapped_column(String)
    sample_count: Mapped[int | None] = mapped_column(Integer)
    asv_count: Mapped[int | None] = mapped_column(Integer)
    asv_table_path: Mapped[str | None] = mapped_column(String)  # BIOM HDF5
    rep_seqs_path: Mapped[str | None] = mapped_column(String)
    pipeline_log_path: Mapped[str | None] = mapped_column(String)
    trunc_len_f: Mapped[int | None] = mapped_column(Integer)
    trunc_len_r: Mapped[int | None] = mapped_column(Integer)
    error_model: Mapped[str | None] = mapped_column(String)  # default/PacBio/Nanopore
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow, onupdate=datetime.utcnow
    )

    upload: Mapped["Upload | None"] = relationship()
    samples: Mapped[list["Sample"]] = relationship(
        back_populates="dataset", cascade="all, delete-orphan"
    )


# -- Samples -------------------------------------------------------------------

class Sample(Base):
    __tablename__ = "samples"
    __table_args__ = (Index("idx_samples_dataset", "dataset_id"),)

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int] = mapped_column(
        ForeignKey("datasets.id", ondelete="CASCADE")
    )
    sample_name: Mapped[str] = mapped_column(String, nullable=False)
    read_count_raw: Mapped[int | None] = mapped_column(Integer)
    read_count_filtered: Mapped[int | None] = mapped_column(Integer)
    read_count_nonchimeric: Mapped[int | None] = mapped_column(Integer)
    asv_count: Mapped[int | None] = mapped_column(Integer)
    asv_table_path: Mapped[str | None] = mapped_column(String)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    dataset: Mapped["Dataset"] = relationship(back_populates="samples")


# -- SourceTracker Runs --------------------------------------------------------

class SourcetrackerRun(Base):
    __tablename__ = "sourcetracker_runs"

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int | None] = mapped_column(ForeignKey("datasets.id"))
    name: Mapped[str] = mapped_column(String, nullable=False)
    status: Mapped[str] = mapped_column(String, default="pending")
    input_biom_path: Mapped[str | None] = mapped_column(String)
    selected_groups: Mapped[str | None] = mapped_column(Text)  # JSON list
    feature_mode: Mapped[str | None] = mapped_column(String)  # asv/otu
    src_depth: Mapped[int | None] = mapped_column(Integer)
    snk_depth: Mapped[int | None] = mapped_column(Integer)
    restarts: Mapped[int | None] = mapped_column(Integer)
    burnin: Mapped[int | None] = mapped_column(Integer)
    draws: Mapped[int | None] = mapped_column(Integer)
    proportions_path: Mapped[str | None] = mapped_column(String)
    stds_path: Mapped[str | None] = mapped_column(String)
    figure_path: Mapped[str | None] = mapped_column(String)
    output_dir: Mapped[str | None] = mapped_column(String)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    dataset: Mapped["Dataset | None"] = relationship()


# -- Pathogen Runs -------------------------------------------------------------

class PathogenRun(Base):
    __tablename__ = "pathogen_runs"

    id: Mapped[int] = mapped_column(primary_key=True)
    dataset_id: Mapped[int | None] = mapped_column(ForeignKey("datasets.id"))
    name: Mapped[str] = mapped_column(String, nullable=False)
    status: Mapped[str] = mapped_column(String, default="pending")
    input_biom_path: Mapped[str | None] = mapped_column(String)
    silva_identity: Mapped[float | None] = mapped_column(Float)
    output_dir: Mapped[str | None] = mapped_column(String)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow)

    dataset: Mapped["Dataset | None"] = relationship()
