"""
MST-Pipeline — Database engine, session management, and initialization.
"""
from contextlib import contextmanager

from sqlalchemy import create_engine, text
from sqlalchemy.orm import Session, sessionmaker

from app.config import DATABASE_URL, DATASET_DIR, UPLOAD_DIR
from app.db.models import Base

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False},
)

SessionLocal = sessionmaker(
    autocommit=False, autoflush=False, bind=engine, expire_on_commit=False
)


def get_db():
    """FastAPI dependency: yields a DB session, auto-closes."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


@contextmanager
def get_session() -> Session:
    """Context manager for Dash callbacks (which bypass FastAPI DI)."""
    db = SessionLocal()
    try:
        yield db
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


def init_db():
    """Create all tables and ensure data directories exist."""
    Base.metadata.create_all(bind=engine)
    for d in [UPLOAD_DIR, DATASET_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # Safe column migration: add asv_table_path to samples if missing
    with engine.connect() as conn:
        cols = [r[1] for r in conn.execute(text("PRAGMA table_info(samples)")).fetchall()]
        if "asv_table_path" not in cols:
            conn.execute(text("ALTER TABLE samples ADD COLUMN asv_table_path VARCHAR"))
            conn.commit()
        if "asv_count" not in cols:
            conn.execute(text("ALTER TABLE samples ADD COLUMN asv_count INTEGER"))
            conn.commit()
