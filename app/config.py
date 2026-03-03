"""
MST-Pipeline Configuration
"""
import os
from pathlib import Path

# --- Paths ---
PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
UPLOAD_DIR = DATA_DIR / "uploads"
DATASET_DIR = DATA_DIR / "datasets"
DB_DIR = PROJECT_DIR / "DB"
CONFIG_DIR = PROJECT_DIR / "config"
R_SCRIPTS_DIR = PROJECT_DIR / "r_scripts"

# --- Database ---
DATABASE_URL = f"sqlite:///{PROJECT_DIR / 'mst.db'}"

# --- Reference files ---
ECOLI_REF = CONFIG_DIR / "ecoli.fas"
SILVA_TRAIN_SET = CONFIG_DIR / "silva_nr99_v138.1_train_set_v4.fa.gz"
SILVA_SPECIES = CONFIG_DIR / "silva_nr99_v138.1_wSpecies_train_set_v4.fa.gz"
SILVA_SPECIES_ASSIGN = CONFIG_DIR / "silva_species_assignment_v138.1_v4.fa.gz"

# --- Source tracking DB files ---
SOURCE_TABLE = DB_DIR / "db_table.csv.gz"
SOURCE_TABLE_NORM = DB_DIR / "db_table.norm.csv.gz"
SOURCE_FASTA = DB_DIR / "db.fasta"
SOURCE_DESIGN = DB_DIR / "MST.design"

# --- Conda ---
CONDA_BASE = Path(os.environ.get("CONDA_BASE", Path.home() / "miniforge3"))
ST_CONDA_ENV = "ST"
DADA2_ENV_NAME = "mst"


def conda_cmd(args: list[str], env_name: str | None = None) -> list[str]:
    """Wrap a command to run inside a conda environment via conda run."""
    env = env_name or ST_CONDA_ENV
    return ["conda", "run", "-n", env, "--no-capture-output"] + args


# --- Pipeline defaults ---
DADA2_DEFAULTS = {
    "trim_left_f": 0,
    "trim_left_r": 0,
    "trunc_len_f": 250,
    "trunc_len_r": 200,
    "min_overlap": 20,
}

# (trunc_r1, trunc_r2) keyed by (region_tag, read_length)
TRUNC_DEFAULTS = {
    ("v34", 250): (0, 220),
    ("v34", 300): (120, 180),
    ("v4", 250): (180, 120),
    ("v4", 300): (180, 120),
    ("v45", 250): (200, 100),
    ("v45", 300): (200, 100),
}

DADA2_LONG_READ_DEFAULTS = {
    "trim_left_f": 0, "trim_left_r": 0,
    "trunc_len_f": 0, "trunc_len_r": 0,  # no truncation
    "min_overlap": 20,
    "max_ee": 10, "min_len": 1000, "max_len": 1600,
    "band_size": 32, "homopolymer_gap_penalty": -1,
}


def is_long_read(variable_region: str | None) -> bool:
    """Return True when variable_region indicates full-length 16S (V1-V9)."""
    return variable_region == "V1-V9"


VREGION_TAG = {"V4": "v4", "V3-V4": "v34", "V4-V5": "v45", "V1-V9": "v19"}
VALID_V4_REGIONS = {"V4", "V3-V4", "V4-V5", "V1-V9"}


def to_relative(abs_path: str | Path) -> str:
    """Convert an absolute path to a path relative to PROJECT_DIR for DB storage."""
    try:
        return str(Path(abs_path).relative_to(PROJECT_DIR))
    except ValueError:
        return str(abs_path)


def to_absolute(rel_path: str | None) -> Path | None:
    """Resolve a DB-stored path (relative or legacy absolute) to an absolute Path."""
    if not rel_path:
        return None
    p = Path(rel_path)
    if p.is_absolute():
        return p  # legacy absolute path — use as-is
    return PROJECT_DIR / p


# --- Server ---
HOST = "0.0.0.0"
PORT = 8050
MAX_THREADS = min(32, max(1, (os.cpu_count() or 1) - 1))
