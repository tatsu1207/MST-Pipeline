"""
MST-Pipeline — DADA2 R script wrapper.

Calls the external run_dada2.R script (same as 16S-Pipeline).
Supports both paired-end and single-end reads.
"""
import json
import logging
import subprocess
from collections.abc import Callable
from pathlib import Path

import pandas as pd

from app.config import DADA2_ENV_NAME, R_SCRIPTS_DIR, conda_cmd


def run_dada2(
    input_dir: Path,
    output_dir: Path,
    sequencing_type: str,
    trim_left_f: int,
    trim_left_r: int,
    trunc_len_f: int,
    trunc_len_r: int,
    min_overlap: int,
    threads: int,
    logger: logging.Logger,
    proc_callback: Callable[[subprocess.Popen], None] | None = None,
    error_model: str = "default",
    band_size: int = 16,
    homopolymer_gap_penalty: int = 0,
    max_ee: float = 5.0,
    min_len: int = 0,
    max_len: int = 0,
) -> dict:
    """Run DADA2 denoising via the R script.

    Returns:
        dict with: asv_table_path, track_reads_path,
                   asv_count, sample_count, track_reads (per-sample dict)
    """
    dada2_dir = output_dir / "dada2"
    dada2_dir.mkdir(parents=True, exist_ok=True)

    is_longread = error_model == "PacBio"

    if is_longread:
        mode = "longread"
        cmd = conda_cmd([
            "Rscript", str(R_SCRIPTS_DIR / "run_dada2.R"),
            "--input_dir", str(input_dir),
            "--output_dir", str(dada2_dir),
            "--mode", mode,
            "--threads", str(threads),
            "--min_len", str(min_len),
            "--max_len", str(max_len),
            "--max_ee", str(int(max_ee)),
            "--band_size", str(band_size),
        ], env_name=DADA2_ENV_NAME)
    else:
        mode = "paired" if sequencing_type == "paired-end" else "single"
        cmd = conda_cmd([
            "Rscript", str(R_SCRIPTS_DIR / "run_dada2.R"),
            "--input_dir", str(input_dir),
            "--output_dir", str(dada2_dir),
            "--mode", mode,
            "--trim_left_f", str(trim_left_f),
            "--trim_left_r", str(trim_left_r),
            "--trunc_len_f", str(trunc_len_f),
            "--trunc_len_r", str(trunc_len_r),
            "--min_overlap", str(min_overlap),
            "--threads", str(threads),
        ], env_name=DADA2_ENV_NAME)

    logger.info(f"Running DADA2 ({mode} mode)...")

    # Stream output line by line
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        start_new_session=True,
    )
    if proc_callback:
        proc_callback(proc)

    last_line = ""
    for line in proc.stdout:
        stripped = line.rstrip()
        if stripped:
            logger.info(f"[DADA2] {stripped}")
            last_line = stripped

    proc.wait()

    if proc.returncode != 0:
        raise RuntimeError(f"DADA2 R script failed (exit code {proc.returncode})")

    # Parse JSON status from last line of output
    r_status = {}
    try:
        r_status = json.loads(last_line)
    except json.JSONDecodeError:
        logger.warning(f"Could not parse DADA2 JSON status: {last_line}")

    if r_status.get("status") == "error":
        raise RuntimeError(f"DADA2 error: {r_status.get('message', 'unknown')}")

    # Verify output files exist
    asv_table_path = dada2_dir / "asv_table.tsv"
    track_reads_path = dada2_dir / "track_reads.tsv"

    for p in [asv_table_path, track_reads_path]:
        if not p.exists():
            raise FileNotFoundError(f"Expected DADA2 output not found: {p}")

    # Parse track_reads for per-sample read counts
    track_df = pd.read_csv(track_reads_path, sep="\t")
    track_reads = {}
    for _, row in track_df.iterrows():
        track_reads[row["sample"]] = {
            "input": int(row.get("input", 0)),
            "filtered": int(row.get("filtered", 0)),
            "denoised": int(row.get("denoised", row.get("denoisedF", 0))),
            "nonchim": int(row.get("nonchim", 0)),
        }

    asv_count = r_status.get("asv_count", 0)
    sample_count = r_status.get("sample_count", len(track_reads))

    logger.info(f"DADA2 complete: {sample_count} samples, {asv_count} ASVs")

    return {
        "asv_table_path": str(asv_table_path),
        "track_reads_path": str(track_reads_path),
        "asv_count": asv_count,
        "sample_count": sample_count,
        "track_reads": track_reads,
    }
