"""
MST-Pipeline — Cutadapt primer/adapter trimming wrapper.

Uses region-specific primers with anchored matching (same as 16S-Pipeline).
"""
import logging
import os
import subprocess
from pathlib import Path

from app.config import is_long_read
from app.pipeline.detect import (
    REGION_PRIMERS,
    _primer_matches,
    _read_fastq_sequences,
    _reverse_complement,
    detect_sequencing_type,
    extract_sample_name,
)


def run_cutadapt(
    fastq_dir: Path,
    output_dir: Path,
    sequencing_type: str,
    variable_region: str | None,
    logger: logging.Logger,
    threads: int = 1,
) -> dict:
    """Run Cutadapt to remove primers from FASTQ files.

    Uses anchored (^) region-specific primers from REGION_PRIMERS.
    If variable_region is None or primers aren't detected, files are symlinked.

    Returns:
        dict with: success (bool), trimmed_dir (str), trim_stats (dict)
    """
    trimmed_dir = output_dir / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    if variable_region and variable_region in REGION_PRIMERS:
        fwd_primer = REGION_PRIMERS[variable_region]["forward"]
        rev_primer = REGION_PRIMERS[variable_region].get("reverse", "")
        logger.info(f"Trimming primers for {variable_region}: F={fwd_primer}, R={rev_primer}")
    else:
        logger.warning(
            f"Variable region '{variable_region}' not recognized or not set. "
            "Skipping primer trimming."
        )
        _symlink_files(fastq_dir, trimmed_dir)
        return {
            "success": True,
            "trimmed_dir": str(trimmed_dir),
            "trim_stats": {"skipped": True},
        }

    is_paired = sequencing_type == "paired-end"
    fastq_files = sorted(
        list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fq.gz"))
    )

    if not _primers_present(fastq_files, fwd_primer, logger):
        logger.info(
            "Primers not detected in reads — files appear already trimmed. "
            "Skipping Cutadapt."
        )
        _symlink_files(fastq_dir, trimmed_dir)
        return {
            "success": True,
            "trimmed_dir": str(trimmed_dir),
            "trim_stats": {"skipped": True, "reason": "primers_absent"},
        }

    trim_stats = {"total_reads": 0, "trimmed_reads": 0, "samples_processed": 0}

    if is_long_read(variable_region):
        _trim_long_read(
            fastq_files, trimmed_dir, fwd_primer, rev_primer, logger, trim_stats, threads
        )
    elif is_paired:
        _trim_paired(
            fastq_files, trimmed_dir, fwd_primer, rev_primer, logger, trim_stats, threads
        )
    else:
        _trim_single(
            fastq_files, trimmed_dir, fwd_primer, logger, trim_stats, threads
        )

    logger.info(
        f"Cutadapt done: {trim_stats['samples_processed']} samples, "
        f"{trim_stats['trimmed_reads']}/{trim_stats['total_reads']} reads trimmed"
    )

    return {
        "success": True,
        "trimmed_dir": str(trimmed_dir),
        "trim_stats": trim_stats,
    }


def _primers_present(
    fastq_files: list[Path],
    fwd_primer: str,
    logger: logging.Logger,
    n_reads: int = 100,
    threshold: float = 0.3,
) -> bool:
    """Check whether the forward primer is present in the first FASTQ file."""
    if not fastq_files:
        return False

    test_file = fastq_files[0]
    sequences = _read_fastq_sequences(test_file, n_reads)
    if not sequences:
        logger.warning(f"Could not read sequences from {test_file.name} for primer check")
        return True  # assume present if we can't check

    matches = sum(1 for seq in sequences if _primer_matches(seq, fwd_primer))
    ratio = matches / len(sequences)
    logger.info(
        f"Primer check: {matches}/{len(sequences)} reads ({ratio:.0%}) "
        f"start with forward primer in {test_file.name}"
    )
    return ratio >= threshold


def _trim_paired(
    fastq_files: list[Path],
    trimmed_dir: Path,
    fwd_primer: str,
    rev_primer: str,
    logger: logging.Logger,
    stats: dict,
    threads: int = 1,
):
    """Trim paired-end files sample by sample."""
    filenames = [f.name for f in fastq_files]
    detection = detect_sequencing_type(filenames)

    for sample_name, info in detection["samples"].items():
        r1_name = info.get("R1")
        r2_name = info.get("R2")
        if not r1_name or not r2_name:
            logger.warning(f"Skipping sample {sample_name}: missing R1 or R2")
            continue

        r1_in = fastq_files[0].parent / r1_name
        r2_in = fastq_files[0].parent / r2_name
        r1_out = trimmed_dir / r1_name
        r2_out = trimmed_dir / r2_name

        # Anchored primers — ^ prevents internal matches
        primer_args = ["-g", f"^{fwd_primer}"]
        if rev_primer:
            primer_args += ["-G", f"^{rev_primer}"]

        cmd = [
            "cutadapt",
            "-j", str(threads),
            *primer_args,
            "--minimum-length", "50",
            "-o", str(r1_out),
            "-p", str(r2_out),
            str(r1_in),
            str(r2_in),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed for {sample_name}: {result.stderr[:300]}")

        stats["samples_processed"] += 1
        _parse_cutadapt_stats(result.stderr, stats)
        logger.info(f"  Trimmed: {sample_name}")


def _trim_single(
    fastq_files: list[Path],
    trimmed_dir: Path,
    fwd_primer: str,
    logger: logging.Logger,
    stats: dict,
    threads: int = 1,
):
    """Trim single-end files."""
    for fq in fastq_files:
        out = trimmed_dir / fq.name
        cmd = [
            "cutadapt",
            "-j", str(threads),
            "-g", f"^{fwd_primer}",
            "--minimum-length", "50",
            "-o", str(out),
            str(fq),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt failed for {fq.name}: {result.stderr[:300]}")

        sample = extract_sample_name(fq.name)
        stats["samples_processed"] += 1
        _parse_cutadapt_stats(result.stderr, stats)
        logger.info(f"  Trimmed: {sample}")


def _trim_long_read(
    fastq_files: list[Path],
    trimmed_dir: Path,
    fwd_primer: str,
    rev_primer: str,
    logger: logging.Logger,
    stats: dict,
    threads: int = 1,
):
    """Trim long-read (full-length 16S) files using cutadapt with --rc.

    Handles both orientations of PacBio/Nanopore reads.
    Trims 27F from 5' and rc(1492R) from 3', with length filtering.
    """
    rev_primer_rc = _reverse_complement(rev_primer) if rev_primer else ""

    for fq in fastq_files:
        out = trimmed_dir / fq.name
        cmd = [
            "cutadapt",
            "-j", str(threads),
            "--rc",  # handle both read orientations
            "-g", f"^{fwd_primer}",  # 27F at 5'
        ]
        if rev_primer_rc:
            cmd += ["-a", f"{rev_primer_rc}$"]  # rc(1492R) at 3'
        cmd += [
            "--minimum-length", "1000",
            "--maximum-length", "1600",
            "--discard-untrimmed",
            "-o", str(out),
            str(fq),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            raise RuntimeError(f"Cutadapt (long-read) failed for {fq.name}: {result.stderr[:300]}")

        sample = extract_sample_name(fq.name)
        stats["samples_processed"] += 1
        _parse_cutadapt_stats(result.stderr, stats)
        logger.info(f"  Trimmed (long-read): {sample}")


def _parse_cutadapt_stats(stderr: str, stats: dict):
    """Parse cutadapt stderr for read counts."""
    for line in stderr.splitlines():
        if "Total reads processed:" in line:
            try:
                n = int(line.split(":")[-1].strip().replace(",", ""))
                stats["total_reads"] += n
            except ValueError:
                pass
        elif "Reads with adapters:" in line:
            try:
                n = int(line.split(":")[-1].strip().split("(")[0].strip().replace(",", ""))
                stats["trimmed_reads"] += n
            except ValueError:
                pass


def _symlink_files(src_dir: Path, dst_dir: Path):
    """Symlink all FASTQ files from src to dst."""
    for ext in ("*.fastq.gz", "*.fq.gz"):
        for f in src_dir.glob(ext):
            link = dst_dir / f.name
            if not link.exists():
                os.symlink(f, link)
