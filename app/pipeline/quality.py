"""
MST-Pipeline — Auto-detection of DADA2 truncation parameters.

Strategy:
  1. Truncate R1 and R2 at 3/4 of read length (removes low-quality tails).
  2. For paired-end, ensure sufficient overlap for merging (extend if needed).
  3. V4 region extraction happens at the source tracking step, not here.
"""
import gzip
import logging
from pathlib import Path

from app.config import is_long_read
from app.pipeline.detect import REGION_PRIMERS, detect_sequencing_type

# Absolute floor — never truncate shorter than this
_MIN_TRUNC = 100


def detect_truncation_params(
    trimmed_dir: Path,
    sequencing_type: str,
    variable_region: str | None,
    min_overlap: int = 20,
    n_reads: int = 5000,
    primer_offset_f: int = 0,
    primer_offset_r: int = 0,
    logger: logging.Logger | None = None,
    **_kwargs,
) -> dict:
    """Auto-detect DADA2 truncation lengths.

    Truncates at 3/4 of read length for both R1 and R2, then ensures
    sufficient overlap for paired-end merging.

    Returns dict with: trunc_len_f, trunc_len_r, details (str)
    """
    log = logger or logging.getLogger(__name__)

    # Long-read mode: no truncation needed
    if is_long_read(variable_region):
        log.info("Long-read mode: skipping truncation (reads processed at full length)")
        return {"trunc_len_f": 0, "trunc_len_r": 0, "details": "Long-read mode: no truncation"}

    # Find FASTQ files
    all_fastq = sorted(trimmed_dir.glob("*.fastq.gz")) + sorted(
        trimmed_dir.glob("*.fq.gz")
    )
    if not all_fastq:
        log.warning("No FASTQ files found for quality profiling")
        return {"trunc_len_f": 0, "trunc_len_r": 0, "details": "No files found"}

    filenames = [f.name for f in all_fastq]
    file_map = {f.name: f for f in all_fastq}
    detection = detect_sequencing_type(filenames)

    r1_files = []
    r2_files = []
    for sample_info in detection["samples"].values():
        if sample_info.get("R1") and sample_info["R1"] in file_map:
            r1_files.append(file_map[sample_info["R1"]])
        if sample_info.get("R2") and sample_info["R2"] in file_map:
            r2_files.append(file_map[sample_info["R2"]])

    if not r1_files:
        r1_files = all_fastq

    is_paired = sequencing_type == "paired-end" and r2_files

    # Detect leading N bases that must be trimmed (e.g. dark-cycle artifacts)
    # Check ALL samples since DADA2 trimLeft applies globally
    trim_left_f = max(_detect_trim_left(f) for f in r1_files)
    if trim_left_f:
        log.info(f"R1: detected {trim_left_f} leading N/low-quality bases → trimLeft={trim_left_f}")

    # Profile quality scores from ALL samples (skip the bases that will be left-trimmed)
    per_sample = max(1, n_reads // len(r1_files))
    log.info(f"Profiling R1 quality scores ({len(r1_files)} files, {per_sample} reads each)...")
    r1_scores = []
    r1_per_sample_lens = []
    for f in r1_files:
        sample_scores = _read_quality_scores(f, per_sample, skip_bases=primer_offset_f + trim_left_f)
        r1_scores.extend(sample_scores)
        r1_per_sample_lens.append(_median_read_length(sample_scores))

    # Warn if samples have very different read lengths (>30bp spread)
    if r1_per_sample_lens and max(r1_per_sample_lens) - min(r1_per_sample_lens) > 30:
        log.warning(
            f"Mixed read lengths detected in R1: {min(r1_per_sample_lens)}-{max(r1_per_sample_lens)}bp. "
            "Consider running samples with different read lengths in separate DADA2 runs."
        )

    # Use the shortest sample's median to avoid discarding shorter reads
    r1_median_len = min(r1_per_sample_lens) if r1_per_sample_lens else 0

    if not is_paired:
        trunc_f = max(_MIN_TRUNC, r1_median_len * 3 // 4)
        r1_pass = _pass_rate_at(r1_scores, trunc_f)
        details = f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤5)"
        if trim_left_f:
            details += f", trim_left_f={trim_left_f}"
        log.info(f"Auto-detected: {details}")
        return {
            "trunc_len_f": trunc_f, "trunc_len_r": 0,
            "trim_left_f": trim_left_f, "trim_left_r": 0,
            "details": details,
        }

    trim_left_r = max(_detect_trim_left(f) for f in r2_files)
    if trim_left_r:
        log.info(f"R2: detected {trim_left_r} leading N/low-quality bases → trimLeft={trim_left_r}")

    per_sample_r2 = max(1, n_reads // len(r2_files))
    log.info(f"Profiling R2 quality scores ({len(r2_files)} files, {per_sample_r2} reads each)...")
    r2_scores = []
    r2_per_sample_lens = []
    for f in r2_files:
        sample_scores = _read_quality_scores(f, per_sample_r2, skip_bases=primer_offset_r + trim_left_r)
        r2_scores.extend(sample_scores)
        r2_per_sample_lens.append(_median_read_length(sample_scores))
    if r2_per_sample_lens and max(r2_per_sample_lens) - min(r2_per_sample_lens) > 30:
        log.warning(
            f"Mixed read lengths detected in R2: {min(r2_per_sample_lens)}-{max(r2_per_sample_lens)}bp. "
            "Consider running samples with different read lengths in separate DADA2 runs."
        )
    r2_median_len = min(r2_per_sample_lens) if r2_per_sample_lens else 0

    # Determine insert length
    insert_len = None
    primer_info = REGION_PRIMERS.get(variable_region, {})
    amplicon_len = _expected_amplicon_length(variable_region)
    if amplicon_len:
        fwd_len = len(primer_info.get("forward", ""))
        rev_len = len(primer_info.get("reverse", ""))
        insert_len = amplicon_len - fwd_len - rev_len
        log.info(
            f"Amplicon ~{amplicon_len}bp, primers {fwd_len}+{rev_len}bp, "
            f"insert ~{insert_len}bp"
        )

    # Truncate at 3/4 of read length
    trunc_f = max(_MIN_TRUNC, r1_median_len * 3 // 4)
    trunc_r = max(_MIN_TRUNC, r2_median_len * 3 // 4)

    # If pass rate is too low at 3/4, shorten until ≥50% pass
    trunc_f = _shorten_for_quality(r1_scores, trunc_f, log, "R1")
    trunc_r = _shorten_for_quality(r2_scores, trunc_r, log, "R2")

    if insert_len:
        min_total = insert_len + min_overlap
        if trunc_f + trunc_r >= min_total:
            overlap = trunc_f + trunc_r - insert_len
            log.info(f"Truncation: R1={trunc_f}, R2={trunc_r}, overlap={overlap}bp")
        else:
            # Not enough overlap — fall back to full read length
            trunc_f = r1_median_len
            trunc_r = r2_median_len
            overlap = trunc_f + trunc_r - insert_len
            log.info(f"3/4 rule insufficient overlap, using full length: R1={trunc_f}, R2={trunc_r}, overlap={overlap}bp")
    else:
        log.info(f"Truncation: R1={trunc_f} (from {r1_median_len}bp), R2={trunc_r} (from {r2_median_len}bp)")

    r1_pass = _pass_rate_at(r1_scores, trunc_f)
    r2_pass = _pass_rate_at(r2_scores, trunc_r)
    combined = r1_pass * r2_pass

    details_parts = [
        f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤5)",
        f"trunc_len_r={trunc_r} ({r2_pass:.0%} pass EE≤5)",
        f"~{combined:.0%} combined retention",
    ]
    if insert_len:
        overlap = trunc_f + trunc_r - insert_len
        details_parts.append(f"overlap={overlap}bp")

    details = ", ".join(details_parts)
    log.info(f"Auto-detected: {details}")

    return {
        "trunc_len_f": trunc_f,
        "trunc_len_r": trunc_r,
        "trim_left_f": trim_left_f,
        "trim_left_r": trim_left_r,
        "details": details,
    }


def _detect_trim_left(fastq_path: Path, n_reads: int = 1000, threshold: float = 0.5) -> int:
    """Detect how many leading bases should be trimmed due to N calls or Q≤2.

    Returns the number of contiguous leading positions where >threshold fraction
    of reads have N or Q≤2.  Typically 0; non-zero for datasets with systematic
    dark-cycle issues (e.g. position 2 of R2 on some Illumina runs).
    """
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    seqs: list[str] = []
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:  # sequence line
                    seqs.append(line.rstrip().upper())
                    if len(seqs) >= n_reads:
                        break
    except Exception:
        return 0

    if not seqs:
        return 0

    min_len = min(len(s) for s in seqs)
    trim = 0
    for pos in range(min(min_len, 10)):  # check first 10 positions
        n_bad = sum(1 for s in seqs if s[pos] == "N") / len(seqs)
        if n_bad >= threshold:
            trim = pos + 1  # must trim through this position
    return trim


def _read_quality_scores(
    fastq_path: Path, n_reads: int = 5000, skip_bases: int = 0
) -> list[list[int]]:
    """Read Phred+33 quality scores sampled evenly across a FASTQ(.gz) file.

    Samples reads from throughout the file (not just the beginning) to get
    a representative quality profile.
    """
    # First pass: count total reads
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    total = 0
    try:
        with opener(fastq_path, "rt") as f:
            for line in f:
                total += 1
        total = total // 4
    except Exception:
        total = 0

    if total == 0:
        return []

    # Compute stride to sample evenly; stride=1 means read every record
    stride = max(1, total // n_reads)

    scores = []
    try:
        with opener(fastq_path, "rt") as f:
            read_idx = 0
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 0:
                    if read_idx % stride == 0:
                        qual = line.rstrip()
                        phred = [ord(c) - 33 for c in qual]
                        if skip_bases:
                            phred = phred[skip_bases:]
                        if phred:
                            scores.append(phred)
                        if len(scores) >= n_reads:
                            break
                    read_idx += 1
    except Exception:
        pass
    return scores


def _median_read_length(quality_scores: list[list[int]]) -> int:
    """Return the median read length from quality score lists."""
    if not quality_scores:
        return 0
    lengths = sorted(len(s) for s in quality_scores)
    return lengths[len(lengths) // 2]


def _shorten_for_quality(
    quality_scores: list[list[int]], trunc_len: int,
    log: logging.Logger, label: str,
    min_rate: float = 0.5, max_ee: float = 5.0,
) -> int:
    """Shorten truncation length until ≥min_rate fraction of reads pass EE≤max_ee."""
    rate = _pass_rate_at(quality_scores, trunc_len, max_ee)
    if rate >= min_rate:
        return trunc_len
    original = trunc_len
    while trunc_len > _MIN_TRUNC:
        trunc_len -= 10
        trunc_len = max(trunc_len, _MIN_TRUNC)
        rate = _pass_rate_at(quality_scores, trunc_len, max_ee)
        if rate >= min_rate:
            break
    log.info(f"{label}: 3/4 rule ({original}bp) had {_pass_rate_at(quality_scores, original, max_ee):.0%} pass rate, shortened to {trunc_len}bp ({rate:.0%} pass)")
    return trunc_len


def _pass_rate_at(quality_scores: list[list[int]], trunc_len: int, max_ee: float = 5.0) -> float:
    """Calculate fraction of ALL reads passing EE<=max_ee at a given truncation length.

    Reads shorter than trunc_len count as failures (DADA2 discards them).
    """
    if not quality_scores or trunc_len <= 0:
        return 0.0
    passing = 0
    for scores in quality_scores:
        if len(scores) >= trunc_len:
            ee = sum(10 ** (-q / 10) for q in scores[:trunc_len])
            if ee <= max_ee:
                passing += 1
    return passing / len(quality_scores)


def _expected_amplicon_length(variable_region: str | None) -> int | None:
    """Get the expected amplicon length (midpoint) for a variable region."""
    if not variable_region or variable_region not in REGION_PRIMERS:
        return None
    lo, hi = REGION_PRIMERS[variable_region]["expected_len"]
    return (lo + hi) // 2
