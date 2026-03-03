"""
MST-Pipeline — Quality profiling and auto-detection of DADA2 truncation parameters.

Strategy:
  1. Determine insert length from the detected variable region and primer lengths.
  2. Compute minimum truncation lengths to guarantee >=20bp overlap for merging.
  3. Profile per-base quality and truncate from the tail where mean quality drops
     below Q20, but never below the overlap-required minimum.
  4. For single-end reads, truncate at the quality drop-off only.
"""
import gzip
import logging
from pathlib import Path

from app.pipeline.detect import REGION_PRIMERS, detect_sequencing_type

# Absolute floor — never truncate shorter than this
_MIN_TRUNC = 100

# Quality threshold: truncate where per-base mean Phred drops below this
_QUAL_FLOOR = 20


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
    """Auto-detect optimal DADA2 truncation lengths.

    For paired-end:
      - Computes insert length from variable region and primers.
      - Sets minimum truncation so that trunc_f + trunc_r >= insert + min_overlap.
      - Truncates tails where mean quality < Q20, R2 first (degrades faster).
      - Never truncates below the overlap-required minimum.

    For single-end:
      - Truncates at quality drop-off only.

    Returns dict with: trunc_len_f, trunc_len_r, details (str)
    """
    log = logger or logging.getLogger(__name__)

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

    # Profile quality scores
    log.info("Profiling R1 quality scores...")
    r1_scores = _read_quality_scores(r1_files[0], n_reads, skip_bases=primer_offset_f)
    r1_median_len = _median_read_length(r1_scores)

    if not is_paired:
        trunc_f = _quality_cutoff(r1_scores)
        r1_pass = _pass_rate_at(r1_scores, trunc_f)
        details = f"trunc_len_f={trunc_f} ({r1_pass:.0%} pass EE≤5)"
        log.info(f"Auto-detected: {details}")
        return {"trunc_len_f": trunc_f, "trunc_len_r": 0, "details": details}

    log.info("Profiling R2 quality scores...")
    r2_scores = _read_quality_scores(r2_files[0], n_reads, skip_bases=primer_offset_r)
    r2_median_len = _median_read_length(r2_scores)

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

    # Compute quality-based truncation (where mean Q drops below threshold)
    qual_trunc_f = _quality_cutoff(r1_scores)
    qual_trunc_r = _quality_cutoff(r2_scores)

    if insert_len:
        # Minimum truncation lengths for desired overlap
        min_total = insert_len + min_overlap

        # Start from quality-based truncation
        trunc_f = qual_trunc_f
        trunc_r = qual_trunc_r

        # If we don't meet the overlap requirement, extend reads back
        if trunc_f + trunc_r < min_total:
            shortfall = min_total - (trunc_f + trunc_r)
            # Extend R1 first (usually better quality), then R2
            extend_f = min(shortfall, r1_median_len - trunc_f)
            extend_f = max(0, extend_f)
            trunc_f += extend_f
            shortfall -= extend_f
            if shortfall > 0:
                extend_r = min(shortfall, r2_median_len - trunc_r)
                extend_r = max(0, extend_r)
                trunc_r += extend_r

        overlap = trunc_f + trunc_r - insert_len
        log.info(f"Quality cutoff: R1={qual_trunc_f}, R2={qual_trunc_r}")
        log.info(f"Overlap-adjusted: R1={trunc_f}, R2={trunc_r}, overlap={overlap}bp")
    else:
        # No region info — just use quality cutoff
        trunc_f = qual_trunc_f
        trunc_r = qual_trunc_r
        log.info(f"No region info, using quality cutoff: R1={trunc_f}, R2={trunc_r}")

    # Clamp to valid range
    trunc_f = max(_MIN_TRUNC, min(trunc_f, r1_median_len))
    trunc_r = max(_MIN_TRUNC, min(trunc_r, r2_median_len))

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
        "details": details,
    }


def _read_quality_scores(
    fastq_path: Path, n_reads: int = 5000, skip_bases: int = 0
) -> list[list[int]]:
    """Read Phred+33 quality scores from first n_reads of a FASTQ(.gz) file."""
    scores = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 0:
                    qual = line.rstrip()
                    phred = [ord(c) - 33 for c in qual]
                    if skip_bases:
                        phred = phred[skip_bases:]
                    if phred:
                        scores.append(phred)
                    if len(scores) >= n_reads:
                        break
    except Exception:
        pass
    return scores


def _median_read_length(quality_scores: list[list[int]]) -> int:
    """Return the median read length from quality score lists."""
    if not quality_scores:
        return 0
    lengths = sorted(len(s) for s in quality_scores)
    return lengths[len(lengths) // 2]


def _quality_cutoff(quality_scores: list[list[int]], q_floor: int = _QUAL_FLOOR) -> int:
    """Find the position where per-base mean quality drops below q_floor.

    Scans from the end of the read backward. Returns the last position
    where a sliding window of 10bp has mean quality >= q_floor.
    """
    if not quality_scores:
        return _MIN_TRUNC

    # Compute per-position mean quality
    max_len = max(len(s) for s in quality_scores)
    mean_qual = []
    for pos in range(max_len):
        quals = [s[pos] for s in quality_scores if pos < len(s)]
        if quals:
            mean_qual.append(sum(quals) / len(quals))
        else:
            mean_qual.append(0)

    # Use a 10bp sliding window from the end to find the quality drop-off
    window = 10
    cutoff = max_len
    for pos in range(max_len - 1, window - 2, -1):
        start = max(0, pos - window + 1)
        win_mean = sum(mean_qual[start:pos + 1]) / (pos + 1 - start)
        if win_mean >= q_floor:
            cutoff = pos + 1  # truncate at this length (1-indexed)
            break
    else:
        # Quality never reaches threshold — use minimum
        cutoff = _MIN_TRUNC

    return max(_MIN_TRUNC, cutoff)


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
