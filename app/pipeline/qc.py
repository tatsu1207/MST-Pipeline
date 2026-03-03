"""
MST-Pipeline — QC stats and primer/case detection.

Extracted from MST/scripts/app_qc.py: pure functions only, no UI code.
"""
import gzip
import os
from collections import Counter

from app.pipeline.detect import (
    PRIMERS,
    _prefix_matches_primer,
    fast_identify_vregion,
    identify_vregion,
)
from app.config import ECOLI_REF, VALID_V4_REGIONS, VREGION_TAG

# -- 12-case matrix description ------------------------------------------------
CASE_TABLE = """
| Case | Region | Length | Primers | cutadapt flags |
|------|--------|--------|---------|----------------|
| 1 | V4 | 250 bp | Included | `-g 515F -G 806R` |
| 2 | V4 | 250 bp | Trimmed | *(no action)* |
| 3 | V4 | 300 bp | Included | `-g 515F -G 806R -a 806R_RC -A 515F_RC` |
| 4 | V4 | 300 bp | Trimmed | `-a 806R_RC -A 515F_RC` |
| 5 | V3-V4 | 250 bp | Included | `-g 515F --overlap 16 -G 806R` |
| 6 | V3-V4 | 250 bp | Trimmed | `-g 515F --overlap 16` |
| 7 | V3-V4 | 300 bp | Included | `-g 515F --overlap 16 -G 806R -A 515F_RC` |
| 8 | V3-V4 | 300 bp | Trimmed | `-g 515F --overlap 16 -A 515F_RC` |
| 9 | V4-V5 | 250 bp | Included | `-g 515F -G 806R --overlap 16 --discard-untrimmed` |
| 10 | V4-V5 | 250 bp | Trimmed | *(min-length only)* |
| 11 | V4-V5 | 300 bp | Included | `-g 515F -G 806R -a 806R_RC --overlap 16 --discard-untrimmed` |
| 12 | V4-V5 | 300 bp | Trimmed | `-a 806R_RC -A 806R --overlap 16 --discard-untrimmed` |
"""


def get_comprehensive_stats(file_path, n_bases=20, limit=50000):
    """Return (prefix, avg_len, read_count, q30_drop, q25_drop, q20_drop, total_reads)."""
    if not file_path or not os.path.exists(file_path):
        return "N/A", 0, 0, 0, 0, 0, 0

    bp_counts = Counter()
    pos_stats = []
    total_len = 0
    read_count = 0
    total_reads = 0

    try:
        with gzip.open(file_path, "rt") as f:
            for i, line in enumerate(f):
                line_type = i % 4
                if line_type == 1:
                    total_reads += 1
                    if read_count < limit:
                        seq = line.strip()
                        s_len = len(seq)
                        total_len += s_len
                        read_count += 1
                        if s_len >= n_bases:
                            bp_counts[seq[:n_bases]] += 1
                elif line_type == 3 and read_count <= limit:
                    q_str = line.strip()
                    for pos, char in enumerate(q_str):
                        q_score = ord(char) - 33
                        if len(pos_stats) <= pos:
                            pos_stats.append([0, 0])
                        pos_stats[pos][0] += q_score
                        pos_stats[pos][1] += 1
    except Exception:
        return "ERR", 0, 0, 0, 0, 0, 0

    if read_count == 0:
        return "N/A", 0, 0, 0, 0, 0, 0

    avg_len = total_len / read_count
    most_common = bp_counts.most_common(1)[0][0] if bp_counts else "N/A"
    averages = [s[0] / s[1] if s[1] > 0 else 0 for s in pos_stats]

    def find_sustained_drop(avg_list, threshold, window=5):
        for idx in range(len(avg_list) - window):
            if sum(avg_list[idx: idx + window]) / window < threshold:
                return idx
        return len(avg_list)

    q30_drop = int(min(avg_len, find_sustained_drop(averages, 30)))
    q25_drop = int(min(avg_len, find_sustained_drop(averages, 25)))
    q20_drop = int(min(avg_len, find_sustained_drop(averages, 20)))

    return most_common, avg_len, read_count, q30_drop, q25_drop, q20_drop, total_reads


def detect_read_length(avg_len: float) -> int:
    """Classify average read length into 250 or 300 bp."""
    return 250 if avg_len <= 275 else 300


def detect_primer_status(r1_prefix: str, vregion: str, r2_prefix: str = None) -> str:
    """Detect whether amplification primers are still in the reads."""
    if vregion in ("V4", "V4-V5"):
        if _prefix_matches_primer(r1_prefix, PRIMERS["515F"]):
            return "Included"
        return "Trimmed"
    if r2_prefix and _prefix_matches_primer(r2_prefix, PRIMERS["806R"]):
        return "Included"
    return "Trimmed"


def detect_r2_pad(r2_path, limit=1000) -> int:
    """Detect dark-cycle pad at the 5' end of R2 reads."""
    if not r2_path or not os.path.exists(r2_path):
        return 0
    n_at_pos = {}
    count = 0
    try:
        with gzip.open(r2_path, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    seq = line.strip()
                    for pos in range(min(5, len(seq))):
                        if seq[pos] == "N":
                            n_at_pos[pos] = n_at_pos.get(pos, 0) + 1
                    count += 1
                    if count >= limit:
                        break
    except Exception:
        return 0
    if count == 0:
        return 0
    pad = 0
    for pos in range(5):
        if n_at_pos.get(pos, 0) / count > 0.5:
            pad = pos + 1
    return pad


def get_cutadapt_args(region: str, length: int, primer_status: str):
    """Return (adapter_args_list, case_number, error_msg).

    For paired-end reads. Includes --minimum-length per case.
    """
    included = primer_status == "Included"
    minlen_v4 = ["--minimum-length", "150:150"]
    minlen_v34 = ["--minimum-length", "40:150"]

    if region == "V4":
        if length == 250:
            if included:
                return ["-g", PRIMERS["515F"], "-G", PRIMERS["806R"]] + minlen_v4, 1, None
            return [], 2, None
        if included:
            return [
                "-g", PRIMERS["515F"], "-G", PRIMERS["806R"],
                "-a", PRIMERS["806R_RC"], "-A", PRIMERS["515F_RC"],
            ] + minlen_v4, 3, None
        return ["-a", PRIMERS["806R_RC"], "-A", PRIMERS["515F_RC"]] + minlen_v4, 4, None

    if region == "V3-V4":
        if length == 250:
            if included:
                return ["-g", PRIMERS["515F"], "--overlap", "16",
                        "-G", PRIMERS["806R"]] + minlen_v34, 5, None
            return ["-g", PRIMERS["515F"], "--overlap", "16"] + minlen_v34, 6, None
        if included:
            return ["-g", PRIMERS["515F"], "--overlap", "16",
                    "-G", PRIMERS["806R"], "-A", PRIMERS["515F_RC"]] + minlen_v34, 7, None
        return ["-g", PRIMERS["515F"], "--overlap", "16",
                "-A", PRIMERS["515F_RC"]] + minlen_v34, 8, None

    if region == "V4-V5":
        if length == 250:
            minlen_v45_250 = ["--minimum-length", "150:150"]
            if included:
                return [
                    "-g", PRIMERS["515F"], "-G", PRIMERS["806R"],
                    "--overlap", "16", "--discard-untrimmed",
                ] + minlen_v45_250, 9, None
            return minlen_v45_250, 10, None
        if included:
            return [
                "-g", PRIMERS["515F"], "-G", PRIMERS["806R"],
                "-a", PRIMERS["806R_RC"], "--times", "2",
                "--overlap", "16", "--discard-untrimmed",
                "--minimum-length", "200:40",
            ], 11, None
        return [
            "-a", PRIMERS["806R_RC"], "-A", PRIMERS["806R"],
            "--overlap", "16", "--discard-untrimmed",
            "--minimum-length", "200:40",
        ], 12, None

    return None, 0, f"Unknown region: {region}"


def get_cutadapt_args_se(region: str, length: int, primer_status: str):
    """Return (adapter_args_list, case_number, error_msg) for single-end reads.

    6 cases: 3 regions x 2 primer states.
    """
    included = primer_status == "Included"
    minlen = ["--minimum-length", "100"]

    if region == "V4":
        if included:
            return ["-g", PRIMERS["515F"]] + minlen, 1, None
        return [] + minlen, 2, None
    if region == "V3-V4":
        if included:
            return ["-g", PRIMERS["515F"], "--overlap", "16"] + minlen, 3, None
        return ["-g", PRIMERS["515F"], "--overlap", "16"] + minlen, 4, None
    if region == "V4-V5":
        if included:
            return ["-g", PRIMERS["515F"], "--overlap", "16", "--discard-untrimmed"] + minlen, 5, None
        return minlen, 6, None

    return None, 0, f"Unknown region: {region}"


def process_qc_sample(s_name, paths, n_bases=20, read_limit=50000, threads=1):
    """Run QC + auto-detection for one sample. Returns a row dict.

    Supports both PE (paths has R1+R2) and SE (paths has R1 only).
    """
    r1_stats = get_comprehensive_stats(paths.get("R1"), n_bases, read_limit)
    r2_path = paths.get("R2")
    is_se = not r2_path or not os.path.exists(r2_path)

    if is_se:
        r2_stats = ("N/A", 0, 0, 0, 0, 0, 0)
        vregion = None
    else:
        r2_stats = get_comprehensive_stats(r2_path, n_bases, read_limit)
        vregion = fast_identify_vregion(r1_stats[0], r2_stats[0])

    if vregion is None:
        if is_se:
            vregion = identify_vregion(paths.get("R1"), None, str(ECOLI_REF), threads=threads)
        else:
            vregion = identify_vregion(paths.get("R1"), r2_path, str(ECOLI_REF), threads=threads)

    det_length = detect_read_length(r1_stats[1])

    if is_se:
        det_primers = detect_primer_status(r1_stats[0], vregion)
        args, case_num, error = get_cutadapt_args_se(vregion, det_length, det_primers)
    else:
        det_primers = detect_primer_status(r1_stats[0], vregion, r2_stats[0])
        args, case_num, error = get_cutadapt_args(vregion, det_length, det_primers)
        r2_pad = detect_r2_pad(r2_path)
        if r2_pad > 0 and "-G" not in (args or []):
            if args is None:
                args = []
            args = ["-U", str(r2_pad)] + args

    return {
        "sample": s_name,
        "sequencing_type": "single-end" if is_se else "paired-end",
        "r1_reads": r1_stats[6], "r1_prefix": r1_stats[0],
        "r1_avg_len": r1_stats[1], "r1_q30": r1_stats[3],
        "r1_q25": r1_stats[4], "r1_q20": r1_stats[5],
        "r2_reads": r2_stats[6], "r2_prefix": r2_stats[0],
        "r2_avg_len": r2_stats[1], "r2_q30": r2_stats[3],
        "r2_q25": r2_stats[4], "r2_q20": r2_stats[5],
        "vregion": vregion, "det_length": det_length,
        "det_primers": det_primers, "det_case": case_num,
        "cutadapt_args": args, "cutadapt_error": error,
    }
