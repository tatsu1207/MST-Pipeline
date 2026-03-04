"""
MST-Pipeline — Auto-detection of sequencing type (SE/PE) and variable region.
"""
import gzip
import logging
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

# -- Filename patterns for paired-end detection --------------------------------

_PE_PATTERNS = [
    (r"^(.+?)_R1(_\d+)?\.fastq\.gz$", r"^(.+?)_R2(_\d+)?\.fastq\.gz$"),
    (r"^(.+?)_R1(_\d+)?\.fq\.gz$", r"^(.+?)_R2(_\d+)?\.fq\.gz$"),
    (r"^(.+?)_R1\.fastq\.gz$", r"^(.+?)_R2\.fastq\.gz$"),
    (r"^(.+?)_R1\.fq\.gz$", r"^(.+?)_R2\.fq\.gz$"),
    (r"^(.+?)_1\.fastq\.gz$", r"^(.+?)_2\.fastq\.gz$"),
    (r"^(.+?)_1\.fq\.gz$", r"^(.+?)_2\.fq\.gz$"),
]


def extract_sample_name(filename: str) -> str:
    """Extract sample name from a FASTQ filename by stripping read direction and extension."""
    name = filename
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    stripped = re.sub(r"_R[12](_\d+)?$", "", name)
    if stripped != name:
        return stripped
    return re.sub(r"_[12]$", "", name)


def detect_sequencing_type(filenames: list[str]) -> dict:
    """Detect single-end vs paired-end from a list of FASTQ filenames.

    Returns dict with keys: type, samples, errors
    """
    if not filenames:
        return {"type": "unknown", "samples": {}, "errors": ["No files provided."]}

    r1_files = {}
    r2_files = {}
    unmatched = []

    for fn in filenames:
        matched = False
        for r1_pat, r2_pat in _PE_PATTERNS:
            m1 = re.match(r1_pat, fn)
            if m1:
                sample = m1.group(1)
                r1_files[sample] = fn
                matched = True
                break
            m2 = re.match(r2_pat, fn)
            if m2:
                sample = m2.group(1)
                r2_files[sample] = fn
                matched = True
                break
        if not matched:
            unmatched.append(fn)

    all_samples = sorted(set(list(r1_files.keys()) + list(r2_files.keys())))
    samples = {}
    errors = []

    if all_samples:
        for sample in all_samples:
            r1 = r1_files.get(sample)
            r2 = r2_files.get(sample)
            if r1 and r2:
                samples[sample] = {"R1": r1, "R2": r2}
            elif r1 and not r2:
                errors.append(f"Sample '{sample}' has R1 but no R2: {r1}")
                samples[sample] = {"R1": r1, "R2": None}
            else:
                errors.append(f"Sample '{sample}' has R2 but no R1: {r2}")
                samples[sample] = {"R1": None, "R2": r2}

        if unmatched:
            errors.append(
                f"{len(unmatched)} file(s) don't match paired-end patterns: "
                f"{unmatched[:3]}{'...' if len(unmatched) > 3 else ''}"
            )

        has_pairs = any(s["R1"] and s["R2"] for s in samples.values())
        has_singles = bool(unmatched) or any(
            (s["R1"] is None) != (s["R2"] is None) for s in samples.values()
        )
        if has_pairs and not has_singles:
            seq_type = "paired-end"
        elif has_singles and not has_pairs:
            seq_type = "single-end"
        else:
            seq_type = "mixed"
    else:
        seq_type = "single-end"
        for fn in filenames:
            sample = extract_sample_name(fn)
            samples[sample] = {"R1": fn, "R2": None}

    return {"type": seq_type, "samples": samples, "errors": errors}


# -- Primer constants (IUPAC) -------------------------------------------------

PRIMERS = {
    "515F": "GTGYCAGCMGCCGCGGTAA",
    "806R": "GACTACNVGGGTWTCTAATCC",
    "515F_RC": "TTACCGCGGCKGCTGRCAC",
    "806R_RC": "GGATTAGAWACCCBNGTAGTC",
    "27F": "AGAGTTTGATCMTGGCTCAG",
    "1492R": "TACGGYTACCTTGTTACGACTT",
}

_COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVNacgtryMkswbdhvn",
                            "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")


def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of an IUPAC sequence."""
    return seq.translate(_COMPLEMENT)[::-1]

REGION_PRIMERS = {
    "V4": {
        "forward": "GTGYCAGCMGCCGCGGTAA",         # 515F
        "reverse": "GGACTACNVGGGTWTCTAATCC",        # 806R + 2bp boundary
        "expected_len": (282, 302),
    },
    "V3-V4": {
        "forward": "CCTACGGGNGGCWGCAG",             # 341F
        "reverse": "GACTACHVGGGTATCTAATCC",          # 806R
        "expected_len": (420, 480),
    },
    "V4-V5": {
        "forward": "GTGYCAGCMGCCGCGGTAA",           # 515F
        "reverse": "CCGYCAATTYMTTTRAGTTT",           # 926R
        "expected_len": (370, 430),
    },
    "V1-V9": {
        "forward": "AGAGTTTGATCMTGGCTCAG",           # 27F
        "reverse": "TACGGYTACCTTGTTACGACTT",          # 1492R
        "expected_len": (1400, 1500),
    },
}

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
    "S": "[GC]", "W": "[AT]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

_IUPAC_SET = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def _iupac_to_regex(seq: str) -> str:
    return "".join(IUPAC.get(c.upper(), c) for c in seq)


def _prefix_matches_primer(prefix: str, primer_seq: str, min_match: int = 15) -> bool:
    """Check if a read prefix starts with a primer (IUPAC-aware)."""
    if not isinstance(prefix, str) or len(prefix) < min_match:
        return False
    pattern = _iupac_to_regex(primer_seq[:min_match])
    return bool(re.match(pattern, prefix[:min_match], re.IGNORECASE))


def _primer_matches(sequence: str, primer: str, max_mismatches: int = 3) -> bool:
    """Check if a sequence starts with a primer, allowing degenerate bases and mismatches."""
    if len(sequence) < len(primer):
        return False
    mismatches = 0
    for seq_base, primer_base in zip(sequence.upper(), primer.upper()):
        allowed = _IUPAC_SET.get(primer_base, primer_base)
        if seq_base not in allowed:
            mismatches += 1
            if mismatches > max_mismatches:
                return False
    return True


def _read_fastq_sequences(fastq_path: Path, n_reads: int = 100) -> list[str]:
    """Read the first n sequences from a (gzipped) FASTQ file."""
    sequences = []
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    try:
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:
                    sequences.append(line.strip())
                    if len(sequences) >= n_reads:
                        break
    except Exception:
        pass
    return sequences


# -- V-region detection --------------------------------------------------------

VREGIONS = [
    ("V1", 69, 99),
    ("V2", 137, 242),
    ("V3", 433, 497),
    ("V4", 576, 682),
    ("V5", 822, 879),
    ("V6", 986, 1043),
    ("V7", 1117, 1173),
    ("V8", 1243, 1294),
    ("V9", 1435, 1465),
]

_REGION_COORDS = {
    "V1-V2": (69, 337),
    "V3-V4": (341, 805),
    "V4": (515, 806),
    "V4-V5": (515, 926),
    "V5-V6": (784, 1061),
}


def fast_identify_vregion(r1_prefix: str, r2_prefix: str) -> str | None:
    """Identify V-region from read prefixes without running bbmap.

    Returns the region string, or None when primers are absent/ambiguous.
    Checks for 27F (full-length 16S) before 515F (short-read regions).
    """
    # Check 27F first — long-read full-length 16S (V1-V9)
    r1_27F = _prefix_matches_primer(r1_prefix, PRIMERS["27F"])
    if r1_27F:
        return "V1-V9"
    # Also check reverse complement of 27F (PacBio reads can be either orientation)
    r1_27F_rc = _prefix_matches_primer(r1_prefix, _reverse_complement(PRIMERS["27F"]))
    if r1_27F_rc:
        return "V1-V9"

    r1_515F = _prefix_matches_primer(r1_prefix, PRIMERS["515F"])
    r2_806R = _prefix_matches_primer(r2_prefix, PRIMERS["806R"])
    if r1_515F and r2_806R:
        return "V4"
    if r1_515F:
        return "V4-V5"
    if r2_806R:
        return "V3-V4"
    return None


def identify_vregion(r1_path, r2_path, ref_path, n_reads=1000, timeout=60, threads=1):
    """Identify variable region using BBMap alignment to E. coli 16S reference.

    Supports both PE (r1+r2) and SE (r1 only, r2_path=None).
    """
    import os
    import shutil

    if not r1_path or not ref_path:
        return "N/A"
    if not os.path.exists(ref_path):
        return "N/A (ref missing)"
    if not shutil.which("bbmap.sh") or not shutil.which("reformat.sh"):
        return "N/A (bbmap not found)"

    tmpdir = tempfile.mkdtemp()
    try:
        # Subsample reads
        if r2_path and os.path.exists(r2_path):
            ref_result = subprocess.run(
                [
                    "reformat.sh",
                    f"in={r1_path}", f"in2={r2_path}",
                    f"out={tmpdir}/sub_R1.fastq.gz", f"out2={tmpdir}/sub_R2.fastq.gz",
                    f"samplereadstarget={n_reads}", "-Xmx1g",
                ],
                capture_output=True, timeout=timeout,
            )
            if ref_result.returncode != 0:
                bbmap_inputs = [(r1_path, "r1"), (r2_path, "r2")]
                use_reads_limit = True
            else:
                bbmap_inputs = [
                    (f"{tmpdir}/sub_R1.fastq.gz", "r1"),
                    (f"{tmpdir}/sub_R2.fastq.gz", "r2"),
                ]
                use_reads_limit = False
        else:
            # Single-end
            ref_result = subprocess.run(
                [
                    "reformat.sh",
                    f"in={r1_path}",
                    f"out={tmpdir}/sub_R1.fastq.gz",
                    f"samplereadstarget={n_reads}", "-Xmx1g",
                ],
                capture_output=True, timeout=timeout,
            )
            if ref_result.returncode != 0:
                bbmap_inputs = [(r1_path, "r1")]
                use_reads_limit = True
            else:
                bbmap_inputs = [(f"{tmpdir}/sub_R1.fastq.gz", "r1")]
                use_reads_limit = False

        all_starts, all_ends = [], []
        for infile, label in bbmap_inputs:
            extra = [f"reads={n_reads}"] if use_reads_limit else []
            bm = subprocess.run(
                [
                    "bbmap.sh", f"ref={ref_path}",
                    f"in={infile}", f"out={tmpdir}/{label}.sam",
                    "nodisk", f"threads={threads}", "ambiguous=best",
                    "minid=0.50", "maxindel=100", "bwr=0.25",
                    "semiperfectmode=f", "strictmaxindel=f", "-Xmx1g",
                ] + extra,
                capture_output=True, timeout=timeout,
            )
            if bm.returncode != 0 or not os.path.exists(f"{tmpdir}/{label}.sam"):
                continue
            starts, ends = _parse_sam(f"{tmpdir}/{label}.sam")
            all_starts.extend(starts)
            all_ends.extend(ends)

        if not all_starts:
            return "N/A (no reads mapped to reference)"

        all_starts.sort()
        all_ends.sort()
        amp_start = all_starts[len(all_starts) // 2]
        amp_end = all_ends[len(all_ends) // 2]

        covered = []
        for name, v_start, v_end in VREGIONS:
            ov_start = max(amp_start, v_start)
            ov_end = min(amp_end, v_end)
            if ov_end > ov_start:
                cov = (ov_end - ov_start) / (v_end - v_start)
                if cov >= 0.5:
                    covered.append(name)

        if not covered:
            return "Unknown"
        if len(covered) == 1:
            return covered[0]
        return f"{covered[0]}-{covered[-1]}"
    except Exception as e:
        return f"N/A ({e})"
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)


def _parse_sam(sam_file):
    starts, ends = [], []
    with open(sam_file) as f:
        for line in f:
            if line.startswith("@"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 10:
                continue
            flag = int(cols[1])
            if flag & 4:
                continue
            pos = int(cols[3])
            cigar = cols[5]
            alen = sum(int(l) for l, op in re.findall(r"(\d+)([MDNX=])", cigar))
            starts.append(pos)
            ends.append(pos + alen - 1)
    return starts, ends


def detect_variable_region(
    fastq_path: Path, n_reads: int = 100, r2_path: Path | None = None,
) -> dict:
    """Detect the 16S variable region from primer sequences in FASTQ reads.

    For PE data, pass *r2_path* to check reverse primers in R2 reads, which
    is essential for disambiguating regions that share the same forward primer
    (e.g. V4 vs V4-V5 both use 515F).

    Returns dict with keys: region, confidence, method, details
    """
    from app.config import ECOLI_REF

    sequences = _read_fastq_sequences(fastq_path, n_reads)
    if not sequences:
        return {
            "region": None, "confidence": 0.0,
            "method": "none", "details": "Could not read sequences from file.",
        }

    # Read R2 sequences when available (PE reverse-primer disambiguation)
    r2_sequences: list[str] = []
    if r2_path and Path(r2_path).exists():
        r2_sequences = _read_fastq_sequences(Path(r2_path), n_reads)

    # Read-length heuristic: if average read length > 1000bp, strongly favor V1-V9
    avg_len = sum(len(s) for s in sequences) / len(sequences) if sequences else 0

    # Count primer matches for each region (forward direction on R1)
    region_primers = {
        "V4": PRIMERS["515F"],
        "V3-V4": "CCTACGGGNGGCWGCAG",
        "V4-V5": PRIMERS["515F"],
        "V1-V9": PRIMERS["27F"],
    }

    region_scores = defaultdict(int)
    for seq in sequences:
        for region_name, primer in region_primers.items():
            if _primer_matches(seq, primer):
                region_scores[region_name] += 1
        # Also check reverse complement of 27F (PacBio reads can be either orientation)
        rc_27f = _reverse_complement(PRIMERS["27F"])
        if _primer_matches(seq, rc_27f):
            region_scores["V1-V9"] += 1

    # Reverse-primer scoring: check R2 start for reverse primers (PE)
    if r2_sequences:
        _rev_primers = {
            "V4": REGION_PRIMERS["V4"]["reverse"],
            "V3-V4": REGION_PRIMERS["V3-V4"]["reverse"],
            "V4-V5": REGION_PRIMERS["V4-V5"]["reverse"],
            "V1-V9": REGION_PRIMERS["V1-V9"]["reverse"],
        }
        for seq in r2_sequences:
            for region_name, rev_primer in _rev_primers.items():
                if _primer_matches(seq, rev_primer):
                    region_scores[region_name] += 1

    # Apply read-length bonus for V1-V9 when reads are long
    if avg_len > 1000 and "V1-V9" in region_scores:
        region_scores["V1-V9"] += len(sequences)
    elif avg_len > 1000 and not region_scores:
        # Long reads but no primer matches — likely V1-V9 with trimmed primers
        return {
            "region": "V1-V9",
            "confidence": 0.6,
            "method": "read_length_heuristic",
            "details": f"Avg read length {avg_len:.0f}bp suggests full-length 16S (V1-V9).",
        }

    if not region_scores:
        # Fallback: BBMap alignment
        bbmap_result = identify_vregion(str(fastq_path), None, str(ECOLI_REF))
        if bbmap_result and bbmap_result not in ("N/A", "Unknown") and not bbmap_result.startswith("N/A"):
            return {
                "region": bbmap_result, "confidence": 0.5,
                "method": "bbmap_alignment",
                "details": f"BBMap detected region: {bbmap_result}",
            }
        return {
            "region": None, "confidence": 0.0,
            "method": "none",
            "details": "No primer matches found and bbmap alignment inconclusive.",
        }

    best_region = max(region_scores, key=region_scores.get)
    best_count = region_scores[best_region]

    # Break ties using expected amplicon length proximity to avg read length.
    # Handles SE disambiguation (e.g. V4 ~250bp vs V4-V5 ~400bp) and cases
    # where R2 is unavailable.
    tied = [r for r, c in region_scores.items() if c == best_count]
    if len(tied) > 1 and avg_len > 0:
        def _len_distance(region):
            info = REGION_PRIMERS.get(region, {})
            lo, hi = info.get("expected_len", (0, 0))
            mid = (lo + hi) / 2
            return abs(avg_len - mid) if mid > 0 else float("inf")
        best_region = min(tied, key=_len_distance)
        best_count = region_scores[best_region]

    confidence = best_count / len(sequences)

    return {
        "region": best_region,
        "confidence": round(confidence, 2),
        "method": "primer_match",
        "details": f"{best_count}/{len(sequences)} reads matched {best_region} primer(s).",
    }
