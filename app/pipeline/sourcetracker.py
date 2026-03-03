"""
MST-Pipeline — Source tracking logic.

Extracted from MST/scripts/app_sourcetracker_ui.py: pure functions only.
Handles feature alignment, group collapsing, and Gibbs sampling dispatch.
"""
import io
import os
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from app.config import SOURCE_TABLE, SOURCE_FASTA, SOURCE_DESIGN, ST_CONDA_ENV

_HERE = Path(__file__).resolve().parent.parent.parent
_GIBBS_SCRIPT = _HERE / "scripts" / "_run_gibbs.py"

# -- Colour palette ------------------------------------------------------------

_FALLBACK_COLORS = [
    "#E07B7B", "#A0785A", "#F0C040", "#5DB0A0", "#9B6BB5",
    "#C8A06E", "#4A90D9", "#90C8F0", "#2060A0", "#60C8C0",
    "#E8A060", "#A0D090", "#D070A0", "#70B0D0",
]

GROUP_COLORS = {
    "pig": "#E07B7B",
    "cow": "#A0785A",
    "chicken": "#F0C040",
    "duck": "#5DB0A0",
    "bat": "#9B6BB5",
    "horse": "#C8A06E",
    "human": "#4A90D9",
    "groundwater": "#90C8F0",
    "seawater": "#2060A0",
    "river": "#60C8C0",
    "Unknown": "#C8C8C8",
}


# -- Data loaders --------------------------------------------------------------

def load_csv_gz(path) -> pd.DataFrame:
    """Load a gzipped CSV with first column as index."""
    if hasattr(path, "seek"):
        path.seek(0)
        data = path.read()
        return pd.read_csv(io.BytesIO(data), compression="gzip", index_col=0)
    with open(path, "rb") as f:
        data = f.read()
    return pd.read_csv(io.BytesIO(data), compression="gzip", index_col=0)


def load_fasta(path) -> dict:
    """Load a FASTA file. Returns {id: sequence}."""
    if hasattr(path, "seek"):
        path.seek(0)
        text = path.read().decode("utf-8", errors="replace")
    else:
        with open(path) as f:
            text = f.read()
    seqs, cur_id, cur_seq = {}, None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if cur_id:
                seqs[cur_id] = "".join(cur_seq)
            cur_id, cur_seq = line[1:].split()[0], []
        elif cur_id:
            cur_seq.append(line)
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


def load_design(path=None) -> dict:
    """Load design file (Sample<tab>Group). Returns {sample: group}."""
    if path is None:
        path = SOURCE_DESIGN
    if hasattr(path, "seek"):
        path.seek(0)
        text = path.read().decode("utf-8", errors="replace")
    else:
        with open(path) as f:
            text = f.read()
    design = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.lower().startswith("sample"):
            continue
        parts = s.split("\t")
        if len(parts) >= 2:
            design[parts[0]] = parts[1]
    return design


# -- Feature alignment ---------------------------------------------------------

def _is_sequence_index(index: pd.Index) -> bool:
    return bool(re.match(r"^[ACGTacgt]{50,}$", str(index[0])))


def _vsearch_map(query_fa: str, db_fa: str, identity: float,
                 threads: int = 16) -> dict:
    uc = query_fa + ".uc"
    subprocess.run(
        ["vsearch", "--usearch_global", query_fa, "--db", db_fa,
         "--id", str(identity), "--uc", uc,
         "--threads", str(threads), "--strand", "plus",
         "--maxaccepts", "1", "--maxrejects", "32", "--quiet"],
        capture_output=True,
    )
    hits = {}
    if os.path.exists(uc):
        with open(uc) as f:
            for line in f:
                p = line.split("\t")
                if p[0] == "H":
                    hits[p[8]] = p[9].strip()
    return hits


def align_features(sink_df, source_df, db_fasta, mode="asv", threads=16):
    """Map sink features to source ASV IDs.

    When sink sequences are longer than source (e.g. V3-V4 sink vs V4 source),
    uses source as query against sink as database so vsearch can find the
    shorter source region within the longer sink amplicon.

    Returns (sink_aligned, source_aligned, (n_matched, n_total)).
    """
    sink_seqbased = _is_sequence_index(sink_df.index)

    with tempfile.TemporaryDirectory() as tmp:
        if sink_seqbased:
            # Write sink sequences
            sq_map = {}
            sink_fa = os.path.join(tmp, "sink.fasta")
            with open(sink_fa, "w") as f:
                for i, seq in enumerate(sink_df.index):
                    sid = f"sq{i}"
                    sq_map[sid] = str(seq)
                    f.write(f">{sid}\n{seq}\n")

            # Write source sequences (only those in source_df)
            src_fa = os.path.join(tmp, "source.fasta")
            with open(src_fa, "w") as f:
                for asv_id, seq in db_fasta.items():
                    if asv_id in source_df.index:
                        f.write(f">{asv_id}\n{seq}\n")

            # Detect length ratio to choose alignment direction
            sink_lens = [len(str(s)) for s in sink_df.index[:20]]
            src_lens = [len(s) for s in db_fasta.values()]
            avg_sink = sum(sink_lens) / len(sink_lens) if sink_lens else 0
            avg_src = sum(src_lens) / len(src_lens) if src_lens else 0

            identity = 1.0 if mode == "asv" else 0.99

            if avg_sink > avg_src * 1.3:
                # Sink seqs are significantly longer than source (e.g. V3-V4 vs V4).
                # Use source as query → sink as db so the shorter source
                # sequence aligns within the longer sink amplicon.
                # Relax identity to 0.97 to tolerate minor boundary differences.
                search_id = min(identity, 0.97)
                raw_hits = _vsearch_map(src_fa, sink_fa, search_id, threads)
                # raw_hits: {source_asv_id: sink_sq_id}
                # Reverse into feature_map: {sink_sequence: source_asv_id}
                feature_map = {}
                for asv_id, sq_id in raw_hits.items():
                    sink_seq = sq_map.get(sq_id)
                    if sink_seq and sink_seq not in feature_map:
                        feature_map[sink_seq] = asv_id
            else:
                # Same-length or sink shorter: original direction
                raw_hits = _vsearch_map(sink_fa, src_fa, identity, threads)
                feature_map = {sq_map[sq_id]: asv_id
                               for sq_id, asv_id in raw_hits.items()}
        else:
            common_ids = sink_df.index.intersection(source_df.index)
            feature_map = {fid: fid for fid in common_ids}

    if not feature_map:
        raise ValueError("No features matched between sink and source tables.")

    matched_source_ids = sorted(set(feature_map.values()))
    source_aligned = source_df.loc[
        source_df.index.isin(matched_source_ids)
    ].T

    new_index = [feature_map.get(str(f)) for f in sink_df.index]
    sink_remapped = sink_df.copy()
    sink_remapped.index = new_index
    sink_remapped = sink_remapped[sink_remapped.index.notna()]
    sink_remapped = sink_remapped.groupby(sink_remapped.index).sum()
    sink_aligned = sink_remapped.T

    common_cols = source_aligned.columns.intersection(sink_aligned.columns)
    source_aligned = source_aligned[common_cols].fillna(0).astype(int)
    sink_aligned = sink_aligned[common_cols].fillna(0).astype(int)

    return sink_aligned, source_aligned, (len(common_cols), len(sink_df.index))


def collapse_by_group(source_aligned: pd.DataFrame, design: dict) -> pd.DataFrame:
    """Sum source counts within each group -> one row per group."""
    cols = source_aligned.columns
    groups = defaultdict(lambda: np.zeros(len(cols), dtype=int))
    for sample in source_aligned.index:
        grp = design.get(sample)
        if grp:
            groups[grp] += source_aligned.loc[sample].values.astype(int)
    if not groups:
        raise ValueError("No source samples matched design entries.")
    return pd.DataFrame(
        {g: v for g, v in groups.items()}, index=cols
    ).T


# -- Gibbs sampling via ST conda env ------------------------------------------

def run_gibbs_subprocess(source_collapsed, sink_aligned,
                         src_depth, snk_depth,
                         restarts, burnin, draws,
                         status_fn=None):
    """Run _run_gibbs.py via ST conda env. Returns (proportions_df, stds_df)."""
    with tempfile.TemporaryDirectory() as tmp:
        src_csv = os.path.join(tmp, "sources.csv")
        snk_csv = os.path.join(tmp, "sinks.csv")
        out_dir = os.path.join(tmp, "results")

        source_collapsed.to_csv(src_csv)
        sink_aligned.to_csv(snk_csv)

        cmd = [
            "conda", "run", "--no-capture-output", "-n", ST_CONDA_ENV,
            "python3", str(_GIBBS_SCRIPT),
            src_csv, snk_csv, out_dir,
            str(src_depth), str(snk_depth),
            str(restarts), str(burnin), str(draws),
        ]

        if status_fn:
            status_fn("Running Gibbs sampling (ST conda env)...")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.stdout and status_fn:
            for line in result.stdout.strip().splitlines():
                status_fn(line)

        if result.returncode != 0:
            raise RuntimeError(
                f"Gibbs sampler failed:\n{result.stdout}\n{result.stderr}"
            )

        mp = pd.read_csv(os.path.join(out_dir, "proportions.csv"), index_col=0)
        mps = pd.read_csv(os.path.join(out_dir, "stds.csv"), index_col=0)

    return mp, mps
