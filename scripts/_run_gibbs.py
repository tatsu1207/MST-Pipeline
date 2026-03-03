#!/usr/bin/env python3
"""SourceTracker2 Gibbs runner — executed inside the ST conda environment.

Usage:
    python3 _run_gibbs.py <source_csv> <sink_csv> <output_dir>
                          <src_depth> <snk_depth> <restarts> <burnin> <draws>

source_csv  : rows = source groups (already collapsed), cols = features
sink_csv    : rows = sink samples,                      cols = features
output_dir  : proportions.csv and stds.csv written here
"""

import sys
import numpy as np
import pandas as pd

# sourcetracker uses removed np.int / np.float aliases (dropped in NumPy 1.24)
np.int   = int
np.float = float
np.bool  = bool

from sourcetracker._sourcetracker import _gibbs, subsample_counts


def subsample_dataframe(df, depth):
    """Row-wise rarefaction wrapper around subsample_counts."""
    rows = {
        idx: subsample_counts(row.values.astype(int), depth)
        for idx, row in df.iterrows()
    }
    return pd.DataFrame(rows, index=df.columns).T.rename_axis(df.index.name)

source_path, sink_path, output_dir = sys.argv[1], sys.argv[2], sys.argv[3]
src_depth = int(sys.argv[4])
snk_depth = int(sys.argv[5])
restarts  = int(sys.argv[6])
burnin    = int(sys.argv[7])
draws     = int(sys.argv[8])

source_df = pd.read_csv(source_path, index_col=0)
sink_df   = pd.read_csv(sink_path,   index_col=0)

# Drop source groups with zero aligned reads before computing depth caps —
# keeping them causes subsample_dataframe to crash when src_depth is 0.
zero_src = source_df.index[source_df.sum(axis=1) == 0].tolist()
if zero_src:
    print(f"WARNING: dropping source groups with 0 aligned reads: {zero_src}", flush=True)
    source_df = source_df[source_df.sum(axis=1) > 0]

if source_df.empty:
    print("ERROR: no source groups have any aligned reads.", flush=True)
    sys.exit(1)

# Cap rarefaction depths to the actual minimum available counts
# (aligned feature counts are lower than the original table totals)
src_min = int(source_df.sum(axis=1).min())
snk_min = int(sink_df.sum(axis=1).min())

if src_depth > src_min:
    print(f"WARNING: source depth {src_depth} > min available {src_min}; "
          f"capping to {src_min}", flush=True)
    src_depth = src_min

if snk_depth > snk_min:
    print(f"WARNING: sink depth {snk_depth} > min available {snk_min}; "
          f"capping to {snk_min}", flush=True)
    snk_depth = snk_min

if src_depth <= 0:
    print("ERROR: source depth is 0 after capping — no aligned reads in source groups.", flush=True)
    sys.exit(1)
if snk_depth <= 0:
    print("ERROR: sink depth is 0 after capping — no aligned reads in sink samples. "
          "This usually means no ASV sequences matched the source reference database.", flush=True)
    sys.exit(1)

print(f"Sources: {source_df.shape}  (depth→{src_depth})", flush=True)
print(f"Sinks:   {sink_df.shape}  (depth→{snk_depth})", flush=True)

# Drop rows below the rarefaction depth before calling subsample_dataframe,
# which would crash on under-depth rows instead of skipping them.
source_df = source_df[source_df.sum(axis=1) >= src_depth]
sink_df   = sink_df[sink_df.sum(axis=1)   >= snk_depth]

source_rare = subsample_dataframe(source_df, src_depth)
sink_rare   = subsample_dataframe(sink_df,   snk_depth)

dropped_src = set(pd.read_csv(source_path, index_col=0).index) - set(source_rare.index)
dropped_snk = set(pd.read_csv(sink_path,   index_col=0).index) - set(sink_rare.index)
if dropped_src:
    print(f"WARNING: dropped source groups (too few reads): {dropped_src}", flush=True)
if dropped_snk:
    print(f"WARNING: dropped sink samples (too few reads): {dropped_snk}", flush=True)

if source_rare.empty:
    print("ERROR: no source groups survived rarefaction.", flush=True)
    sys.exit(1)
if sink_rare.empty:
    print("ERROR: no sink samples survived rarefaction.", flush=True)
    sys.exit(1)

# Align columns (should already match, but just in case)
common = source_rare.columns.intersection(sink_rare.columns)
source_rare = source_rare[common]
sink_rare   = sink_rare[common]

print(f"Rarefied: {source_rare.shape[0]} sources, "
      f"{sink_rare.shape[0]} sinks, {len(common)} features", flush=True)
print(f"Running Gibbs: restarts={restarts} burnin={burnin} draws={draws}", flush=True)

mpm, mps = _gibbs(
    source_rare, sink_rare,
    alpha1=0.001, alpha2=0.1, beta=10,
    restarts=restarts, draws_per_restart=draws,
    burnin=burnin, delay=2,  # must be >1; delay=1 means (x%1)==1 is never true
)

import os
os.makedirs(output_dir, exist_ok=True)
mpm.to_csv(f"{output_dir}/proportions.csv")
mps.to_csv(f"{output_dir}/stds.csv")
print("Done.", flush=True)
