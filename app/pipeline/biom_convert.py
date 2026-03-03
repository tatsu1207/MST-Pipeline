"""
MST-Pipeline — BIOM HDF5 conversion utilities.

Converts DADA2 CSV.gz ASV tables to/from HDF5 BIOM format with sequences
stored as observation metadata.
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from biom import Table


def tsv_to_biom(
    tsv_path: Path,
    output_dir: Path,
    logger: logging.Logger,
) -> Path:
    """Convert a DADA2 ASV table (TSV) to BIOM format.

    The input TSV has columns: ASV_ID, sequence, Sample1, Sample2, ...
    The output BIOM stores ASV_ID as observation IDs, sequences as
    observation metadata, and sample abundances as the data matrix.

    Returns the path to the .biom file.
    """
    df = pd.read_csv(tsv_path, sep="\t")

    asv_ids = df["ASV_ID"].values
    sequences = df["sequence"].values
    sample_columns = [c for c in df.columns if c not in ("ASV_ID", "sequence")]
    data = df[sample_columns].values.astype(np.float64)

    obs_metadata = [{"sequence": seq} for seq in sequences]

    table = Table(
        data,
        observation_ids=asv_ids,
        sample_ids=sample_columns,
        observation_metadata=obs_metadata,
        type="OTU table",
    )

    import h5py
    biom_path = output_dir / "asv_table.biom"
    with h5py.File(biom_path, "w") as f:
        table.to_hdf5(f, generated_by="MST-Pipeline-DADA2")

    logger.info(
        f"BIOM table written: {biom_path} "
        f"({len(asv_ids)} ASVs x {len(sample_columns)} samples)"
    )
    return biom_path


def biom_to_dataframe(biom_path: Path) -> pd.DataFrame:
    """Load a BIOM file back to a features x samples DataFrame."""
    from biom import load_table

    table = load_table(str(biom_path))
    df = pd.DataFrame(
        table.to_dataframe(dense=True).values,
        index=table.ids(axis="observation"),
        columns=table.ids(axis="sample"),
    )
    return df


def extract_sequences_from_biom(biom_path: Path) -> dict:
    """Extract ASV_ID -> sequence mapping from BIOM observation metadata."""
    from biom import load_table

    table = load_table(str(biom_path))
    seqs = {}
    for obs_id in table.ids(axis="observation"):
        md = table.metadata(obs_id, axis="observation")
        if md and "sequence" in md:
            seqs[obs_id] = md["sequence"]
    return seqs


def extract_v4_from_biom(biom_path: Path, logger: logging.Logger | None = None) -> Path | None:
    """Extract V4-region ASVs from a BIOM file using primer matching.

    Filters ASVs whose sequences start with the V4 forward primer (515F).
    Returns path to the filtered BIOM file, or None if no V4 ASVs found.
    """
    from biom import load_table
    from app.pipeline.detect import REGION_PRIMERS, _primer_matches

    table = load_table(str(biom_path))
    v4_primer = REGION_PRIMERS["V4"]["forward"]

    keep_ids = []
    for obs_id in table.ids(axis="observation"):
        md = table.metadata(obs_id, axis="observation")
        if not md or "sequence" not in md:
            continue
        seq = md["sequence"]
        if _primer_matches(seq, v4_primer, max_mismatches=3):
            keep_ids.append(obs_id)

    if not keep_ids:
        if logger:
            logger.info("No V4 ASVs found; skipping V4 extraction.")
        return None

    filtered = table.filter(keep_ids, axis="observation", inplace=False)

    import h5py
    v4_path = biom_path.parent / "asv_table.biom"
    with h5py.File(v4_path, "w") as f:
        filtered.to_hdf5(f, generated_by="MST-Pipeline-V4-Extract")

    if logger:
        logger.info(
            f"V4 extraction: kept {len(keep_ids)}/{len(table.ids(axis='observation'))} ASVs"
        )
    return v4_path


def biom_to_csv_gz(biom_path: Path, output_path: Path | None = None) -> Path:
    """Convert a BIOM file back to CSV.gz format with Sequence column."""
    from biom import load_table

    table = load_table(str(biom_path))
    obs_ids = table.ids(axis="observation")
    sample_ids = table.ids(axis="sample")

    data = {}
    sequences = []
    for obs_id in obs_ids:
        md = table.metadata(obs_id, axis="observation")
        seq = md.get("sequence", obs_id) if md else obs_id
        sequences.append(seq)

    data["Sequence"] = sequences
    for sid in sample_ids:
        data[sid] = table.data(sid, axis="sample", dense=True).astype(int)

    df = pd.DataFrame(data)
    if output_path is None:
        output_path = biom_path.parent / "asv_table.csv.gz"
    df.to_csv(output_path, index=False, compression="gzip")
    return output_path
