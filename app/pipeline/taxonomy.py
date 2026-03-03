"""
MST-Pipeline — SILVA taxonomy classification via vsearch.
"""
import gzip
import os
import subprocess
import tempfile

from app.config import SILVA_TRAIN_SET

# Module-level singleton for decompressed SILVA path
_silva_decompressed: str | None = None


def _get_silva_path() -> str:
    """Decompress SILVA gz once; cache the path for the process lifetime."""
    global _silva_decompressed
    if _silva_decompressed and os.path.exists(_silva_decompressed):
        return _silva_decompressed
    tmp = tempfile.mktemp(suffix=".fasta", prefix="silva_mst_")
    with gzip.open(str(SILVA_TRAIN_SET), "rt") as gz_in, open(tmp, "w") as out:
        out.write(gz_in.read())
    _silva_decompressed = tmp
    return tmp


def classify(sequences: list[str], identity: float = 0.97,
             threads: int = 4) -> dict[str, str | None]:
    """Classify sequences against SILVA using vsearch.

    SILVA headers are 'Bacteria;Phylum;...;Genus;' — the last non-empty
    field is the genus.

    Returns {sequence: genus_str | None}.
    """
    silva_path = _get_silva_path()

    with tempfile.TemporaryDirectory() as tmp:
        q_fa = os.path.join(tmp, "query.fasta")
        id_to_seq = {}
        with open(q_fa, "w") as fh:
            for i, seq in enumerate(sequences):
                sid = f"sq{i}"
                id_to_seq[sid] = seq
                fh.write(f">{sid}\n{seq.upper()}\n")

        uc = os.path.join(tmp, "out.uc")
        subprocess.run(
            ["vsearch", "--usearch_global", q_fa,
             "--db", silva_path,
             "--id", str(identity),
             "--uc", uc,
             "--threads", str(threads),
             "--strand", "plus",
             "--maxaccepts", "1",
             "--maxrejects", "32",
             "--quiet"],
            capture_output=True,
        )

        genus_map = {}
        if os.path.exists(uc):
            with open(uc) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 10 or parts[0] != "H":
                        continue
                    sid, tax_str = parts[8], parts[9]
                    taxa = [t.strip() for t in tax_str.split(";") if t.strip()]
                    genus_map[sid] = taxa[-1] if taxa else None

    return {id_to_seq[sid]: g for sid, g in genus_map.items()}


def classify_full_taxonomy(sequences: list[str], identity: float = 0.97,
                           threads: int = 4) -> dict[str, list[str]]:
    """Classify sequences and return full taxonomy ranks.

    Returns {sequence: [Domain, Phylum, Class, Order, Family, Genus]}.
    """
    silva_path = _get_silva_path()

    with tempfile.TemporaryDirectory() as tmp:
        q_fa = os.path.join(tmp, "query.fasta")
        id_to_seq = {}
        with open(q_fa, "w") as fh:
            for i, seq in enumerate(sequences):
                sid = f"sq{i}"
                id_to_seq[sid] = seq
                fh.write(f">{sid}\n{seq.upper()}\n")

        uc = os.path.join(tmp, "out.uc")
        subprocess.run(
            ["vsearch", "--usearch_global", q_fa,
             "--db", silva_path,
             "--id", str(identity),
             "--uc", uc,
             "--threads", str(threads),
             "--strand", "plus",
             "--maxaccepts", "1",
             "--maxrejects", "32",
             "--quiet"],
            capture_output=True,
        )

        tax_map = {}
        if os.path.exists(uc):
            with open(uc) as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 10 or parts[0] != "H":
                        continue
                    sid, tax_str = parts[8], parts[9]
                    taxa = [t.strip() for t in tax_str.split(";") if t.strip()]
                    tax_map[sid] = taxa

    return {id_to_seq[sid]: ranks for sid, ranks in tax_map.items()}
