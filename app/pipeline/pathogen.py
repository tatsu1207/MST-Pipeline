"""
MST-Pipeline — WHO pathogen dictionary and detection logic.

Extracted from MST/scripts/app_pathogen_ui.py: pure functions only.
"""
import pandas as pd

# -- Pathogen dictionary -------------------------------------------------------

PATHOGENS = {
    # WHO Critical Priority
    "Acinetobacter": "critical",
    "Pseudomonas": "critical",
    "Klebsiella": "critical",
    "Escherichia": "critical",
    "Enterobacter": "critical",
    "Serratia": "critical",
    "Proteus": "critical",
    "Morganella": "critical",
    "Providencia": "critical",
    # WHO High Priority
    "Enterococcus": "high",
    "Staphylococcus": "high",
    "Helicobacter": "high",
    "Campylobacter": "high",
    "Salmonella": "high",
    "Neisseria": "high",
    # WHO Medium Priority
    "Streptococcus": "medium",
    "Haemophilus": "medium",
    "Shigella": "medium",
    # Other notable pathogens
    "Clostridioides": "other",
    "Clostridium": "other",
    "Vibrio": "other",
    "Listeria": "other",
    "Yersinia": "other",
    "Legionella": "other",
    "Brucella": "other",
    "Francisella": "other",
    "Bordetella": "other",
    "Mycobacterium": "other",
    "Arcobacter": "other",
    "Aliarcobacter": "other",
    # Waterborne & fecal-associated pathogens (MST-relevant)
    "Aeromonas": "other",
    "Plesiomonas": "other",
    "Leptospira": "other",
    "Bacteroides": "other",
    "Treponema": "other",
    "Prevotella": "other",
    "Fusobacterium": "other",
    "Veillonella": "other",
    "Edwardsiella": "other",
    "Citrobacter": "other",
    "Hafnia": "other",
    "Cronobacter": "other",
}

_PRIORITY_ORDER = ["critical", "high", "medium", "other"]
_PRIORITY_LABELS = {
    "critical": "WHO Critical",
    "high": "WHO High",
    "medium": "WHO Medium",
    "other": "Other Notable",
}


def detect(asv_df: pd.DataFrame, seq_to_genus: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Group ASV counts by pathogenic genus, compute relative abundance.

    Args:
        asv_df: features x samples DataFrame (index = sequences)
        seq_to_genus: {sequence: genus_str | None} from taxonomy classification

    Returns:
        (count_df, ra_df) — both indexed samples x genera.
        Returns empty DataFrames if no pathogens found.
    """
    raw_genera = pd.Series(
        [seq_to_genus.get(s) for s in asv_df.index], index=asv_df.index
    )
    # SILVA uses qualified names like "Clostridium sensu stricto 1";
    # match on the first word so they still hit the pathogen list.
    genera = raw_genera.map(
        lambda g: g.split()[0] if isinstance(g, str) and g.split()[0] in PATHOGENS else g
    )
    mask = genera.isin(PATHOGENS)
    if not mask.any():
        return pd.DataFrame(), pd.DataFrame()

    grouped = asv_df.loc[mask].copy()
    grouped.index = genera[mask]
    grouped = grouped.groupby(grouped.index).sum()  # genera x samples

    total = asv_df.sum(axis=0)
    ra = grouped.div(total, axis=1) * 100

    return grouped.T, ra.T  # samples x genera


def genus_colors(genera: list[str]) -> dict[str, str]:
    """Assign colors to genera based on their priority group.

    Uses Plotly-compatible hex colors.
    """
    _COLOR_FAMILIES = {
        "critical": ["#ef5350", "#e53935", "#c62828", "#b71c1c", "#d32f2f",
                      "#f44336", "#e57373", "#ff8a80", "#ff5252"],
        "high": ["#ff9800", "#f57c00", "#ef6c00", "#e65100", "#fb8c00", "#ffa726"],
        "medium": ["#fdd835", "#f9a825", "#f57f17"],
        "other": ["#42a5f5", "#1e88e5", "#1565c0", "#0d47a1", "#2196f3",
                   "#64b5f6", "#90caf9", "#bbdefb", "#1976d2", "#0277bd", "#01579b"],
    }

    by_priority = {}
    for g in genera:
        by_priority.setdefault(PATHOGENS.get(g, "other"), []).append(g)

    colors = {}
    for p, gens in by_priority.items():
        family = _COLOR_FAMILIES.get(p, _COLOR_FAMILIES["other"])
        for i, g in enumerate(gens):
            colors[g] = family[i % len(family)]
    return colors


def build_summary(ra_df: pd.DataFrame) -> list[dict]:
    """Build a detection summary table from relative abundance data."""
    summary = []
    for genus in sorted(
        ra_df.columns,
        key=lambda g: (_PRIORITY_ORDER.index(PATHOGENS.get(g, "other")), g),
    ):
        p = PATHOGENS.get(genus, "other")
        summary.append({
            "Genus": genus,
            "Priority": _PRIORITY_LABELS[p],
            "Max RA (%)": round(float(ra_df[genus].max()), 4),
            "Mean RA (%)": round(float(ra_df[genus].mean()), 4),
            "Samples detected": int((ra_df[genus] > 0).sum()),
        })
    return summary
