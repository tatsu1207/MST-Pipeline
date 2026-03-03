"""
MST-Pipeline — Unified single-page layout.

Combines File Manager (upload/register), DADA2 pipeline, Source Tracking,
and Taxonomy/Pathogen detection into one top-to-bottom workflow.
"""
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path

import dash_bootstrap_components as dbc
import dash_uploader as du
import plotly.graph_objects as go
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update, dash_table

from app.config import UPLOAD_DIR, to_absolute, to_relative
from app.dashboard.app import app as dash_app
from app.db.database import SessionLocal
from app.db.models import FastqFile, Upload
from app.pipeline.detect import (
    REGION_PRIMERS,
    _primer_matches,
    _read_fastq_sequences,
    detect_sequencing_type,
    detect_variable_region,
    extract_sample_name,
)

MAX_CPUS = os.cpu_count() or 1

STEPS = ["cutadapt", "auto_trunc", "dada2", "taxonomy", "split_samples"]
STEP_LABELS = {
    "cutadapt": "Cutadapt",
    "auto_trunc": "Auto-Trunc",
    "dada2": "DADA2",
    "taxonomy": "Taxonomy",
    "split_samples": "Split Samples",
}


# ══════════════════════════════════════════════════════════════════════════════
# Layout
# ══════════════════════════════════════════════════════════════════════════════


def get_layout():
    """Build the unified single-page layout dynamically."""
    from app.pipeline.sourcetracker import load_design, GROUP_COLORS

    # -- Source DB stats ----------------------------------------------------------
    source_db_stats = _get_source_db_stats()

    # -- Source group checkboxes -----------------------------------------------
    design = load_design()
    groups = sorted(set(design.values()))

    group_checkboxes = []
    for g in groups:
        color = GROUP_COLORS.get(g, "#888888")
        group_checkboxes.append(
            dbc.Checklist(
                options=[{"label": html.Span([
                    html.Span(style={"display": "inline-block", "width": "12px", "height": "12px",
                                     "backgroundColor": color, "borderRadius": "50%",
                                     "marginRight": "6px"}),
                    g
                ]), "value": g}],
                value=[g],
                id={"type": "group-check", "index": g},
                inline=True,
                className="me-3 mb-2",
            )
        )

    return dbc.Container(
        [
            html.H3("MST Pipeline", className="mb-4"),

            # ── Section A: Upload FASTQ Files ────────────────────────────────
            dbc.Card(
                [
                    dbc.CardHeader(html.H5("Upload FASTQ Files", className="mb-0")),
                    dbc.CardBody(
                        [
                            du.Upload(
                                id="fm-du-upload",
                                text="Drag & drop FASTQ files here or click to browse",
                                text_completed="Uploaded: ",
                                max_file_size=5120,
                                chunk_size=10,
                                max_files=100,
                                filetypes=["gz", "fastq", "fq"],
                                cancel_button=True,
                                pause_button=True,
                            ),
                            html.Div(id="fm-upload-file-status", className="mt-2"),
                            dbc.Button(
                                "Register Upload",
                                id="fm-btn-register",
                                color="success",
                                className="mt-2",
                                disabled=True,
                            ),
                            html.Div(id="fm-upload-status", className="mt-3"),
                        ]
                    ),
                ],
                className="mb-4",
            ),
            # Hidden stores for upload/file management
            dcc.Store(id="fm-staging-dir", data=None),
            dcc.Store(id="fm-upload-trigger", data=0),
            dcc.Store(id="fm-sort-field", data="sample_name"),
            dcc.Store(id="fm-sort-asc", data=True),
            dcc.Store(id="fm-checked-samples", data=[]),
            dcc.ConfirmDialog(id="fm-confirm-delete", message=""),

            # ── Section A2: Import BIOM File ─────────────────────────────────
            dbc.Card(
                [
                    dbc.CardHeader(html.H5("Import Pre-processed BIOM File", className="mb-0")),
                    dbc.CardBody(
                        [
                            html.P(
                                "Upload a BIOM file from an external DADA2 run to skip "
                                "directly to SourceTracker and Pathogen detection.",
                                className="text-muted small",
                            ),
                            dcc.Upload(
                                id="biom-upload",
                                children=html.Div([
                                    "Drag & drop a ",
                                    html.B(".biom"),
                                    " file here, or ",
                                    html.A("click to browse"),
                                ]),
                                style={
                                    "width": "100%",
                                    "height": "60px",
                                    "lineHeight": "60px",
                                    "borderWidth": "1px",
                                    "borderStyle": "dashed",
                                    "borderRadius": "5px",
                                    "textAlign": "center",
                                    "cursor": "pointer",
                                },
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        dbc.Input(
                                            id="biom-name",
                                            placeholder="Dataset name (optional)",
                                            value="",
                                            size="sm",
                                        ),
                                        width=4,
                                        className="mt-2",
                                    ),
                                ],
                            ),
                            html.Div(id="biom-import-status", className="mt-2"),
                        ]
                    ),
                ],
                className="mb-4",
            ),
            dcc.Store(id="biom-taxonomy-dataset-id"),

            # ── Section B: Registered Samples + DADA2 ────────────────────────
            dbc.Card(
                [
                    dbc.CardHeader(
                        dbc.Row(
                            [
                                dbc.Col(html.H5("Registered Samples", className="mb-0"), width="auto"),
                                dbc.Col(
                                    dbc.Button(
                                        "Delete Selected",
                                        id="fm-btn-delete-selected",
                                        color="danger",
                                        size="sm",
                                        disabled=True,
                                        className="me-2",
                                    ),
                                    width="auto",
                                    className="ms-auto",
                                ),
                                dbc.Col(
                                    dbc.ButtonGroup(
                                        [
                                            dbc.Button("Sample", id="fm-sort-sample", color="primary", size="sm", outline=True),
                                            dbc.Button("Type", id="fm-sort-type", color="primary", size="sm", outline=True),
                                            dbc.Button("Reads", id="fm-sort-reads", color="primary", size="sm", outline=True),
                                            dbc.Button("Region", id="fm-sort-region", color="primary", size="sm", outline=True),
                                            dbc.Button("Read Len", id="fm-sort-readlen", color="primary", size="sm", outline=True),
                                            dbc.Button("Date", id="fm-sort-date", color="primary", size="sm", outline=True),
                                        ],
                                        size="sm",
                                    ),
                                    width="auto",
                                ),
                            ],
                            align="center",
                            className="g-0",
                        ),
                    ),
                    dbc.CardBody(
                        [
                            dbc.Row(
                                [
                                    dbc.Col(dbc.Input(id="fm-f-sample", placeholder="Sample", debounce=True, size="sm"), width=2),
                                    dbc.Col(dbc.Input(id="fm-f-type", placeholder="Type", debounce=True, size="sm"), width=1),
                                    dbc.Col(width=1),
                                    dbc.Col(dbc.Input(id="fm-f-region", placeholder="Region", debounce=True, size="sm"), width=1),
                                    dbc.Col(width=1),
                                    dbc.Col(width=1),
                                ],
                                className="mb-2 g-1",
                            ),
                            html.Div(id="fm-files-table", children="No files registered yet."),
                        ]
                    ),
                    # DADA2 controls within same card
                    dbc.CardFooter(
                        [
                            html.Div(id="sample-selection-summary", className="mb-2"),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dbc.Label("CPU Threads", className="fw-bold"),
                                            dbc.InputGroup(
                                                [
                                                    dbc.Input(
                                                        id="input-threads",
                                                        type="text",
                                                        value="default",
                                                        placeholder="default",
                                                        style={"maxWidth": "120px"},
                                                    ),
                                                    dbc.InputGroupText(
                                                        f"default = samples x 2 (max {MAX_CPUS} cores)",
                                                        className="text-muted small",
                                                    ),
                                                ],
                                            ),
                                        ],
                                        width=5,
                                    ),
                                    dbc.Col(
                                        html.Div(
                                            [
                                                dbc.Label("Error Model", className="fw-bold"),
                                                dbc.Select(
                                                    id="input-error-model",
                                                    options=[
                                                        {"label": "Default (short reads)", "value": "default"},
                                                        {"label": "PacBio HiFi / CCS", "value": "PacBio"},
                                                    ],
                                                    value="default",
                                                ),
                                            ],
                                            id="error-model-container",
                                            style={"display": "none"},
                                        ),
                                        width=3,
                                    ),
                                    dbc.Col(
                                        dbc.Button(
                                            "Start DADA2",
                                            id="btn-start-pipeline",
                                            color="success",
                                            disabled=True,
                                            className="mt-4",
                                        ),
                                        width="auto",
                                        className="ms-auto",
                                    ),
                                ],
                                align="end",
                            ),
                            html.Div(id="pipeline-launch-error", className="mt-2"),
                            html.Div(id="pipeline-detection-summary"),
                            # DADA2 progress section
                            html.Div(
                                id="pipeline-progress-section",
                                children=[
                                    html.Div(
                                        [
                                            html.H5("Progress", className="d-inline"),
                                            dbc.Button(
                                                "Cancel Pipeline",
                                                id="btn-cancel-pipeline",
                                                color="danger",
                                                size="sm",
                                                className="ms-3",
                                            ),
                                        ],
                                        className="d-flex align-items-center mb-2",
                                    ),
                                    html.Div(id="pipeline-status-badge", className="mb-3"),
                                    html.Div(
                                        dbc.Progress(
                                            id="pipeline-progress-bar",
                                            value=0,
                                            striped=True,
                                            animated=True,
                                            className="mb-3",
                                        ),
                                        id="pipeline-progress-bar-wrapper",
                                    ),
                                    html.Div(id="pipeline-steps-display", className="mb-3"),
                                    html.H6("Pipeline Log"),
                                    html.Pre(
                                        id="pipeline-log-viewer",
                                        style={
                                            "maxHeight": "400px",
                                            "overflowY": "auto",
                                            "backgroundColor": "#f5f5f5",
                                            "color": "#333",
                                            "padding": "10px",
                                            "borderRadius": "5px",
                                            "fontSize": "0.8rem",
                                            "whiteSpace": "pre-wrap",
                                        },
                                    ),
                                    html.Div(
                                        id="pipeline-download-section",
                                        style={"display": "none"},
                                    ),
                                ],
                                style={"display": "none"},
                            ),
                        ]
                    ),
                ],
                className="mb-4",
            ),
            # Hidden stores for DADA2
            dcc.Store(id="pipeline-active-dataset", storage_type="session"),
            dcc.Store(id="file-ids-map-store"),
            dcc.Interval(id="pipeline-poll", interval=3000, disabled=True, n_intervals=0),

            # ── Section D: Pipeline History ───────────────────────────────────
            dbc.Card(
                [
                    dbc.CardHeader(html.H5("Pipeline History", className="mb-0")),
                    dbc.CardBody(
                        html.Div(
                            id="pipeline-history-table",
                            style={"maxHeight": "300px", "overflowY": "auto"},
                        ),
                    ),
                ],
                className="mb-4",
            ),
            dcc.Download(id="download-history-log"),
            dcc.Download(id="download-sample-asv"),

            # ── Section C: Analysis (Source Tracking + Pathogen) ──────────────
            dbc.Card(
                [
                    dbc.CardHeader(html.H5("Analysis", className="mb-0")),
                    dbc.CardBody(
                        [
                            # Selection summary (populated from history checkboxes)
                            html.Div(
                                id="an-selection-summary",
                                children=html.P(
                                    "Select samples from Pipeline History above.",
                                    className="text-muted",
                                ),
                                className="mb-3",
                            ),
                            # Pathogen Detection (uses pre-computed taxonomy from DADA2)
                            dbc.Card(dbc.CardBody([
                                html.H6("Pathogen Detection"),
                                html.P("Uses taxonomy assigned during DADA2 processing.", className="text-muted small mb-0"),
                            ]), className="mb-3"),
                            # Hidden inputs to keep callbacks happy
                            html.Div([
                                dbc.Input(id="tax-identity", type="hidden", value=97),
                                dbc.Input(id="tax-threads", type="hidden", value=4),
                            ], style={"display": "none"}),
                            # Source Tracking params
                            dbc.Card(dbc.CardBody([
                                html.H6("Source Tracking"),
                                dbc.Row([
                                    dbc.Col(
                                        _build_source_db_info(source_db_stats),
                                        md=6,
                                    ),
                                    dbc.Col([
                                        html.P("Source Groups:", className="text-muted mb-1"),
                                        html.Div(group_checkboxes, className="d-flex flex-wrap mb-2") if group_checkboxes else html.P("No groups found.", className="text-muted"),
                                    ], md=6),
                                ], className="mb-2"),
                                dbc.Row([
                                    dbc.Col([
                                        dbc.Label("Mode", className="small"),
                                        dbc.Select(
                                            id="st-mode",
                                            options=[
                                                {"label": "ASV (100%)", "value": "asv"},
                                                {"label": "OTU (99%)", "value": "otu"},
                                            ],
                                            value="asv",
                                        ),
                                    ], md=2),
                                    dbc.Col([
                                        dbc.Label("Source Depth", className="small"),
                                        dbc.Input(id="st-src-depth", type="number", value=1000, min=100),
                                    ], md=2),
                                    dbc.Col([
                                        dbc.Label("Sink Depth", className="small"),
                                        dbc.Input(id="st-snk-depth", type="number", value=1000, min=100),
                                    ], md=2),
                                    dbc.Col([
                                        dbc.Label("Restarts", className="small"),
                                        dbc.Input(id="st-restarts", type="number", value=10, min=1),
                                    ], md=2),
                                    dbc.Col([
                                        dbc.Label("Burnin", className="small"),
                                        dbc.Input(id="st-burnin", type="number", value=100, min=10),
                                    ], md=1),
                                    dbc.Col([
                                        dbc.Label("Draws", className="small"),
                                        dbc.Input(id="st-draws", type="number", value=1, min=1),
                                    ], md=1),
                                    dbc.Col([
                                        dbc.Label("Threads", className="small"),
                                        dbc.Input(id="st-threads", type="number", value=4, min=1),
                                    ], md=1),
                                ]),
                            ]), className="mb-3"),
                            dbc.Button(
                                "Run Analysis", id="run-analysis-btn",
                                color="primary", size="lg", className="mb-3",
                                disabled=True,
                            ),
                            # Pathogen progress + results
                            html.Div(id="tax-progress-area", className="mb-3"),
                            html.Div(id="pathogen-results-area", className="mb-3"),
                            # Source Tracking progress + results
                            html.Div(id="st-progress-area", className="mb-3"),
                            html.Div(id="st-results-area"),
                        ]
                    ),
                ],
                className="mb-4",
            ),
            # Hidden stores for analysis
            dcc.Store(id="an-checked-history", data=[]),
            dcc.Store(id="st-run-id"),
            dcc.Store(id="tax-run-id"),
            dcc.Interval(id="st-poll", interval=3000, disabled=True),
            dcc.Interval(id="tax-poll", interval=3000, disabled=True),
        ],
        fluid=True,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════


def _human_size(size_mb):
    """Format size in MB to human-readable string."""
    if size_mb is None:
        return "?"
    if size_mb < 1:
        return f"{size_mb * 1024:.0f} KB"
    if size_mb >= 1024:
        return f"{size_mb / 1024:.1f} GB"
    return f"{size_mb:.1f} MB"


def _count_reads(fastq_path: Path) -> int | None:
    import gzip
    try:
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        count = 0
        with opener(fastq_path, "rt") as f:
            for _ in f:
                count += 1
        return count // 4
    except Exception:
        return None


def _avg_read_length(fastq_path: Path, n_reads: int = 200) -> int | None:
    import gzip
    try:
        opener = gzip.open if str(fastq_path).endswith(".gz") else open
        lengths = []
        with opener(fastq_path, "rt") as f:
            line_num = 0
            for line in f:
                line_num += 1
                if line_num % 4 == 2:
                    lengths.append(len(line.strip()))
                    if len(lengths) >= n_reads:
                        break
        return round(sum(lengths) / len(lengths)) if lengths else None
    except Exception:
        return None


def _build_files_table(sort_by="sample_name", filters=None, ascending=True, checked_samples=None):
    """Build a sample-level table of all registered FASTQ files.

    Also returns a file_ids_map {sample_name: [file_id, ...]} for DADA2 launch.
    """
    if filters is None:
        filters = {}
    if checked_samples is None:
        checked_samples = set()
    db = SessionLocal()
    try:
        query = (
            db.query(FastqFile, Upload.created_at, Upload.variable_region, Upload.primers_detected)
            .join(Upload, FastqFile.upload_id == Upload.id)
        )

        results = query.all()
        if not results:
            return "No files registered yet.", {}

        sample_map = defaultdict(lambda: {
            "files": [], "total_size": 0.0, "date": None,
            "region": None, "upload_id": None, "primers_detected": None,
        })
        for f, upload_date, region, primers_detected in results:
            s = sample_map[f.sample_name]
            s["files"].append(f)
            s["total_size"] += f.file_size_mb or 0
            s["upload_id"] = f.upload_id
            if upload_date:
                s["date"] = upload_date
            if region:
                s["region"] = region
            if primers_detected is not None:
                s["primers_detected"] = primers_detected

        sample_rows = []
        file_ids_map = {}
        for sample_name, info in sample_map.items():
            files = sorted(info["files"], key=lambda x: x.filename)
            file_ids_map[sample_name] = [f.id for f in files]
            directions = set(f.read_direction for f in files)
            if "R1" in directions and "R2" in directions:
                direction = "PE"
            elif directions == {"single"}:
                direction = "SE"
            else:
                direction = ", ".join(sorted(directions))

            total_reads = sum(f.read_count or 0 for f in files)
            avg_lengths = [f.avg_read_length for f in files if f.avg_read_length]
            avg_len = round(sum(avg_lengths) / len(avg_lengths)) if avg_lengths else None

            sample_rows.append({
                "sample_name": sample_name,
                "direction": direction,
                "total_reads": total_reads,
                "region": info["region"],
                "avg_len": avg_len,
                "date": info["date"],
                "n_files": len(files),
                "primers_detected": info["primers_detected"],
            })

        # Per-column filters
        def _matches(r):
            if filters.get("sample") and filters["sample"] not in r["sample_name"].lower():
                return False
            if filters.get("type") and filters["type"] not in r["direction"].lower():
                return False
            if filters.get("region") and filters["region"] not in (r["region"] or "").lower():
                return False
            return True

        sample_rows = [r for r in sample_rows if _matches(r)]

        if not sample_rows:
            return html.P("No matching samples.", className="text-muted"), file_ids_map

        # Sort
        reverse = not ascending
        sort_keys = {
            "sample_name": lambda r: r["sample_name"].lower(),
            "direction": lambda r: r["direction"],
            "total_reads": lambda r: r["total_reads"],
            "region": lambda r: r["region"] or "",
            "avg_read_length": lambda r: r["avg_len"] or 0,
            "created_at": lambda r: r["date"] or "",
        }
        numeric_fields = {"total_reads", "avg_read_length", "created_at"}
        key_fn = sort_keys.get(sort_by, sort_keys["sample_name"])
        effective_reverse = (not reverse) if sort_by in numeric_fields else reverse
        sample_rows.sort(key=key_fn, reverse=effective_reverse)

        rows = []
        for r in sample_rows:
            date_str = r["date"].strftime("%Y-%m-%d") if r["date"] else ""
            reads_str = f"{r['total_reads']:,}" if r["total_reads"] else "\u2014"
            is_checked = r["sample_name"] in checked_samples
            pd_val = r["primers_detected"]
            if pd_val is True:
                primers_badge = dbc.Badge("Yes", color="success")
            elif pd_val is False:
                primers_badge = dbc.Badge("No", color="info")
            else:
                primers_badge = dbc.Badge("?", color="secondary")
            rows.append(
                html.Tr(
                    [
                        html.Td(
                            dbc.Checkbox(
                                id={"type": "fm-check", "index": r["sample_name"]},
                                value=is_checked,
                            ),
                            style={"width": "30px"},
                        ),
                        html.Td(r["sample_name"]),
                        html.Td(
                            dbc.Badge(
                                r["direction"],
                                color="primary" if r["direction"] == "PE" else "info",
                            )
                        ),
                        html.Td(reads_str),
                        html.Td(r["region"] or "\u2014", className="small"),
                        html.Td(primers_badge),
                        html.Td(f"{r['avg_len']} bp" if r["avg_len"] else "\u2014", className="small"),
                        html.Td(date_str, className="text-muted small"),
                    ]
                )
            )

        total_files = sum(r["n_files"] for r in sample_rows)
        total_reads = sum(r["total_reads"] for r in sample_rows)

        arrow = " ^" if ascending else " v"
        def _th(label, field):
            suffix = arrow if sort_by == field else ""
            return html.Th(f"{label}{suffix}")

        table_div = html.Div(
            [
                html.P(
                    f"{len(sample_rows)} samples, {total_files} files, {total_reads:,} reads total",
                    className="text-muted mb-2",
                ),
                dbc.Table(
                    [
                        html.Thead(
                            html.Tr(
                                [
                                    html.Th(
                                        dbc.Checkbox(id="fm-check-all", value=False),
                                        style={"width": "30px"},
                                    ),
                                    _th("Sample", "sample_name"),
                                    _th("Type", "direction"),
                                    _th("# Reads", "total_reads"),
                                    _th("Region", "region"),
                                    html.Th("Primers"),
                                    _th("Avg Read Len", "avg_read_length"),
                                    _th("Date", "created_at"),
                                ]
                            )
                        ),
                        html.Tbody(rows),
                    ],
                    bordered=True,
                    color="light",
                    hover=True,
                    size="sm",
                ),
            ]
        )
        return table_div, file_ids_map
    finally:
        db.close()


def _download_filename(ds, prefix: str, ext: str) -> str:
    if not ds:
        return f"{prefix}.{ext}"
    name = re.sub(r"[^\w\-.]", "_", ds.name or prefix).strip("_") or prefix
    region = ds.variable_region or ""
    date = ds.created_at.strftime("%Y-%m-%d") if ds.created_at else ""
    parts = [p for p in (name, region, date) if p]
    return f"{'_'.join(parts)}.{ext}"


def _dir_size(path: Path) -> int:
    if not path.is_dir():
        return 0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())


def _format_size(nbytes: int) -> str:
    if nbytes == 0:
        return "\u2014"
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024:
            return f"{nbytes:.1f} {unit}" if unit != "B" else f"{nbytes} B"
        nbytes /= 1024
    return f"{nbytes:.1f} TB"


def _get_source_db_stats():
    """Compute summary statistics for the source microbiome database."""
    from app.pipeline.sourcetracker import load_csv_gz, load_design, GROUP_COLORS
    from app.config import SOURCE_TABLE, SOURCE_DESIGN

    try:
        source_df = load_csv_gz(str(SOURCE_TABLE))
        design = load_design(str(SOURCE_DESIGN))
    except Exception:
        return None

    # source_df: rows=ASVs, columns=samples
    group_samples = defaultdict(list)
    for sample in source_df.columns:
        grp = design.get(sample)
        if grp:
            group_samples[grp].append(sample)

    stats = []
    for grp in sorted(group_samples.keys()):
        samples = group_samples[grp]
        sub = source_df[samples]
        n_samples = len(samples)
        avg_reads = int(sub.sum(axis=0).mean())
        n_asvs = int((sub.sum(axis=1) > 0).sum())
        color = GROUP_COLORS.get(grp, "#888888")
        stats.append({
            "group": grp,
            "n_samples": n_samples,
            "avg_reads": avg_reads,
            "n_asvs": n_asvs,
            "color": color,
        })

    return {
        "groups": stats,
        "total_asvs": source_df.shape[0],
        "total_samples": sum(len(v) for v in group_samples.values()),
    }


def _build_source_db_info(stats):
    """Build a compact display of source microbiome DB statistics."""
    if not stats:
        return html.P("Source DB not available.", className="text-muted small mb-2")

    rows = []
    for g in stats["groups"]:
        rows.append(
            html.Tr([
                html.Td(
                    html.Span([
                        html.Span(style={
                            "display": "inline-block", "width": "10px", "height": "10px",
                            "backgroundColor": g["color"], "borderRadius": "50%",
                            "marginRight": "5px",
                        }),
                        g["group"],
                    ]),
                    className="small",
                ),
                html.Td(str(g["n_samples"]), className="small text-end"),
                html.Td(f"{g['avg_reads']:,}", className="small text-end"),
                html.Td(f"{g['n_asvs']:,}", className="small text-end"),
            ])
        )
    # Total row
    rows.append(
        html.Tr([
            html.Td(html.Strong("Total"), className="small"),
            html.Td(html.Strong(str(stats["total_samples"])), className="small text-end"),
            html.Td("", className="small"),
            html.Td(html.Strong(f"{stats['total_asvs']:,}"), className="small text-end"),
        ], style={"borderTop": "2px solid #555"})
    )

    return html.Div([
        html.P("Source Microbiome DB:", className="text-muted mb-1 small fw-bold"),
        dbc.Table(
            [
                html.Thead(html.Tr([
                    html.Th("Source", className="small"),
                    html.Th("Samples", className="small text-end"),
                    html.Th("Avg Reads", className="small text-end"),
                    html.Th("ASVs", className="small text-end"),
                ])),
                html.Tbody(rows),
            ],
            bordered=True, size="sm", className="mb-2",
            style={"fontSize": "0.8rem"},
        ),
    ], className="mb-3")


def _build_history_table(checked_history=None):
    from app.db.database import get_session
    from app.db.models import Dataset, Sample
    from app.config import DATASET_DIR

    with get_session() as db:
        datasets = (
            db.query(Dataset)
            .order_by(Dataset.created_at.desc())
            .all()
        )
        ds_ids = [d.id for d in datasets]
        all_samples = []
        if ds_ids:
            all_samples = (
                db.query(Sample)
                .filter(Sample.dataset_id.in_(ds_ids))
                .order_by(Sample.sample_name)
                .all()
            )

    if not datasets:
        return html.P("No pipeline runs yet.", className="text-muted")

    checked_set = set(checked_history) if checked_history else set()

    # Build dataset lookup
    ds_map = {d.id: d for d in datasets}

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }

    rows = []
    for s in all_samples:
        d = ds_map.get(s.dataset_id)
        if not d:
            continue

        raw = f"{s.read_count_raw:,}" if s.read_count_raw else "\u2014"
        filt = f"{s.read_count_filtered:,}" if s.read_count_filtered else "\u2014"
        nonchim = f"{s.read_count_nonchimeric:,}" if s.read_count_nonchimeric else "\u2014"

        badge = dbc.Badge(
            d.status.upper(),
            color=badge_color.get(d.status, "secondary"),
        )
        if d.status == "complete" and not d.taxonomy_path:
            status_badge = html.Span([badge, dbc.Badge("No taxonomy", color="warning", className="ms-1")])
        else:
            status_badge = badge

        is_complete = d.status == "complete"
        check_key = f"{d.id}:{s.sample_name}"
        is_checked = check_key in checked_set

        checkbox_td = html.Td(
            dbc.Checkbox(
                id={"type": "hist-check", "index": check_key},
                value=is_checked,
            ) if is_complete else "",
            style={"width": "30px"},
        )

        is_active = d.status in ("pending", "processing")
        action_btns = []
        if is_complete and s.asv_table_path:
            action_btns.append(dbc.Button(
                "\u2b07",
                id={"type": "btn-sample-dl", "index": s.id},
                color="info",
                size="sm",
                outline=True,
                title="Download ASV table",
                style={"padding": "0 6px", "lineHeight": "1.4"},
                className="me-1",
            ))
        action_btns.append(dbc.Button(
            "x",
            id={"type": "btn-sample-del", "index": s.id},
            color="danger",
            size="sm",
            outline=True,
            disabled=is_active,
            title="Cannot delete while pipeline is running" if is_active else "Delete sample",
            style={"padding": "0 6px", "lineHeight": "1.4"},
        ))

        rows.append(
            html.Tr(
                [
                    checkbox_td,
                    html.Td(s.sample_name),
                    html.Td(status_badge),
                    html.Td(raw),
                    html.Td(filt),
                    html.Td(nonchim),
                    html.Td(f"{s.asv_count:,}" if s.asv_count else "\u2014"),
                    html.Td(d.variable_region or "\u2014", className="small"),
                    html.Td(
                        d.created_at.strftime("%Y-%m-%d") if d.created_at else "",
                        className="text-muted small",
                    ),
                    html.Td(html.Div(action_btns, className="d-flex")),
                ]
            )
        )

    # Also show datasets with no samples yet (pending/processing/failed)
    ds_with_samples = {s.dataset_id for s in all_samples}
    for d in datasets:
        if d.id in ds_with_samples:
            continue

        ds_dir = DATASET_DIR / str(d.id)
        log_exists = (
            d.status in ("failed", "cancelled")
            and (ds_dir / "pipeline.log").exists()
        )
        if log_exists:
            status_el = dbc.Button(
                d.status.upper(),
                id={"type": "btn-history-log", "index": d.id},
                color=badge_color.get(d.status, "secondary"),
                size="sm",
                title="Download pipeline log",
            )
        else:
            status_el = dbc.Badge(
                d.status.upper(),
                color=badge_color.get(d.status, "secondary"),
            )

        is_active = d.status in ("pending", "processing")
        del_btn = dbc.Button(
            "x",
            id={"type": "btn-history-del", "index": d.id},
            color="danger",
            size="sm",
            outline=True,
            disabled=is_active,
            title="Cannot delete a running pipeline" if is_active else "Delete",
        )

        rows.append(
            html.Tr(
                [
                    html.Td("", style={"width": "30px"}),
                    html.Td("\u2014"),
                    html.Td(status_el),
                    html.Td("\u2014"),
                    html.Td("\u2014"),
                    html.Td("\u2014"),
                    html.Td("\u2014"),
                    html.Td(d.variable_region or "\u2014", className="small"),
                    html.Td(
                        d.created_at.strftime("%Y-%m-%d") if d.created_at else "",
                        className="text-muted small",
                    ),
                    html.Td(del_btn),
                ],
                className="text-muted",
            )
        )

    if not rows:
        return html.P("No pipeline runs yet.", className="text-muted")

    return dbc.Table(
        [
            html.Thead(
                html.Tr(
                    [
                        html.Th(
                            dbc.Checkbox(id="hist-check-all", value=False),
                            style={"width": "30px"},
                        ),
                        html.Th("Sample"),
                        html.Th("Status"),
                        html.Th("Raw Reads"),
                        html.Th("Filtered"),
                        html.Th("Non-chimeric"),
                        html.Th("ASVs"),
                        html.Th("Region"),
                        html.Th("Date"),
                        html.Th(""),
                    ]
                )
            ),
            html.Tbody(rows),
        ],
        bordered=True,
        hover=True,
        size="sm",
    )


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — Upload (Section A)
# ══════════════════════════════════════════════════════════════════════════════


@du.callback(
    output=Output("fm-staging-dir", "data"),
    id="fm-du-upload",
)
def on_file_uploaded(file_paths):
    if not file_paths:
        return no_update
    latest = Path(file_paths[0])
    return str(latest.parent)


@dash_app.callback(
    Output("fm-upload-file-status", "children"),
    Output("fm-btn-register", "disabled"),
    Output("fm-staging-dir", "data", allow_duplicate=True),
    Input("fm-staging-dir", "data"),
    Input("url", "pathname"),
    prevent_initial_call=True,
)
def update_upload_status(staging_dir, _pathname):
    # Check stored path first
    if staging_dir:
        staging = Path(staging_dir)
        if staging.exists():
            all_files = [
                f.name for f in sorted(staging.iterdir())
                if f.is_file() and (f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"))
            ]
            if all_files:
                n = len(all_files)
                return (
                    html.Small(f"{n} file(s) transferred: {', '.join(all_files)}", className="text-success"),
                    False,
                    no_update,
                )

    # Scan for any orphaned staging subdirectory (e.g. after disconnect)
    for child in UPLOAD_DIR.iterdir():
        if child.is_dir():
            fastqs = [
                f.name for f in sorted(child.iterdir())
                if f.is_file() and (f.name.endswith(".fastq.gz") or f.name.endswith(".fq.gz"))
            ]
            if fastqs:
                n = len(fastqs)
                return (
                    html.Small(f"{n} file(s) ready to register: {', '.join(fastqs)}", className="text-info"),
                    False,
                    str(child),
                )

    return "", True, no_update


@dash_app.callback(
    Output("fm-upload-status", "children"),
    Output("fm-upload-trigger", "data"),
    Output("fm-staging-dir", "data", allow_duplicate=True),
    Output("fm-upload-file-status", "children", allow_duplicate=True),
    Output("fm-btn-register", "disabled", allow_duplicate=True),
    Input("fm-btn-register", "n_clicks"),
    State("fm-staging-dir", "data"),
    State("fm-upload-trigger", "data"),
    prevent_initial_call=True,
)
def on_register_upload(n_clicks, staging_dir_path, trigger):
    # Find staging directory: use stored path, or scan UPLOAD_DIR for any UUID subdirectory
    staging_dir = None
    if staging_dir_path:
        p = Path(staging_dir_path)
        if p.exists():
            staging_dir = p

    if staging_dir is None:
        # Scan for any staging subdirectory left by dash_uploader
        for child in UPLOAD_DIR.iterdir():
            if child.is_dir() and any(
                f.name.endswith((".fastq.gz", ".fq.gz")) for f in child.iterdir() if f.is_file()
            ):
                staging_dir = child
                break

    if staging_dir is None:
        return dbc.Alert("No upload in progress.", color="warning"), no_update, no_update, no_update, no_update

    staged_files = []
    for fpath in sorted(staging_dir.iterdir()):
        if fpath.is_dir():
            continue
        if fpath.name.endswith(".fastq.gz") or fpath.name.endswith(".fq.gz"):
            staged_files.append(fpath)
        else:
            fpath.unlink(missing_ok=True)

    if not staged_files:
        shutil.rmtree(staging_dir, ignore_errors=True)
        return dbc.Alert("No FASTQ files found in upload.", color="warning"), no_update, None, no_update, no_update

    saved_filenames = [f.name for f in staged_files]

    db = SessionLocal()
    try:
        existing = (
            db.query(FastqFile.filename)
            .filter(FastqFile.filename.in_(saved_filenames))
            .all()
        )
        dupes = [row.filename for row in existing]
    finally:
        db.close()

    if dupes:
        shutil.rmtree(staging_dir, ignore_errors=True)
        dupe_list = html.Ul([html.Li(f, className="text-warning") for f in sorted(dupes)])
        return (
            dbc.Alert(
                [
                    html.P("These filenames are already registered:", className="mb-1 fw-bold"),
                    dupe_list,
                    html.P("Please rename the files and upload again.", className="mb-0 mt-2"),
                ],
                color="danger",
            ),
            no_update,
            None,
            "",
            True,
        )

    total_size = 0.0
    for fpath in staged_files:
        dest = UPLOAD_DIR / fpath.name
        fpath.rename(dest)
        total_size += dest.stat().st_size / (1024 * 1024)
    shutil.rmtree(staging_dir, ignore_errors=True)

    detection = detect_sequencing_type(saved_filenames)

    variable_region = None
    first_sample = next(iter(detection["samples"].values()), {})
    r1_name = first_sample.get("R1")
    if r1_name:
        try:
            region_result = detect_variable_region(UPLOAD_DIR / r1_name)
            variable_region = region_result["region"]
        except Exception:
            pass

    primers_detected = None
    if variable_region and r1_name:
        try:
            fwd_primer = REGION_PRIMERS[variable_region]["forward"]
            seqs = _read_fastq_sequences(UPLOAD_DIR / r1_name, n_reads=100)
            if seqs:
                matches = sum(1 for s in seqs if _primer_matches(s, fwd_primer))
                primers_detected = (matches / len(seqs)) >= 0.30
        except Exception:
            pass

    db = SessionLocal()
    try:
        upload = Upload(
            upload_dir=str(UPLOAD_DIR),
            sequencing_type=detection["type"],
            variable_region=variable_region,
            primers_detected=primers_detected,
            study=None,
            total_files=len(saved_filenames),
            total_size_mb=round(total_size, 2),
            status="uploaded",
        )
        db.add(upload)
        db.flush()

        for filename in saved_filenames:
            sample_name = extract_sample_name(filename)
            sample_info = detection["samples"].get(sample_name, {})
            if sample_info.get("R1") == filename:
                read_direction = "R1"
            elif sample_info.get("R2") == filename:
                read_direction = "R2"
            else:
                read_direction = "single"

            fpath = UPLOAD_DIR / filename
            avg_len = _avg_read_length(fpath)
            n_reads = _count_reads(fpath)
            db.add(
                FastqFile(
                    upload_id=upload.id,
                    sample_name=sample_name,
                    filename=filename,
                    file_path=to_relative(fpath),
                    read_direction=read_direction,
                    file_size_mb=round(fpath.stat().st_size / (1024 * 1024), 2),
                    read_count=n_reads,
                    avg_read_length=avg_len,
                )
            )

        db.commit()
    except Exception as e:
        db.rollback()
        return dbc.Alert(f"Database error: {e}", color="danger"), no_update, None, no_update, no_update
    finally:
        db.close()

    region_str = f", region: {variable_region}" if variable_region else ""
    if primers_detected is True:
        primer_str = ", primers: detected"
    elif primers_detected is False:
        primer_str = ", primers: not detected"
    else:
        primer_str = ""
    msg = (
        f"Registered {len(saved_filenames)} files "
        f"({_human_size(total_size)}), {detection['type']}{region_str}{primer_str}"
    )
    warnings = detection.get("errors", [])
    alert_children = [html.P(msg, className="mb-0")]
    if warnings:
        alert_children.append(
            html.Small(
                " | ".join(warnings), className="text-warning d-block mt-1"
            )
        )

    return (
        dbc.Alert(alert_children, color="success"),
        (trigger or 0) + 1,
        None,
        "",
        True,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — File Table Sort/Filter/Check (Section B, file management)
# ══════════════════════════════════════════════════════════════════════════════


_SORT_BUTTONS = [
    ("fm-sort-sample", "sample_name"),
    ("fm-sort-type", "direction"),
    ("fm-sort-reads", "total_reads"),
    ("fm-sort-region", "region"),
    ("fm-sort-readlen", "avg_read_length"),
    ("fm-sort-date", "created_at"),
]


@dash_app.callback(
    Output("fm-sort-field", "data"),
    Output("fm-sort-asc", "data"),
    *[Output(btn_id, "outline") for btn_id, _ in _SORT_BUTTONS],
    *[Input(btn_id, "n_clicks") for btn_id, _ in _SORT_BUTTONS],
    State("fm-sort-field", "data"),
    State("fm-sort-asc", "data"),
    prevent_initial_call=True,
)
def on_sort_click(*args):
    current_field = args[-2]
    current_asc = args[-1]
    triggered = ctx.triggered_id
    sort_map = {btn_id: field for btn_id, field in _SORT_BUTTONS}
    new_field = sort_map.get(triggered, "sample_name")

    if new_field == current_field:
        new_asc = not current_asc
    else:
        new_asc = True

    outlines = [new_field != field for _, field in _SORT_BUTTONS]
    return new_field, new_asc, *outlines


@dash_app.callback(
    Output("fm-files-table", "children"),
    Output("file-ids-map-store", "data"),
    Input("fm-upload-trigger", "data"),
    Input("fm-sort-field", "data"),
    Input("fm-sort-asc", "data"),
    Input("fm-f-sample", "value"),
    Input("fm-f-type", "value"),
    Input("fm-f-region", "value"),
    State("fm-checked-samples", "data"),
)
def refresh_files_table(trigger, sort_by, sort_asc, f_sample, f_type, f_region, checked):
    """Rebuild the files table whenever data, sort, or filter changes."""
    filters = {
        "sample": (f_sample or "").strip().lower(),
        "type": (f_type or "").strip().lower(),
        "region": (f_region or "").strip().lower(),
    }
    checked_set = set(checked) if checked else set()
    result = _build_files_table(sort_by or "sample_name", filters, sort_asc if sort_asc is not None else True, checked_samples=checked_set)
    if isinstance(result, tuple):
        return result
    # String or component without file_ids_map
    return result, {}


@dash_app.callback(
    Output("fm-checked-samples", "data"),
    Output("fm-btn-delete-selected", "disabled"),
    Output("fm-btn-delete-selected", "children"),
    Input({"type": "fm-check", "index": ALL}, "value"),
    State({"type": "fm-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_row_check(values, ids):
    checked = [
        id_dict["index"] for id_dict, val in zip(ids, values) if val
    ]
    n = len(checked)
    if n:
        return checked, False, f"Delete Selected ({n})"
    return checked, True, "Delete Selected"


@dash_app.callback(
    Output({"type": "fm-check", "index": ALL}, "value"),
    Input("fm-check-all", "value"),
    State({"type": "fm-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_check_all(select_all, ids):
    return [bool(select_all)] * len(ids)


@dash_app.callback(
    Output("fm-confirm-delete", "displayed"),
    Output("fm-confirm-delete", "message"),
    Input("fm-btn-delete-selected", "n_clicks"),
    State("fm-checked-samples", "data"),
    prevent_initial_call=True,
)
def on_delete_click(n_clicks, checked):
    if not n_clicks or not checked:
        return False, ""
    n = len(checked)
    names = ", ".join(sorted(checked)[:10])
    suffix = f"... and {n - 10} more" if n > 10 else ""
    return True, (
        f"Permanently delete {n} sample(s) and their FASTQ files?\n\n"
        f"{names}{suffix}\n\n"
        "This cannot be undone."
    )


@dash_app.callback(
    Output("fm-upload-trigger", "data", allow_duplicate=True),
    Output("fm-checked-samples", "data", allow_duplicate=True),
    Input("fm-confirm-delete", "submit_n_clicks"),
    State("fm-checked-samples", "data"),
    State("fm-upload-trigger", "data"),
    prevent_initial_call=True,
)
def on_delete_confirmed(submit_n_clicks, checked, trigger):
    if not submit_n_clicks or not checked:
        return no_update, no_update

    checked_set = set(checked)
    db = SessionLocal()
    try:
        fastq_records = (
            db.query(FastqFile)
            .filter(FastqFile.sample_name.in_(checked_set))
            .all()
        )
        affected_upload_ids = set()
        for ff in fastq_records:
            affected_upload_ids.add(ff.upload_id)
            p = to_absolute(ff.file_path)
            if p:
                p.unlink(missing_ok=True)
            db.delete(ff)

        db.flush()

        for uid in affected_upload_ids:
            remaining = (
                db.query(FastqFile)
                .filter(FastqFile.upload_id == uid)
                .count()
            )
            if remaining == 0:
                upload = db.query(Upload).filter(Upload.id == uid).first()
                if upload:
                    db.delete(upload)

        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()

    return (trigger or 0) + 1, []


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — DADA2 Pipeline (Section B, pipeline part)
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("sample-selection-summary", "children"),
    Output("btn-start-pipeline", "disabled"),
    Output("error-model-container", "style"),
    Output("input-error-model", "value"),
    Input("fm-checked-samples", "data"),
    State("file-ids-map-store", "data"),
)
def update_selection_summary(checked_samples, file_ids_map):
    hidden = {"display": "none"}
    if not checked_samples or not file_ids_map:
        return "", True, hidden, "default"

    from app.db.database import get_session
    from app.db.models import FastqFile as FF, Upload as UL

    with get_session() as db:
        results = (
            db.query(FF, UL.variable_region)
            .join(UL, FF.upload_id == UL.id)
            .filter(FF.sample_name.in_(checked_samples))
            .all()
        )

    sample_info = defaultdict(lambda: {"directions": set(), "region": None})
    for f, region in results:
        s = sample_info[f.sample_name]
        s["directions"].add(f.read_direction)
        if region:
            s["region"] = region

    n = len(sample_info)
    if n == 0:
        return "", True, hidden, "default"

    types = set()
    regions = set()
    for info in sample_info.values():
        dirs = info["directions"]
        if "R1" in dirs and "R2" in dirs:
            types.add("PE")
        elif dirs == {"single"}:
            types.add("SE")
        else:
            types.add(", ".join(sorted(dirs)))
        if info["region"]:
            regions.add(info["region"])

    badges = []
    for t in sorted(types):
        badges.append(dbc.Badge(
            t, color="success" if t == "PE" else "info", className="me-1",
        ))
    for r in sorted(regions):
        badges.append(dbc.Badge(r, color="primary", className="me-1"))

    warnings = []
    if len(types) > 1:
        warnings.append(dbc.Badge("Mixed SE/PE", color="warning", className="me-1"))
    blocked = False
    if len(regions) > 1:
        warnings.append(dbc.Badge("Mixed regions — cannot run DADA2", color="danger", className="me-1"))
        blocked = True

    summary = html.Div([
        html.Span(f"{n} sample{'s' if n != 1 else ''} selected  ", className="me-2"),
        *badges,
        *warnings,
    ])

    # Show error model dropdown when V1-V9 (long-read) is detected
    has_long_read = "V1-V9" in regions
    error_model_style = {"display": "block"} if has_long_read else hidden
    error_model_value = "PacBio" if has_long_read else "default"

    return summary, blocked, error_model_style, error_model_value


@dash_app.callback(
    Output("pipeline-active-dataset", "data"),
    Output("pipeline-progress-section", "style"),
    Output("pipeline-poll", "disabled"),
    Output("pipeline-launch-error", "children"),
    Output("pipeline-detection-summary", "children"),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("btn-start-pipeline", "n_clicks"),
    State("fm-checked-samples", "data"),
    State("file-ids-map-store", "data"),
    State("input-threads", "value"),
    State("input-error-model", "value"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def on_start_pipeline(n_clicks, checked_samples, file_ids_map, threads_val, error_model_val, checked_history):
    if not n_clicks or not checked_samples or not file_ids_map:
        return no_update, no_update, no_update, no_update, no_update, no_update

    no_change = (no_update, no_update, no_update)

    threads_override = None
    use_default_threads = True
    if threads_val and str(threads_val).strip().lower() != "default":
        use_default_threads = False
        try:
            threads_override = max(1, int(threads_val))
        except ValueError:
            return (
                *no_change,
                dbc.Alert(
                    f"Invalid CPU threads value: '{threads_val}'. Enter a number or 'default'.",
                    color="danger",
                ),
                no_update, no_update,
            )

    from datetime import datetime as _dt
    dataset_name = f"Run_{_dt.utcnow().strftime('%Y%m%d_%H%M%S')}"

    all_file_ids = []
    for sample_name in checked_samples:
        ids = file_ids_map.get(sample_name, [])
        all_file_ids.extend(ids)

    if not all_file_ids:
        return (
            *no_change,
            dbc.Alert("No files found for selected samples.", color="danger"),
            no_update, no_update,
        )

    n_samples = len(checked_samples)

    from app.db.database import get_session
    from app.db.models import FastqFile as FF, Upload as UL

    with get_session() as db:
        sel_results = (
            db.query(FF, UL.variable_region)
            .join(UL, FF.upload_id == UL.id)
            .filter(FF.sample_name.in_(checked_samples))
            .all()
        )

    sel_info = defaultdict(lambda: {"directions": set(), "region": None})
    for f, region in sel_results:
        s = sel_info[f.sample_name]
        s["directions"].add(f.read_direction)
        if region:
            s["region"] = region

    types = set()
    regions = set()
    for info in sel_info.values():
        dirs = info["directions"]
        if "R1" in dirs and "R2" in dirs:
            types.add("PE")
        elif dirs == {"single"}:
            types.add("SE")
        else:
            types.add(", ".join(sorted(dirs)))
        if info["region"]:
            regions.add(info["region"])

    if len(types) == 1:
        type_str = next(iter(types))
        type_color = "success" if type_str == "PE" else "info"
    else:
        type_str = "MIXED"
        type_color = "warning"

    if len(regions) == 1:
        region_str = next(iter(regions))
        region_color = "primary"
    elif len(regions) > 1:
        region_str = "Mixed regions"
        region_color = "warning"
    else:
        region_str = "Region unknown"
        region_color = "secondary"

    is_v19 = region_str == "V1-V9"
    if is_v19:
        trunc_note = f"Long-read mode: no truncation, error model = {error_model_val or 'PacBio'}"
    else:
        trunc_note = "Truncation parameters will be auto-detected from quality profiles."

    card_items = [
        html.Div([
            dbc.Badge(type_str.upper(), color=type_color, className="me-2 fs-6"),
            dbc.Badge(region_str, color=region_color, className="me-2 fs-6"),
            dbc.Badge("Long-read", color="info", className="me-2 fs-6") if is_v19 else None,
            html.Span(f"{n_samples} sample(s)  |  {len(all_file_ids)} file(s)"),
        ]),
        html.Small(
            trunc_note,
            className="text-muted d-block mt-2",
        ),
    ]
    # Filter out None badges
    card_items[0].children = [c for c in card_items[0].children if c is not None]
    summary_card = dbc.Card(dbc.CardBody(card_items), className="mb-3")

    if use_default_threads:
        threads = min(n_samples * 2, max(1, MAX_CPUS - 1))
    else:
        threads = threads_override

    error_model = error_model_val or "default"

    from app.pipeline.runner import launch_dada2_pipeline

    try:
        dataset_id = launch_dada2_pipeline(
            file_ids=all_file_ids,
            name=dataset_name,
            threads=threads,
            error_model=error_model,
        )
    except Exception as e:
        return (
            *no_change,
            dbc.Alert(f"Pipeline launch failed: {e}", color="danger"),
            summary_card,
            no_update,
        )

    return (
        dataset_id,
        {"display": "block"},
        False,
        "",
        summary_card,
        _build_history_table(checked_history),
    )


@dash_app.callback(
    Output("pipeline-launch-error", "children", allow_duplicate=True),
    Output("pipeline-progress-section", "style", allow_duplicate=True),
    Output("pipeline-status-badge", "children", allow_duplicate=True),
    Output("pipeline-progress-bar", "value", allow_duplicate=True),
    Output("pipeline-progress-bar", "label", allow_duplicate=True),
    Output("pipeline-progress-bar", "animated", allow_duplicate=True),
    Output("pipeline-log-viewer", "children", allow_duplicate=True),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("btn-cancel-pipeline", "n_clicks"),
    State("pipeline-active-dataset", "data"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def on_cancel_pipeline(n_clicks, dataset_id, checked_history):
    if not n_clicks or not dataset_id:
        return (no_update,) * 9

    import time
    from app.pipeline.runner import cancel_pipeline, get_pipeline_status

    ok = cancel_pipeline(dataset_id)
    if not ok:
        return (
            dbc.Alert("No running pipeline found to cancel.", color="secondary"),
            {"display": "none"},
            "", 0, "", False, "",
            True,
            _build_history_table(checked_history),
        )

    time.sleep(0.5)

    status = get_pipeline_status(dataset_id)
    badge = dbc.Badge("CANCELLED", color="warning", className="fs-6 p-2")
    log_tail = status.get("log_tail", "Pipeline cancelled by user.")

    return (
        dbc.Alert("Pipeline cancelled.", color="warning"),
        {"display": "block"},
        badge,
        0,
        "Cancelled",
        False,
        log_tail,
        True,
        _build_history_table(checked_history),
    )


@dash_app.callback(
    Output("pipeline-status-badge", "children"),
    Output("pipeline-progress-bar", "value"),
    Output("pipeline-progress-bar", "label"),
    Output("pipeline-steps-display", "children"),
    Output("pipeline-log-viewer", "children"),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Output("pipeline-download-section", "style"),
    Output("pipeline-progress-bar-wrapper", "style"),
    Input("pipeline-poll", "n_intervals"),
    State("pipeline-active-dataset", "data"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def on_poll(n_intervals, dataset_id, checked_history):
    if not dataset_id:
        return (no_update,) * 9

    from app.pipeline.runner import get_pipeline_status

    status = get_pipeline_status(dataset_id)
    db_status = status["status"]
    pct = status["progress_pct"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }.get(db_status, "secondary")
    badge = dbc.Badge(
        db_status.upper(), color=badge_color, className="fs-6 p-2"
    )

    completed = set(status.get("steps_completed", []))
    current = status.get("current_step")
    step_items = []
    for step in STEPS:
        if step in completed and step != current:
            icon = "  "
            color = "text-success"
        elif step == current and db_status == "processing":
            icon = "  "
            color = "text-primary"
        elif step == current and db_status == "failed":
            icon = "  "
            color = "text-danger"
        else:
            icon = "  "
            color = "text-muted"
        step_items.append(
            html.Div(
                f"{icon} {STEP_LABELS.get(step, step)}",
                className=f"{color} mb-1",
            )
        )

    auto_info = None
    if status.get("auto_trunc_details"):
        auto_info = dbc.Alert(
            [
                html.Strong("Auto-detected: "),
                html.Span(status["auto_trunc_details"]),
            ],
            color="info",
            className="mt-2 py-2 small",
        )

    stop_poll = db_status in ("complete", "failed", "cancelled")
    history = _build_history_table(checked_history)

    steps_display = html.Div(step_items + ([auto_info] if auto_info else []))

    download_style = (
        {"display": "block"} if db_status == "complete" else {"display": "none"}
    )

    finished = db_status in ("complete", "failed", "cancelled")
    bar_style = {"display": "none"} if finished else {"display": "block"}

    return (
        badge,
        pct,
        f"{pct}%",
        steps_display,
        status.get("log_tail", ""),
        stop_poll,
        history,
        download_style,
        bar_style,
    )


@dash_app.callback(
    Output("pipeline-history-table", "children"),
    Input("pipeline-active-dataset", "id"),
    State("an-checked-history", "data"),
)
def load_history(_, checked_history):
    return _build_history_table(checked_history)


@dash_app.callback(
    Output("pipeline-active-dataset", "data", allow_duplicate=True),
    Output("pipeline-progress-section", "style", allow_duplicate=True),
    Output("pipeline-poll", "disabled", allow_duplicate=True),
    Output("pipeline-status-badge", "children", allow_duplicate=True),
    Output("pipeline-progress-bar", "value", allow_duplicate=True),
    Output("pipeline-progress-bar", "label", allow_duplicate=True),
    Output("pipeline-steps-display", "children", allow_duplicate=True),
    Output("pipeline-log-viewer", "children", allow_duplicate=True),
    Output("pipeline-download-section", "style", allow_duplicate=True),
    Output("pipeline-progress-bar-wrapper", "style", allow_duplicate=True),
    Input("pipeline-active-dataset", "id"),
    State("pipeline-active-dataset", "data"),
    prevent_initial_call="initial_duplicate",
)
def restore_progress_on_load(_, dataset_id):
    from app.db.database import get_session
    from app.db.models import Dataset
    from app.pipeline.runner import get_pipeline_status

    no_all = (no_update,) * 10

    if dataset_id:
        with get_session() as db:
            ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if not ds or ds.status not in ("pending", "processing"):
                return (None, {"display": "none"}, True,
                        no_update, no_update, no_update, no_update,
                        no_update, no_update, no_update)
    else:
        with get_session() as db:
            running = (
                db.query(Dataset)
                .filter(Dataset.status.in_(["pending", "processing"]))
                .order_by(Dataset.created_at.desc())
                .first()
            )
            if running:
                dataset_id = running.id

    if not dataset_id:
        return no_all

    status = get_pipeline_status(dataset_id)
    db_status = status["status"]
    pct = status["progress_pct"]

    badge_color = {
        "pending": "secondary",
        "processing": "primary",
        "complete": "success",
        "failed": "danger",
        "cancelled": "warning",
    }.get(db_status, "secondary")
    badge = dbc.Badge(db_status.upper(), color=badge_color, className="fs-6 p-2")

    completed = set(status.get("steps_completed", []))
    current = status.get("current_step")
    step_items = []
    for step in STEPS:
        if step in completed and step != current:
            icon = "  "
            color = "text-success"
        elif step == current and db_status == "processing":
            icon = "  "
            color = "text-primary"
        elif step == current and db_status == "failed":
            icon = "  "
            color = "text-danger"
        else:
            icon = "  "
            color = "text-muted"
        step_items.append(
            html.Div(
                f"{icon} {STEP_LABELS.get(step, step)}",
                className=f"{color} mb-1",
            )
        )
    steps_display = html.Div(step_items)
    log_tail = status.get("log_tail", "")
    bar_style = {"display": "block"}

    keep_polling = db_status in ("pending", "processing")

    return (
        dataset_id,
        {"display": "block"},
        not keep_polling,
        badge,
        pct,
        f"{pct}%",
        steps_display,
        log_tail,
        {"display": "block"} if db_status == "complete" else {"display": "none"},
        {"display": "none"} if not keep_polling else bar_style,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — Download / History (Section D)
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("download-history-log", "data"),
    Input({"type": "btn-history-log", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_history_log_download(n_clicks_list):
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()

    log_path = DATASET_DIR / str(dataset_id) / "pipeline.log"
    if not log_path.exists():
        return no_update

    return dcc.send_file(
        str(log_path),
        filename=_download_filename(ds, "pipeline_log", "txt"),
    )


@dash_app.callback(
    Output("download-sample-asv", "data"),
    Input({"type": "btn-sample-dl", "index": ALL}, "n_clicks"),
    prevent_initial_call=True,
)
def on_sample_asv_download(n_clicks_list):
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    sample_id = triggered["index"]

    from app.db.database import get_session
    from app.db.models import Sample

    with get_session() as db:
        sample = db.query(Sample).filter(Sample.id == sample_id).first()

    if not sample or not sample.asv_table_path:
        return no_update

    p = to_absolute(sample.asv_table_path)
    if not p or not p.exists():
        return no_update

    return dcc.send_file(str(p), filename=f"{sample.sample_name}_asv_table.csv.gz")


@dash_app.callback(
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input({"type": "btn-history-del", "index": ALL}, "n_clicks"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def on_history_delete(n_clicks_list, checked_history):
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    dataset_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset

    with get_session() as db:
        ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not ds:
            return _build_history_table(checked_history)

        if ds.status in ("pending", "processing"):
            return no_update

        db.delete(ds)
        db.commit()

    ds_path = DATASET_DIR / str(dataset_id)
    if ds_path.exists():
        shutil.rmtree(ds_path, ignore_errors=True)

    return _build_history_table(checked_history)


@dash_app.callback(
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input({"type": "btn-sample-del", "index": ALL}, "n_clicks"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def on_sample_delete(n_clicks_list, checked_history):
    if not any(n_clicks_list):
        return no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update
    sample_id = triggered["index"]

    from app.config import DATASET_DIR
    from app.db.database import get_session
    from app.db.models import Dataset, Sample

    with get_session() as db:
        sample = db.query(Sample).filter(Sample.id == sample_id).first()
        if not sample:
            return _build_history_table(checked_history)

        dataset = db.query(Dataset).filter(Dataset.id == sample.dataset_id).first()
        if dataset and dataset.status in ("pending", "processing"):
            return no_update

        # Delete the sample's csv.gz file
        if sample.asv_table_path:
            p = to_absolute(sample.asv_table_path)
            if p and p.exists():
                p.unlink()

        dataset_id = sample.dataset_id
        db.delete(sample)
        db.flush()

        # If no samples remain, delete the dataset and its directory
        remaining = db.query(Sample).filter(Sample.dataset_id == dataset_id).count()
        if remaining == 0 and dataset:
            db.delete(dataset)
            db.flush()

    if remaining == 0:
        ds_path = DATASET_DIR / str(dataset_id)
        if ds_path.exists():
            shutil.rmtree(ds_path, ignore_errors=True)

    return _build_history_table(checked_history)


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — History Checkboxes + Analysis Selection
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("an-checked-history", "data"),
    Input({"type": "hist-check", "index": ALL}, "value"),
    State({"type": "hist-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_hist_row_check(values, ids):
    """Sync history table checkboxes → store."""
    checked = [
        id_dict["index"] for id_dict, val in zip(ids, values) if val
    ]
    return checked


@dash_app.callback(
    Output({"type": "hist-check", "index": ALL}, "value"),
    Input("hist-check-all", "value"),
    State({"type": "hist-check", "index": ALL}, "id"),
    prevent_initial_call=True,
)
def on_hist_check_all(select_all, ids):
    """Toggle all history checkboxes."""
    return [bool(select_all)] * len(ids)


@dash_app.callback(
    Output("an-selection-summary", "children"),
    Output("run-analysis-btn", "disabled"),
    Input("an-checked-history", "data"),
)
def update_analysis_selection(checked_items):
    """Update analysis selection summary from history checkboxes."""
    if not checked_items:
        return html.P("Select samples from Pipeline History above.", className="text-muted"), True

    by_dataset = defaultdict(list)
    for item in checked_items:
        parts = item.split(":", 1)
        if len(parts) == 2:
            try:
                by_dataset[int(parts[0])].append(parts[1])
            except ValueError:
                pass

    if not by_dataset:
        return html.P("Select samples from Pipeline History above.", className="text-muted"), True

    n = sum(len(v) for v in by_dataset.values())
    return (
        html.Div([
            dbc.Badge(f"{n} sample{'s' if n != 1 else ''} selected", color="primary", className="me-2 fs-6"),
        ]),
        False,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — Analysis: Source Tracking + Pathogen (Section C)
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("st-run-id", "data"),
    Output("tax-run-id", "data"),
    Output("st-poll", "disabled", allow_duplicate=True),
    Output("tax-poll", "disabled", allow_duplicate=True),
    Output("st-progress-area", "children", allow_duplicate=True),
    Output("tax-progress-area", "children", allow_duplicate=True),
    Output("biom-taxonomy-dataset-id", "data", allow_duplicate=True),
    Input("run-analysis-btn", "n_clicks"),
    State("an-checked-history", "data"),
    State("st-mode", "value"),
    State("st-src-depth", "value"),
    State("st-snk-depth", "value"),
    State("st-restarts", "value"),
    State("st-burnin", "value"),
    State("st-draws", "value"),
    State({"type": "group-check", "index": ALL}, "value"),
    State({"type": "group-check", "index": ALL}, "id"),
    State("st-threads", "value"),
    State("tax-identity", "value"),
    State("tax-threads", "value"),
    prevent_initial_call=True,
)
def launch_analysis(n_clicks, checked_history,
                    mode, src_depth, snk_depth, restarts, burnin, draws,
                    group_values, group_ids, st_threads, identity, tax_threads):
    if not n_clicks:
        return (no_update,) * 7

    from app.pipeline.runner import launch_sourcetracker_run, launch_pathogen_run
    from app.db.database import get_session
    from app.db.models import Dataset, Sample

    # Parse checked history items → dataset_id + sample names
    by_dataset = defaultdict(list)
    for item in (checked_history or []):
        parts = item.split(":", 1)
        if len(parts) == 2:
            try:
                by_dataset[int(parts[0])].append(parts[1])
            except ValueError:
                pass

    if not by_dataset:
        err = dbc.Alert("No samples selected. Check samples in Pipeline History above.", color="danger")
        return no_update, no_update, no_update, no_update, err, err, no_update

    # Query per-sample ASV table paths across all selected datasets
    with get_session() as db:
        sample_table_paths = []
        has_taxonomy = True
        for dataset_id, sample_names in by_dataset.items():
            ds = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if ds and not ds.taxonomy_path:
                has_taxonomy = False
            samples = (
                db.query(Sample)
                .filter(Sample.dataset_id == dataset_id, Sample.sample_name.in_(sample_names))
                .all()
            )
            sample_table_paths.extend(str(to_absolute(s.asv_table_path)) for s in samples if s.asv_table_path)

    if not sample_table_paths:
        err = dbc.Alert("No per-sample ASV tables found for selected samples.", color="danger")
        return no_update, no_update, no_update, no_update, err, err, no_update

    # Read checked source groups from pattern-matching checkboxes
    selected_groups = []
    for id_dict, vals in zip(group_ids, group_values):
        if vals:  # vals is a list like ["pig"] when checked, [] when unchecked
            selected_groups.append(id_dict["index"])

    # Launch SourceTracker
    st_run_id = launch_sourcetracker_run(
        sample_table_paths=sample_table_paths,
        selected_groups=selected_groups,
        feature_mode=mode or "asv",
        src_depth=src_depth or 1000,
        snk_depth=snk_depth or 1000,
        restarts=restarts or 10,
        burnin=burnin or 100,
        draws=draws or 1,
        dataset_id=dataset_id,
        threads=st_threads or 4,
    )

    st_progress = dbc.Card(dbc.CardBody([
        html.H6("Source Tracking Running..."),
        dbc.Progress(value=0, striped=True, animated=True),
    ]), className="mb-3")

    # Launch Pathogen only if taxonomy is available
    if has_taxonomy:
        tax_run_id = launch_pathogen_run(
            sample_table_paths=sample_table_paths,
            silva_identity=(identity or 97) / 100.0,
            dataset_id=dataset_id,
            threads=tax_threads or 4,
        )
        tax_progress = dbc.Card(dbc.CardBody([
            html.H6("Pathogen Detection Running..."),
            dbc.Progress(value=0, striped=True, animated=True),
        ]), className="mb-3")
    else:
        tax_run_id = None
        tax_progress = dbc.Alert(
            [
                html.Strong("Pathogen detection requires taxonomy. "),
                "Selected dataset has no taxonomy data. ",
                "Run taxonomy assignment first, then re-run analysis.",
                html.Br(),
                dbc.Button(
                    "Run Taxonomy",
                    id="btn-run-taxonomy-biom",
                    color="warning",
                    size="sm",
                    className="mt-2",
                ),
            ],
            color="warning",
            className="mt-2",
        )

    biom_tax_ds = dataset_id if not has_taxonomy else no_update
    return st_run_id, tax_run_id, False, not has_taxonomy, st_progress, tax_progress, biom_tax_ds


@dash_app.callback(
    Output("st-progress-area", "children"),
    Output("st-results-area", "children"),
    Output("st-poll", "disabled"),
    Input("st-poll", "n_intervals"),
    State("st-run-id", "data"),
    prevent_initial_call=True,
)
def poll_st_status(n_intervals, run_id):
    if not run_id:
        return no_update, no_update, True

    from app.pipeline.runner import get_sourcetracker_status

    status = get_sourcetracker_status(run_id)
    pct = status.get("progress_pct", 0)
    db_status = status.get("status", "unknown")

    if db_status == "complete":
        from app.pipeline.sourcetracker import GROUP_COLORS
        import pandas as pd
        from app.config import DATASET_DIR

        run_dir = DATASET_DIR / f"st_{run_id}"
        mp_path = run_dir / "proportions.csv"

        low_samples = status.get("low_alignment_samples", [])
        sample_pcts = status.get("sample_alignment_pct", {})
        warn_banner = None
        if low_samples:
            detail_lines = [
                f"{s}: {sample_pcts.get(s, 0)}% reads retained"
                for s in low_samples
            ]
            warn_banner = dbc.Alert([
                html.Strong(
                    f"\u26a0 {len(low_samples)} sample(s) excluded "
                    f"due to low alignment (<10% reads matched):"
                ),
                html.Ul([html.Li(line) for line in detail_lines], className="mb-0 mt-2"),
            ], color="warning", className="mb-3")

        results_area = []
        if mp_path.exists():
            mp = pd.read_csv(mp_path, index_col=0)

            fig = go.Figure()
            for col in mp.columns:
                color = GROUP_COLORS.get(col, "#888888")
                fig.add_trace(go.Bar(
                    name=col,
                    x=mp.index.tolist(),
                    y=mp[col].values,
                    marker_color=color,
                    text=[f"{v:.0%}" if v > 0.05 else "" for v in mp[col].values],
                    textposition="inside",
                    textfont_color="white",
                    hovertemplate=f"{col} %{{y:.1%}}<extra></extra>",
                ))

            fig.update_layout(
                barmode="stack",
                title="Source Contribution (SourceTracker2 Gibbs)",
                yaxis_title="Proportion",
                yaxis_tickformat=".0%",
                template="plotly_white",
                legend_title="Source Group",
                height=400,
            )

            results_area = dbc.Card(dbc.CardBody([
                html.H6("Source Tracking Results", className="text-success"),
                dcc.Graph(figure=fig),
                html.H6("Mixing Proportions", className="mt-3"),
                dash_table.DataTable(
                    data=mp.round(4).reset_index().to_dict("records"),
                    columns=[{"name": c, "id": c} for c in ["index"] + list(mp.columns)],
                    style_header={"backgroundColor": "#2c3e50", "color": "white"},
                    style_cell={"backgroundColor": "white", "color": "#333", "textAlign": "center"},
                ),
            ]), className="mb-3")
        else:
            results_area = dbc.Alert("Results not found", color="warning")

        progress = dbc.Card(dbc.CardBody([
            html.H6("Source Tracking Complete"),
            dbc.Progress(value=100, color="success"),
        ]), className="mb-3")

        combined = [warn_banner, results_area] if warn_banner else [results_area]
        return progress, html.Div(combined), True

    if db_status == "failed":
        log_tail = status.get("log_tail", "")
        # Extract last meaningful error line
        err_lines = [l for l in log_tail.splitlines() if "ERROR" in l or "failed" in l.lower()]
        err_detail = err_lines[-1] if err_lines else "Check pipeline logs for details."
        return (
            dbc.Alert([
                html.Strong("Source Tracking failed."),
                html.Pre(err_detail, className="mt-2 mb-0 small",
                         style={"whiteSpace": "pre-wrap"}),
            ], color="danger"),
            no_update, True,
        )

    progress = dbc.Card(dbc.CardBody([
        html.H6("Source Tracking Running..."),
        dbc.Progress(value=pct, striped=True, animated=True, label=f"{pct}%"),
        html.P(f"Step: {status.get('current_step', '')}", className="mt-2 text-muted small"),
    ]), className="mb-3")

    return progress, no_update, False


@dash_app.callback(
    Output("tax-progress-area", "children"),
    Output("pathogen-results-area", "children"),
    Output("tax-poll", "disabled"),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("tax-poll", "n_intervals"),
    State("tax-run-id", "data"),
    State("biom-taxonomy-dataset-id", "data"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def poll_tax_status(n_intervals, run_id, biom_tax_ds_id, checked_history):
    # If there's a pathogen run, poll it
    if run_id:
        return _poll_pathogen(run_id) + (no_update,)

    # Otherwise, if there's a taxonomy-for-BIOM run, poll that
    if biom_tax_ds_id:
        return _poll_biom_taxonomy(biom_tax_ds_id, checked_history)

    return no_update, no_update, True, no_update


def _poll_pathogen(run_id):
    """Poll pathogen detection status. Returns (progress, results, disable_poll)."""
    from app.pipeline.runner import get_pathogen_status

    status = get_pathogen_status(run_id)
    pct = status.get("progress_pct", 0)
    db_status = status.get("status", "unknown")

    if db_status == "complete":
        import pandas as pd
        from app.config import DATASET_DIR
        from app.pipeline.pathogen import PATHOGENS, _PRIORITY_ORDER, genus_colors

        run_dir = DATASET_DIR / f"pathogen_{run_id}"
        ra_path = run_dir / "pathogen_ra.csv"
        summary_path = run_dir / "pathogen_summary.csv"

        results = []

        if ra_path.exists():
            ra_df = pd.read_csv(ra_path, index_col=0)

            genera = sorted(
                ra_df.columns,
                key=lambda g: (_PRIORITY_ORDER.index(PATHOGENS.get(g, "other")), g),
            )
            colors = genus_colors(genera)

            fig = go.Figure()
            for genus in genera:
                fig.add_trace(go.Bar(
                    name=genus,
                    x=ra_df.index.tolist(),
                    y=ra_df[genus].values,
                    marker_color=colors.get(genus, "#888"),
                    hovertemplate="%{y:.3f}%<extra>%{fullData.name}</extra>",
                ))

            fig.update_layout(
                barmode="stack",
                title="Pathogenic Bacteria  -  Relative Abundance (%)",
                yaxis_title="Relative Abundance (%)",
                template="plotly_white",
                legend_title="Genus",
                height=400,
            )

            results.append(dcc.Graph(figure=fig))

        if summary_path.exists():
            summary_df = pd.read_csv(summary_path)
            results.append(html.H6("Detection Summary", className="mt-3"))
            results.append(
                dash_table.DataTable(
                    data=summary_df.to_dict("records"),
                    columns=[{"name": c, "id": c} for c in summary_df.columns],
                    style_header={"backgroundColor": "#2c3e50", "color": "white"},
                    style_cell={"backgroundColor": "white", "color": "#333", "textAlign": "center"},
                )
            )

        if not results:
            results = [dbc.Alert("No pathogenic bacteria detected in this dataset.", color="info")]

        progress = dbc.Card(dbc.CardBody([
            html.H6("Pathogen Detection Complete"),
            dbc.Progress(value=100, color="success"),
        ]), className="mb-3")

        return (
            progress,
            dbc.Card(dbc.CardBody([html.H6("Pathogen Detection Results")] + results), className="mb-3"),
            True,
        )

    if db_status == "failed":
        log_tail = status.get("log_tail", "")
        err_lines = [l for l in log_tail.splitlines() if "ERROR" in l or "failed" in l.lower()]
        err_detail = err_lines[-1] if err_lines else "Check pipeline logs for details."
        return (
            dbc.Alert([
                html.Strong("Pathogen detection failed."),
                html.Pre(err_detail, className="mt-2 mb-0 small",
                         style={"whiteSpace": "pre-wrap"}),
            ], color="danger"),
            no_update, True,
        )

    progress = dbc.Card(dbc.CardBody([
        html.H6("Pathogen Detection Running..."),
        dbc.Progress(value=pct, striped=True, animated=True, label=f"{pct}%"),
        html.P(f"Step: {status.get('current_step', '')}", className="mt-2 text-muted small"),
    ]), className="mb-3")

    return progress, no_update, False


def _poll_biom_taxonomy(dataset_id, checked_history):
    """Poll taxonomy-for-BIOM status. Returns (progress, results, disable_poll, history)."""
    from app.pipeline.runner import get_taxonomy_status as _get_tax_status

    status = _get_tax_status(dataset_id)
    tax_status = status.get("status", "idle")
    pct = status.get("progress_pct", 0)

    if tax_status == "complete":
        progress = dbc.Card(dbc.CardBody([
            html.H6("Taxonomy Assignment Complete", className="text-success"),
            dbc.Progress(value=100, color="success"),
            html.P("You can now re-run Analysis to include Pathogen detection.", className="mt-2 small"),
        ]), className="mb-3")
        return progress, no_update, True, _build_history_table(checked_history)

    if tax_status == "failed":
        return (
            dbc.Alert("Taxonomy assignment failed. Check logs.", color="danger"),
            no_update, True, no_update,
        )

    if tax_status == "idle":
        return no_update, no_update, True, no_update

    progress = dbc.Card(dbc.CardBody([
        html.H6("Taxonomy Assignment Running..."),
        dbc.Progress(value=pct, striped=True, animated=True, label=f"{pct}%"),
        html.P(f"Step: {status.get('current_step', '')}", className="mt-2 text-muted small"),
    ]), className="mb-3")

    return progress, no_update, False, no_update


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — BIOM Import
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("biom-import-status", "children"),
    Output("pipeline-history-table", "children", allow_duplicate=True),
    Input("biom-upload", "contents"),
    State("biom-upload", "filename"),
    State("biom-name", "value"),
    State("an-checked-history", "data"),
    prevent_initial_call=True,
)
def import_biom_file(contents, filename, name, checked_history):
    """Handle BIOM file upload and import."""
    if contents is None:
        return no_update, no_update

    import base64
    import tempfile

    from app.pipeline.biom_import import import_biom

    # Decode base64 content
    try:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
    except Exception:
        return dbc.Alert("Failed to decode uploaded file.", color="danger"), no_update

    # Write to temp file
    suffix = ".biom"
    if filename and "." in filename:
        suffix = "." + filename.rsplit(".", 1)[-1]

    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
        tmp.write(decoded)
        tmp_path = tmp.name

    try:
        dataset_name = (name or "").strip() or (filename or "BIOM Import")
        result = import_biom(tmp_path, name=dataset_name)
    except Exception as e:
        Path(tmp_path).unlink(missing_ok=True)
        return dbc.Alert(f"BIOM import failed: {e}", color="danger"), no_update

    Path(tmp_path).unlink(missing_ok=True)

    tax_note = "Taxonomy found" if result["has_taxonomy"] else "No taxonomy (run Taxonomy before Pathogen detection)"
    region = result["variable_region"] or "Unknown"

    alert = dbc.Alert(
        [
            html.Strong("BIOM import successful! "),
            html.Br(),
            f"{result['sample_count']} samples, {result['asv_count']} ASVs, "
            f"region: {region}. {tax_note}.",
        ],
        color="success",
        dismissable=True,
    )

    return alert, _build_history_table(checked_history)


# ══════════════════════════════════════════════════════════════════════════════
# Callbacks — Run Taxonomy for BIOM import
# ══════════════════════════════════════════════════════════════════════════════


@dash_app.callback(
    Output("tax-progress-area", "children", allow_duplicate=True),
    Output("tax-poll", "disabled", allow_duplicate=True),
    Input("btn-run-taxonomy-biom", "n_clicks"),
    State("biom-taxonomy-dataset-id", "data"),
    State("tax-threads", "value"),
    prevent_initial_call=True,
)
def launch_taxonomy_for_biom(n_clicks, dataset_id, tax_threads):
    """Launch taxonomy assignment for a BIOM-imported dataset."""
    if not n_clicks or not dataset_id:
        return no_update, no_update

    from app.pipeline.runner import launch_taxonomy_for_dataset

    try:
        launch_taxonomy_for_dataset(dataset_id, threads=tax_threads or 4)
    except ValueError as e:
        return dbc.Alert(str(e), color="danger"), True

    progress = dbc.Card(dbc.CardBody([
        html.H6("Taxonomy Assignment Running..."),
        dbc.Progress(id="biom-tax-progress-bar", value=0, striped=True, animated=True),
        html.P(id="biom-tax-step-label", className="mt-2 text-muted small"),
    ]), className="mb-3")

    return progress, False


