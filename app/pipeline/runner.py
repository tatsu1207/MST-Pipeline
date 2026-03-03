"""
MST-Pipeline — Background pipeline orchestrator.

Launches DADA2, SourceTracker, and Pathogen pipelines in background threads
and tracks progress via status.json files.
"""
import json
import logging
import os
import signal
import subprocess
import threading
from datetime import datetime
from pathlib import Path

from app.config import DADA2_DEFAULTS, DADA2_LONG_READ_DEFAULTS, DATASET_DIR, MAX_THREADS, is_long_read, to_absolute, to_relative

# Active pipeline threads keyed by dataset_id or "st-{run_id}" or "pathogen-{run_id}"
_running_pipelines: dict[int | str, threading.Thread] = {}

# Cancellation infrastructure
_cancel_events: dict[int | str, threading.Event] = {}
_active_procs: dict[int | str, subprocess.Popen] = {}


class PipelineCancelled(Exception):
    pass


def _check_cancel(key: int | str):
    event = _cancel_events.get(key)
    if event and event.is_set():
        raise PipelineCancelled()


def cancel_pipeline(dataset_id: int) -> bool:
    """Cancel a running DADA2 pipeline."""
    event = _cancel_events.get(dataset_id)
    if not event:
        return False
    event.set()
    proc = _active_procs.get(dataset_id)
    if proc and proc.poll() is None:
        def _kill():
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except (OSError, ProcessLookupError):
                pass
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except (OSError, ProcessLookupError):
                    pass
        threading.Thread(target=_kill, daemon=True).start()
    return True


def _update_status(output_dir: Path, step: str, pct: int, completed: list[str], **extra):
    status_file = output_dir / "status.json"
    data = {
        "current_step": step,
        "progress_pct": pct,
        "steps_completed": completed,
        "updated_at": datetime.utcnow().isoformat(),
    }
    if status_file.exists():
        try:
            existing = json.loads(status_file.read_text())
            data["started_at"] = existing.get("started_at")
        except Exception:
            pass
    else:
        data["started_at"] = datetime.utcnow().isoformat()
    data.update(extra)
    status_file.write_text(json.dumps(data, indent=2))


def get_pipeline_status(dataset_id: int) -> dict:
    """Get current status for a DADA2 pipeline."""
    from app.db.database import SessionLocal
    from app.db.models import Dataset

    output_dir = DATASET_DIR / str(dataset_id)
    status_file = output_dir / "status.json"

    step_info = {"current_step": None, "progress_pct": 0, "steps_completed": []}
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    log_tail = ""
    log_file = output_dir / "pipeline.log"
    if log_file.exists():
        try:
            lines = log_file.read_text().splitlines()
            log_tail = "\n".join(lines[-50:])
        except Exception:
            pass

    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        db_status = dataset.status if dataset else "unknown"
    finally:
        db.close()

    thread = _running_pipelines.get(dataset_id)
    thread_alive = thread is not None and thread.is_alive()

    if db_status == "processing" and not thread_alive:
        db = SessionLocal()
        try:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if dataset and dataset.status == "processing":
                if (output_dir / "samples").is_dir():
                    dataset.status = "complete"
                    db_status = "complete"
                else:
                    dataset.status = "failed"
                    db_status = "failed"
                db.commit()
        finally:
            db.close()

    return {
        "dataset_id": dataset_id,
        "status": db_status,
        "current_step": step_info.get("current_step"),
        "progress_pct": step_info.get("progress_pct", 0),
        "steps_completed": step_info.get("steps_completed", []),
        "log_tail": log_tail,
        "thread_alive": thread_alive,
    }


# -- DADA2 Pipeline ------------------------------------------------------------

def launch_dada2_pipeline(
    file_ids: list[int],
    name: str,
    threads: int | None = None,
    error_model: str = "default",
) -> int:
    """Create a Dataset record and launch DADA2 in a background thread.

    Accepts a list of FastqFile IDs (selected by user on the DADA2 page).
    Truncation lengths are auto-detected from quality profiles during the run.
    No taxonomy step. After BIOM conversion, auto-extracts V4 region.
    Returns dataset_id.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, Upload, FastqFile

    db = SessionLocal()
    try:
        fastq_files = db.query(FastqFile).filter(FastqFile.id.in_(file_ids)).all()
        if not fastq_files:
            raise ValueError("No FASTQ files found for the given IDs")

        # Determine sequencing type and variable region from the uploads
        upload_ids = set(f.upload_id for f in fastq_files)
        uploads = db.query(Upload).filter(Upload.id.in_(upload_ids)).all()
        upload_map = {u.id: u for u in uploads}

        seq_types = set(u.sequencing_type for u in uploads if u.sequencing_type)
        regions = set(u.variable_region for u in uploads if u.variable_region)
        seq_type = seq_types.pop() if len(seq_types) == 1 else "paired-end"
        variable_region = regions.pop() if len(regions) == 1 else None

        # Auto-compute threads: samples x 2 capped at CPUs - 1
        n_samples = len(set(f.sample_name for f in fastq_files))
        if threads is None:
            max_cpus = max(1, (os.cpu_count() or 1) - 1)
            threads = min(max(1, n_samples * 2), max_cpus)

        dataset = Dataset(
            upload_id=None,
            name=name,
            status="pending",
            sequencing_type=seq_type,
            variable_region=variable_region,
            trunc_len_f=0,
            trunc_len_r=0,
            error_model=error_model,
        )
        db.add(dataset)
        db.flush()
        dataset_id = dataset.id

        output_dir = DATASET_DIR / str(dataset_id)
        output_dir.mkdir(parents=True, exist_ok=True)
        dataset.pipeline_log_path = to_relative(output_dir / "pipeline.log")
        db.commit()
    finally:
        db.close()

    _cancel_events[dataset_id] = threading.Event()
    resolved_threads = threads

    t = threading.Thread(
        target=_run_dada2_pipeline,
        args=(dataset_id, file_ids, resolved_threads),
        daemon=True,
        name=f"dada2-{dataset_id}",
    )
    _running_pipelines[dataset_id] = t
    t.start()
    return dataset_id


def _run_dada2_pipeline(dataset_id: int, file_ids: list[int], threads: int):
    """Execute full DADA2 pipeline in background thread.

    Uses file_ids (FastqFile IDs) selected by the user, not upload_id.
    Steps: cutadapt → auto-trunc → DADA2 → BIOM.
    V4 extraction from BIOM is done separately after completion.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, FastqFile, Sample

    output_dir = DATASET_DIR / str(dataset_id)
    completed_steps = []

    logger = logging.getLogger(f"dada2-{dataset_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(output_dir / "pipeline.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        dataset.status = "processing"
        db.commit()

        fastq_files = db.query(FastqFile).filter(FastqFile.id.in_(file_ids)).all()

        seq_type = dataset.sequencing_type or "paired-end"
        logger.info(f"Pipeline started for dataset {dataset_id}")
        logger.info(f"  Type: {seq_type}, Files: {len(fastq_files)}, Threads: {threads}")

        # Build sample paths
        sample_paths = {}
        for ff in fastq_files:
            sname = ff.sample_name
            abs_path = str(to_absolute(ff.file_path))
            if sname not in sample_paths:
                sample_paths[sname] = {}
            if ff.read_direction == "R1" or ff.read_direction == "single":
                sample_paths[sname]["R1"] = abs_path
            elif ff.read_direction == "R2":
                sample_paths[sname]["R2"] = abs_path

        # Step 1: Cutadapt primer trimming
        _check_cancel(dataset_id)
        _update_status(output_dir, "cutadapt", 10, completed_steps)
        logger.info("Running cutadapt primer trimming...")

        # Build a temporary directory with symlinks to the selected FASTQ files
        fastq_dir = output_dir / "input_fastq"
        fastq_dir.mkdir(parents=True, exist_ok=True)
        for ff in fastq_files:
            abs_fp = to_absolute(ff.file_path)
            link = fastq_dir / abs_fp.name
            if not link.exists():
                os.symlink(abs_fp, link)

        from app.pipeline.trim import run_cutadapt
        trim_result = run_cutadapt(
            fastq_dir=fastq_dir,
            output_dir=output_dir,
            sequencing_type=seq_type,
            variable_region=dataset.variable_region,
            logger=logger,
            threads=threads,
        )
        trim_dir = trim_result["trimmed_dir"]

        completed_steps.append("cutadapt")

        # Step 3: Auto-detect truncation parameters from quality profiles
        _check_cancel(dataset_id)
        _update_status(output_dir, "auto_trunc", 30, completed_steps)
        logger.info("Auto-detecting truncation parameters from quality profiles...")

        from app.pipeline.quality import detect_truncation_params
        trimmed_dir = Path(trim_dir)

        # Select defaults based on long-read vs short-read mode
        variable_region = dataset.variable_region
        if is_long_read(variable_region):
            defaults = DADA2_LONG_READ_DEFAULTS
            error_model = dataset.error_model or "PacBio"
            logger.info(f"Long-read mode: error_model={error_model}")
        else:
            defaults = DADA2_DEFAULTS
            error_model = dataset.error_model or "default"

        min_overlap = defaults["min_overlap"]
        trim_left_f = defaults["trim_left_f"]
        trim_left_r = defaults["trim_left_r"]

        auto = detect_truncation_params(
            trimmed_dir=trimmed_dir,
            sequencing_type=seq_type,
            variable_region=variable_region,
            min_overlap=min_overlap,
            logger=logger,
        )
        trunc_f = auto["trunc_len_f"]
        trunc_r = auto["trunc_len_r"]
        dataset.trunc_len_f = trunc_f
        dataset.trunc_len_r = trunc_r
        db.commit()
        logger.info(f"Auto-detected truncation: {auto['details']}")
        _update_status(
            output_dir, "auto_trunc", 35, completed_steps,
            auto_trunc_len_f=trunc_f,
            auto_trunc_len_r=trunc_r,
            auto_trunc_details=auto["details"],
        )
        completed_steps.append("auto_trunc")

        # Step 4: DADA2
        _check_cancel(dataset_id)
        _update_status(output_dir, "dada2", 40, completed_steps)
        logger.info("Running DADA2...")

        from app.pipeline.dada2 import run_dada2

        # Long-read specific DADA2 parameters
        band_size = defaults.get("band_size", 16)
        homopolymer_gap_penalty = defaults.get("homopolymer_gap_penalty", 0)
        max_ee = defaults.get("max_ee", 5.0)
        min_len = defaults.get("min_len", 0)
        max_len = defaults.get("max_len", 0)

        dada2_result = run_dada2(
            input_dir=trimmed_dir,
            output_dir=output_dir,
            sequencing_type=seq_type,
            trim_left_f=trim_left_f,
            trim_left_r=trim_left_r,
            trunc_len_f=trunc_f,
            trunc_len_r=trunc_r,
            min_overlap=min_overlap,
            threads=threads,
            logger=logger,
            proc_callback=lambda proc: _active_procs.__setitem__(dataset_id, proc),
            error_model=error_model,
            band_size=band_size,
            homopolymer_gap_penalty=homopolymer_gap_penalty,
            max_ee=max_ee,
            min_len=min_len,
            max_len=max_len,
        )

        _check_cancel(dataset_id)
        completed_steps.append("dada2")

        # Step 5: Taxonomy assignment
        _update_status(output_dir, "taxonomy", 65, completed_steps)
        logger.info("Running taxonomy assignment...")

        from app.pipeline.dada2 import run_taxonomy
        taxonomy_result = run_taxonomy(
            asv_table_path=Path(dada2_result["asv_table_path"]),
            output_dir=output_dir / "dada2",
            threads=threads,
            logger=logger,
            proc_callback=lambda proc: _active_procs.__setitem__(dataset_id, proc),
        )
        dataset.taxonomy_path = to_relative(taxonomy_result["taxonomy_path"])
        db.commit()

        _check_cancel(dataset_id)
        completed_steps.append("taxonomy")
        _update_status(output_dir, "split_samples", 85, completed_steps)

        # Step 6: Split into per-sample ASV tables
        logger.info("Splitting ASV table into per-sample files...")
        import pandas as pd

        tsv_path = Path(dada2_result["asv_table_path"])
        asv_df = pd.read_csv(tsv_path, sep="\t")
        asv_df = asv_df.set_index("sequence").drop(columns=["ASV_ID"], errors="ignore")

        samples_dir = output_dir / "samples"
        samples_dir.mkdir(parents=True, exist_ok=True)

        # Update dataset
        dataset.asv_table_path = None
        dataset.asv_count = dada2_result["asv_count"]
        dataset.sample_count = dada2_result["sample_count"]

        # Create Sample rows with per-sample csv.gz files
        for sname in asv_df.columns:
            col = asv_df[sname]
            col_nonzero = col[col > 0]
            sample_path = samples_dir / f"{sname}.csv.gz"
            col_nonzero.to_frame(name=sname).to_csv(
                sample_path, compression="gzip"
            )

            track = dada2_result["track_reads"].get(sname, {})
            sample = Sample(
                dataset_id=dataset_id,
                sample_name=sname,
                read_count_raw=track.get("input"),
                read_count_filtered=track.get("filtered"),
                read_count_nonchimeric=track.get("nonchim"),
                asv_count=len(col_nonzero),
                asv_table_path=to_relative(sample_path),
            )
            db.add(sample)

        logger.info(f"  Wrote {len(asv_df.columns)} per-sample csv.gz files")

        dataset.status = "complete"
        db.commit()
        completed_steps.append("split_samples")
        _update_status(output_dir, "complete", 100, completed_steps)
        logger.info("Pipeline completed successfully!")

    except PipelineCancelled:
        logger.info("Pipeline cancelled by user.")
        try:
            dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
            if dataset:
                dataset.status = "cancelled"
                db.commit()
        except Exception:
            pass
        _update_status(output_dir, "cancelled", 0, completed_steps)

    except Exception as e:
        cancel_event = _cancel_events.get(dataset_id)
        if cancel_event and cancel_event.is_set():
            logger.info("Pipeline cancelled by user.")
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if dataset:
                    dataset.status = "cancelled"
                    db.commit()
            except Exception:
                pass
        else:
            logger.error(f"Pipeline failed: {e}", exc_info=True)
            try:
                dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
                if dataset:
                    dataset.status = "failed"
                    db.commit()
            except Exception:
                pass
            _update_status(output_dir, "failed", 0, completed_steps)
    finally:
        db.close()
        _running_pipelines.pop(dataset_id, None)
        _cancel_events.pop(dataset_id, None)
        _active_procs.pop(dataset_id, None)
        logger.handlers.clear()


# -- SourceTracker Pipeline ----------------------------------------------------

def launch_sourcetracker_run(
    sample_table_paths: list[str],
    selected_groups: list[str],
    feature_mode: str = "asv",
    src_depth: int = 1000,
    snk_depth: int = 1000,
    restarts: int = 10,
    burnin: int = 100,
    draws: int = 1,
    dataset_id: int | None = None,
    name: str = "SourceTracker Run",
    threads: int = 4,
) -> int:
    """Create a SourcetrackerRun record and launch in background.
    Returns run_id.
    """
    import json as _json
    from app.db.database import SessionLocal
    from app.db.models import SourcetrackerRun

    db = SessionLocal()
    try:
        run = SourcetrackerRun(
            dataset_id=dataset_id,
            name=name,
            status="pending",
            input_biom_path=None,
            selected_groups=_json.dumps(selected_groups),
            feature_mode=feature_mode,
            src_depth=src_depth,
            snk_depth=snk_depth,
            restarts=restarts,
            burnin=burnin,
            draws=draws,
        )
        db.add(run)
        db.flush()
        run_id = run.id

        run_dir = DATASET_DIR / f"st_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)
        run.output_dir = to_relative(run_dir)
        db.commit()
    finally:
        db.close()

    key = f"st-{run_id}"
    _cancel_events[key] = threading.Event()

    t = threading.Thread(
        target=_run_sourcetracker,
        args=(run_id, sample_table_paths, threads),
        daemon=True,
        name=f"st-{run_id}",
    )
    _running_pipelines[key] = t
    t.start()
    return run_id


def get_sourcetracker_status(run_id: int) -> dict:
    from app.db.database import SessionLocal
    from app.db.models import SourcetrackerRun

    run_dir = DATASET_DIR / f"st_{run_id}"
    status_file = run_dir / "status.json"

    step_info = {}
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    log_tail = ""
    log_file = run_dir / "pipeline.log"
    if log_file.exists():
        try:
            lines = log_file.read_text().splitlines()
            log_tail = "\n".join(lines[-30:])
        except Exception:
            pass

    db = SessionLocal()
    try:
        run = db.query(SourcetrackerRun).filter(SourcetrackerRun.id == run_id).first()
        db_status = run.status if run else "unknown"
    finally:
        db.close()

    return {
        "run_id": run_id,
        "status": db_status,
        "current_step": step_info.get("current_step"),
        "progress_pct": step_info.get("progress_pct", 0),
        "log_tail": log_tail,
        "alignment_matched": step_info.get("alignment_matched"),
        "alignment_total": step_info.get("alignment_total"),
        "sample_alignment_pct": step_info.get("sample_alignment_pct", {}),
        "low_alignment_samples": step_info.get("low_alignment_samples", []),
    }


def _run_sourcetracker(run_id: int, sample_table_paths: list[str], threads: int = 4):
    """Execute SourceTracker pipeline in background thread."""
    import json as _json
    import pandas as pd
    from app.db.database import SessionLocal
    from app.db.models import SourcetrackerRun

    run_dir = DATASET_DIR / f"st_{run_id}"
    key = f"st-{run_id}"

    logger = logging.getLogger(f"st-{run_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(run_dir / "pipeline.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    db = SessionLocal()
    try:
        run = db.query(SourcetrackerRun).filter(SourcetrackerRun.id == run_id).first()
        run.status = "processing"
        db.commit()

        selected_groups = _json.loads(run.selected_groups)
        logger.info(f"SourceTracker run {run_id} started")
        logger.info(f"  Groups: {selected_groups}")
        logger.info(f"  Mode: {run.feature_mode}, src_depth={run.src_depth}, snk_depth={run.snk_depth}")

        from app.pipeline.sourcetracker import (
            load_csv_gz, load_fasta, load_design,
            align_features, collapse_by_group, run_gibbs_subprocess,
            extract_v4_from_full_length,
        )
        from app.config import SOURCE_TABLE, SOURCE_FASTA, SOURCE_DESIGN

        _update_status(run_dir, "loading", 5, [])

        # Load and merge per-sample csv.gz files
        dfs = [load_csv_gz(p) for p in sample_table_paths]
        sink_df = pd.concat(dfs, axis=1).fillna(0).astype(int)
        logger.info(f"  Merged {len(sample_table_paths)} sample tables → {sink_df.shape[0]} ASVs x {sink_df.shape[1]} samples")

        # If the dataset used full-length 16S, extract V4 sub-region before alignment
        from app.db.models import Dataset
        dataset_vregion = None
        if run.dataset_id:
            ds = db.query(Dataset).filter(Dataset.id == run.dataset_id).first()
            if ds:
                dataset_vregion = ds.variable_region
        if is_long_read(dataset_vregion):
            logger.info("Full-length 16S detected — extracting V4 sub-region for source alignment...")
            _update_status(run_dir, "v4_extraction", 10, [])
            sink_df = extract_v4_from_full_length(sink_df, logger=logger)

        source_df = load_csv_gz(str(SOURCE_TABLE))
        db_fasta = load_fasta(str(SOURCE_FASTA))
        design = load_design(str(SOURCE_DESIGN))

        design_filtered = {s: g for s, g in design.items() if g in selected_groups}

        _update_status(run_dir, "aligning", 20, [])
        logger.info("Aligning features...")

        _check_cancel(key)
        sink_al, src_al, (n_m, n_t) = align_features(
            sink_df, source_df, db_fasta,
            mode=run.feature_mode, threads=threads,
        )
        logger.info(f"  Matched {n_m}/{n_t} features")

        # Per-sample read retention after alignment
        total_before = sink_df.sum(axis=0)          # reads per sample before
        total_after  = sink_al.sum(axis=1)           # reads per sample after (sink_al is samples x features)
        sample_pct = {}
        low_samples = []
        for sample in total_before.index:
            before = int(total_before[sample])
            after  = int(total_after.get(sample, 0))
            pct = (after / before * 100) if before > 0 else 0
            sample_pct[sample] = round(pct, 1)
            logger.info(f"  Alignment — {sample}: {after}/{before} reads retained ({pct:.1f}%)")
            if pct < 10:
                low_samples.append(sample)

        if low_samples:
            logger.warning(
                f"Low alignment rate for {len(low_samples)} sample(s): "
                + ", ".join(low_samples)
                + ". Excluding from Gibbs sampling."
            )
            sink_al = sink_al.drop(index=low_samples, errors="ignore")

        if sink_al.empty:
            raise ValueError(
                "All samples had <10% read alignment to the source database. "
                "Cannot run SourceTracker."
            )

        _update_status(run_dir, "aligning", 30, ["aligning"],
                       alignment_matched=n_m, alignment_total=n_t,
                       sample_alignment_pct=sample_pct,
                       low_alignment_samples=low_samples)

        _update_status(run_dir, "collapsing", 35, ["aligning"])
        logger.info("Collapsing by group...")
        src_col = collapse_by_group(src_al, design_filtered)

        _update_status(run_dir, "gibbs", 40, ["aligning", "collapsing"])
        logger.info(f"Running Gibbs sampling on {len(sink_al)} sample(s)...")

        _check_cancel(key)
        mp, mps = run_gibbs_subprocess(
            src_col, sink_al,
            src_depth=run.src_depth,
            snk_depth=run.snk_depth,
            restarts=run.restarts,
            burnin=run.burnin,
            draws=run.draws,
            status_fn=lambda msg: logger.info(f"  Gibbs: {msg}"),
        )

        # Save results
        mp.to_csv(run_dir / "proportions.csv")
        mps.to_csv(run_dir / "stds.csv")

        run.proportions_path = to_relative(run_dir / "proportions.csv")
        run.stds_path = to_relative(run_dir / "stds.csv")
        run.status = "complete"
        db.commit()

        _update_status(run_dir, "complete", 100, ["aligning", "collapsing", "gibbs"],
                       alignment_matched=n_m, alignment_total=n_t,
                       sample_alignment_pct=sample_pct,
                       low_alignment_samples=low_samples)
        logger.info("SourceTracker completed successfully!")

    except PipelineCancelled:
        logger.info("SourceTracker cancelled.")
        try:
            run = db.query(SourcetrackerRun).filter(SourcetrackerRun.id == run_id).first()
            if run:
                run.status = "failed"
                db.commit()
        except Exception:
            pass
    except Exception as e:
        logger.error(f"SourceTracker failed: {e}", exc_info=True)
        try:
            run = db.query(SourcetrackerRun).filter(SourcetrackerRun.id == run_id).first()
            if run:
                run.status = "failed"
                db.commit()
        except Exception:
            pass
        _update_status(run_dir, "failed", 0, [])
    finally:
        db.close()
        _running_pipelines.pop(key, None)
        _cancel_events.pop(key, None)
        logger.handlers.clear()


# -- Taxonomy-only Pipeline (for BIOM imports) ---------------------------------

def launch_taxonomy_for_dataset(dataset_id: int, threads: int = 4) -> int:
    """Launch taxonomy assignment for an existing Dataset that lacks taxonomy.

    Collects all unique ASV sequences from per-sample csv.gz files,
    writes an asv_table.tsv + rep_seqs.fasta, then runs the taxonomy R script.
    Returns dataset_id.
    """
    from app.db.database import SessionLocal
    from app.db.models import Dataset, Sample

    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found")
        if dataset.taxonomy_path:
            raise ValueError(f"Dataset {dataset_id} already has taxonomy")
        if dataset.status != "complete":
            raise ValueError(f"Dataset {dataset_id} is not complete (status={dataset.status})")

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset_id, Sample.asv_table_path.isnot(None))
            .all()
        )
        if not samples:
            raise ValueError(f"No sample ASV tables found for dataset {dataset_id}")

        sample_paths = [str(to_absolute(s.asv_table_path)) for s in samples]
    finally:
        db.close()

    key = f"tax-{dataset_id}"
    _cancel_events[key] = threading.Event()

    t = threading.Thread(
        target=_run_taxonomy_for_dataset,
        args=(dataset_id, sample_paths, threads),
        daemon=True,
        name=key,
    )
    _running_pipelines[key] = t
    t.start()
    return dataset_id


def get_taxonomy_status(dataset_id: int) -> dict:
    """Get current status for a taxonomy-only pipeline."""
    from app.db.database import SessionLocal
    from app.db.models import Dataset

    output_dir = DATASET_DIR / str(dataset_id)
    status_file = output_dir / "taxonomy_status.json"

    step_info = {}
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        has_taxonomy = bool(dataset.taxonomy_path) if dataset else False
    finally:
        db.close()

    key = f"tax-{dataset_id}"
    thread = _running_pipelines.get(key)
    thread_alive = thread is not None and thread.is_alive()

    status = "idle"
    if thread_alive:
        status = "running"
    elif has_taxonomy:
        status = "complete"
    elif step_info.get("current_step") == "failed":
        status = "failed"

    return {
        "dataset_id": dataset_id,
        "status": status,
        "current_step": step_info.get("current_step"),
        "progress_pct": step_info.get("progress_pct", 0),
        "has_taxonomy": has_taxonomy,
        "thread_alive": thread_alive,
    }


def _run_taxonomy_for_dataset(dataset_id: int, sample_paths: list[str], threads: int):
    """Run taxonomy assignment in a background thread for an imported dataset."""
    import pandas as pd
    from app.db.database import SessionLocal
    from app.db.models import Dataset
    from app.pipeline.sourcetracker import load_csv_gz

    output_dir = DATASET_DIR / str(dataset_id)
    key = f"tax-{dataset_id}"

    logger = logging.getLogger(key)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(output_dir / "taxonomy.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    # Use a separate status file to avoid clobbering DADA2 status.json
    def _tax_status(step, pct):
        status_file = output_dir / "taxonomy_status.json"
        status_file.write_text(json.dumps({
            "current_step": step,
            "progress_pct": pct,
        }))

    db = SessionLocal()
    try:
        _tax_status("loading", 10)
        logger.info(f"Taxonomy pipeline started for dataset {dataset_id}")

        # Load and merge all per-sample csv.gz to collect unique sequences
        dfs = [load_csv_gz(p) for p in sample_paths]
        merged = pd.concat(dfs, axis=1).fillna(0).astype(int)
        sequences = list(merged.index)
        logger.info(f"Collected {len(sequences)} unique ASV sequences")

        # Write asv_table.tsv + rep_seqs.fasta
        dada2_dir = output_dir / "dada2"
        dada2_dir.mkdir(parents=True, exist_ok=True)

        asv_rows = []
        with open(dada2_dir / "rep_seqs.fasta", "w") as fasta:
            for i, seq in enumerate(sequences, 1):
                asv_id = f"ASV_{i}"
                fasta.write(f">{asv_id}\n{seq}\n")
                row = {"ASV_ID": asv_id, "sequence": seq}
                for col in merged.columns:
                    row[col] = int(merged.loc[seq, col])
                asv_rows.append(row)

        asv_df = pd.DataFrame(asv_rows)
        asv_table_path = dada2_dir / "asv_table.tsv"
        asv_df.to_csv(asv_table_path, sep="\t", index=False)
        logger.info(f"Wrote {len(sequences)} sequences to {asv_table_path}")

        # Run taxonomy R script
        _check_cancel(key)
        _tax_status("taxonomy", 30)
        logger.info("Running taxonomy assignment...")

        from app.pipeline.dada2 import run_taxonomy
        taxonomy_result = run_taxonomy(
            asv_table_path=asv_table_path,
            output_dir=dada2_dir,
            threads=threads,
            logger=logger,
            proc_callback=lambda proc: _active_procs.__setitem__(key, proc),
        )

        # Update dataset
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        dataset.taxonomy_path = to_relative(taxonomy_result["taxonomy_path"])
        db.commit()

        _tax_status("complete", 100)
        logger.info("Taxonomy assignment complete!")

    except PipelineCancelled:
        logger.info("Taxonomy pipeline cancelled.")
        _tax_status("failed", 0)
    except Exception as e:
        logger.error(f"Taxonomy pipeline failed: {e}", exc_info=True)
        _tax_status("failed", 0)
    finally:
        db.close()
        _running_pipelines.pop(key, None)
        _cancel_events.pop(key, None)
        _active_procs.pop(key, None)
        logger.handlers.clear()


# -- Pathogen Pipeline ---------------------------------------------------------

def launch_pathogen_run(
    sample_table_paths: list[str],
    silva_identity: float = 0.97,
    dataset_id: int | None = None,
    name: str = "Pathogen Run",
    threads: int = 4,
) -> int:
    """Create a PathogenRun record and launch in background. Returns run_id."""
    from app.db.database import SessionLocal
    from app.db.models import PathogenRun

    db = SessionLocal()
    try:
        run = PathogenRun(
            dataset_id=dataset_id,
            name=name,
            status="pending",
            input_biom_path=None,
            silva_identity=silva_identity,
        )
        db.add(run)
        db.flush()
        run_id = run.id

        run_dir = DATASET_DIR / f"pathogen_{run_id}"
        run_dir.mkdir(parents=True, exist_ok=True)
        run.output_dir = to_relative(run_dir)
        db.commit()
    finally:
        db.close()

    key = f"pathogen-{run_id}"
    _cancel_events[key] = threading.Event()

    t = threading.Thread(
        target=_run_pathogen,
        args=(run_id, sample_table_paths, threads),
        daemon=True,
        name=f"pathogen-{run_id}",
    )
    _running_pipelines[key] = t
    t.start()
    return run_id


def get_pathogen_status(run_id: int) -> dict:
    from app.db.database import SessionLocal
    from app.db.models import PathogenRun

    run_dir = DATASET_DIR / f"pathogen_{run_id}"
    status_file = run_dir / "status.json"

    step_info = {}
    if status_file.exists():
        try:
            step_info = json.loads(status_file.read_text())
        except Exception:
            pass

    log_tail = ""
    log_file = run_dir / "pipeline.log"
    if log_file.exists():
        try:
            lines = log_file.read_text().splitlines()
            log_tail = "\n".join(lines[-30:])
        except Exception:
            pass

    db = SessionLocal()
    try:
        run = db.query(PathogenRun).filter(PathogenRun.id == run_id).first()
        db_status = run.status if run else "unknown"
    finally:
        db.close()

    return {
        "run_id": run_id,
        "status": db_status,
        "current_step": step_info.get("current_step"),
        "progress_pct": step_info.get("progress_pct", 0),
        "log_tail": log_tail,
    }


def _run_pathogen(run_id: int, sample_table_paths: list[str], threads: int = 4):
    """Execute pathogen detection pipeline in background thread."""
    import pandas as pd
    from app.db.database import SessionLocal
    from app.db.models import PathogenRun

    run_dir = DATASET_DIR / f"pathogen_{run_id}"
    key = f"pathogen-{run_id}"

    logger = logging.getLogger(f"pathogen-{run_id}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh_handler = logging.FileHandler(run_dir / "pipeline.log")
    fh_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh_handler)

    db = SessionLocal()
    try:
        run = db.query(PathogenRun).filter(PathogenRun.id == run_id).first()
        run.status = "processing"
        db.commit()

        logger.info(f"Pathogen detection run {run_id} started")

        from app.pipeline.sourcetracker import load_csv_gz
        from app.pipeline.pathogen import detect, build_summary

        _update_status(run_dir, "loading", 10, [])

        # Load and merge per-sample csv.gz files
        dfs = [load_csv_gz(p) for p in sample_table_paths]
        asv_df = pd.concat(dfs, axis=1).fillna(0).astype(int)
        logger.info(f"  Merged {len(sample_table_paths)} sample tables → {asv_df.shape[0]} ASVs x {asv_df.shape[1]} samples")

        sequences = list(asv_df.index)
        logger.info(f"Loaded {len(sequences)} ASVs x {asv_df.shape[1]} samples")

        _check_cancel(key)
        _update_status(run_dir, "classifying", 30, [])

        # Load pre-computed taxonomy from DADA2 pipeline
        from app.db.models import Dataset
        taxonomy_path = None
        if run.dataset_id:
            ds = db.query(Dataset).filter(Dataset.id == run.dataset_id).first()
            if ds and ds.taxonomy_path:
                taxonomy_path = to_absolute(ds.taxonomy_path)

        if not taxonomy_path or not taxonomy_path.exists():
            raise RuntimeError("No taxonomy data found. Please re-run DADA2 for this dataset.")

        logger.info("Loading pre-computed taxonomy...")
        tax_df = pd.read_csv(taxonomy_path, sep="\t")

        # Build ASV_ID → sequence mapping from rep_seqs.fasta
        rep_seqs_path = taxonomy_path.parent / "rep_seqs.fasta"
        asv_to_seq = {}
        with open(rep_seqs_path) as fh:
            asv_id = None
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    asv_id = line[1:]
                elif asv_id:
                    asv_to_seq[asv_id] = line
                    asv_id = None

        # Build sequence → genus mapping
        seq_to_genus = {}
        for _, row in tax_df.iterrows():
            seq = asv_to_seq.get(row["ASV_ID"])
            if seq:
                genus = row.get("Genus")
                seq_to_genus[seq] = genus if pd.notna(genus) else None
        n_classified = sum(1 for g in seq_to_genus.values() if g)
        logger.info(f"Loaded {n_classified}/{len(sequences)} genus assignments")

        _check_cancel(key)
        _update_status(run_dir, "detecting", 70, ["classifying"])
        logger.info("Detecting pathogens...")

        count_df, ra_df = detect(asv_df, seq_to_genus)

        if ra_df.empty:
            logger.info("No pathogenic bacteria detected")
        else:
            logger.info(f"Detected {ra_df.shape[1]} pathogenic genera")
            ra_df.to_csv(run_dir / "pathogen_ra.csv")
            count_df.to_csv(run_dir / "pathogen_counts.csv")

            summary = build_summary(ra_df)
            import pandas as pd
            pd.DataFrame(summary).to_csv(run_dir / "pathogen_summary.csv", index=False)

        run.status = "complete"
        db.commit()
        _update_status(run_dir, "complete", 100, ["classifying", "detecting"])
        logger.info("Pathogen detection completed!")

    except PipelineCancelled:
        logger.info("Pathogen detection cancelled.")
        try:
            run = db.query(PathogenRun).filter(PathogenRun.id == run_id).first()
            if run:
                run.status = "failed"
                db.commit()
        except Exception:
            pass
    except Exception as e:
        logger.error(f"Pathogen detection failed: {e}", exc_info=True)
        try:
            run = db.query(PathogenRun).filter(PathogenRun.id == run_id).first()
            if run:
                run.status = "failed"
                db.commit()
        except Exception:
            pass
        _update_status(run_dir, "failed", 0, [])
    finally:
        db.close()
        _running_pipelines.pop(key, None)
        _cancel_events.pop(key, None)
        logger.handlers.clear()
