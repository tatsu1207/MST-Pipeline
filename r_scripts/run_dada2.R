#!/usr/bin/env Rscript
# ============================================================================
# MicrobiomeDash — DADA2 Denoising Pipeline
#
# Supports paired-end and single-end modes.
# Outputs: asv_table.tsv, rep_seqs.fasta, track_reads.tsv
# Prints JSON status to stdout on completion.
# ============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(optparse)
  library(jsonlite)
})

# ── CLI Arguments ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--input_dir",   type="character", help="Directory with FASTQ files"),
  make_option("--output_dir",  type="character", help="Output directory"),
  make_option("--mode",        type="character", default="paired", help="paired or single"),
  make_option("--trim_left_f", type="integer",   default=0,   help="Trim N bases from 5' of forward reads"),
  make_option("--trim_left_r", type="integer",   default=0,   help="Trim N bases from 5' of reverse reads"),
  make_option("--trunc_len_f", type="integer",   default=0,   help="Truncate forward reads (0=no truncation)"),
  make_option("--trunc_len_r", type="integer",   default=0,   help="Truncate reverse reads (0=no truncation)"),
  make_option("--min_overlap", type="integer",   default=12,  help="Min overlap for merging (PE only)"),
  make_option("--threads",     type="integer",   default=1,   help="Number of threads")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input_dir) || is.null(opt$output_dir)) {
  stop("--input_dir and --output_dir are required")
}

is_paired <- (opt$mode == "paired")

# ── Main Pipeline ─────────────────────────────────────────────────────────────

tryCatch({

  dir.create(opt$output_dir, recursive=TRUE, showWarnings=FALSE)

  # Discover FASTQ files
  if (is_paired) {
    fnFs <- sort(list.files(opt$input_dir, pattern="_R1([_.].*)?\\.(fastq|fq)(\\.gz)?$|_1\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))
    fnRs <- sort(list.files(opt$input_dir, pattern="_R2([_.].*)?\\.(fastq|fq)(\\.gz)?$|_2\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))

    if (length(fnFs) == 0) stop("No forward (R1) FASTQ files found")
    if (length(fnFs) != length(fnRs)) stop("Mismatched number of forward and reverse files")

    # Extract sample names from forward reads
    sample.names <- sub("_R1([_.].*)?\\.(fastq|fq)(\\.gz)?$|_1\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
    cat("Found", length(fnFs), "paired-end samples\n")
  } else {
    fnFs <- sort(list.files(opt$input_dir, pattern="\\.(fastq|fq)(\\.gz)?$", full.names=TRUE))
    if (length(fnFs) == 0) stop("No FASTQ files found")
    sample.names <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(fnFs))
    # Also strip _R1 suffixes if present in single-end mode
    sample.names <- sub("_R1[_.]?.*$|_1$", "", sample.names)
    cat("Found", length(fnFs), "single-end samples\n")
  }

  # ── Step 1: Filter and Trim ──────────────────────────────────────────────

  filt_dir <- file.path(opt$output_dir, "filtered")
  dir.create(filt_dir, recursive=TRUE, showWarnings=FALSE)

  filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))

  # Use the --threads value for DADA2's multithread parameter.
  # When threads > 1 on Linux, DADA2 uses mclapply (forked parallelism).
  # (On Windows/WSL2 with fork issues, pass --threads 1 from the caller.)
  use_mt <- ifelse(opt$threads > 1, opt$threads, FALSE)

  if (is_paired) {
    filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
    cat("Filtering and trimming (paired-end)...\n")
    filt_out <- filterAndTrim(
      fnFs, filtFs, fnRs, filtRs,
      trimLeft  = c(opt$trim_left_f, opt$trim_left_r),
      truncLen  = c(opt$trunc_len_f, opt$trunc_len_r),
      maxN = 0, maxEE = c(5, 5), truncQ = 2,
      rm.phix = TRUE, compress = TRUE,
      multithread = use_mt
    )
  } else {
    cat("Filtering and trimming (single-end)...\n")
    filt_out <- filterAndTrim(
      fnFs, filtFs,
      trimLeft  = opt$trim_left_f,
      truncLen  = opt$trunc_len_f,
      maxN = 0, maxEE = 5, truncQ = 2,
      rm.phix = TRUE, compress = TRUE,
      multithread = use_mt
    )
  }

  # Remove samples that had 0 reads pass filtering
  keep <- filt_out[, "reads.out"] > 0
  if (sum(keep) == 0) stop("No reads passed filtering for any sample")

  cat(sum(keep), "of", length(keep), "samples passed filtering\n")
  sample.names <- sample.names[keep]
  filtFs <- filtFs[keep]
  if (is_paired) filtRs <- filtRs[keep]

  # ── Step 2: Learn Error Rates ────────────────────────────────────────────

  cat("Learning error rates (forward)...\n")
  errF <- learnErrors(filtFs, multithread=use_mt)

  if (is_paired) {
    cat("Learning error rates (reverse)...\n")
    errR <- learnErrors(filtRs, multithread=use_mt)
  }

  # ── Step 3: Denoise ──────────────────────────────────────────────────────

  cat("Denoising forward reads...\n")
  dadaFs <- dada(filtFs, err=errF, multithread=use_mt)

  if (is_paired) {
    cat("Denoising reverse reads...\n")
    dadaRs <- dada(filtRs, err=errR, multithread=use_mt)
  }

  # ── Step 4: Merge (PE only) ──────────────────────────────────────────────

  if (is_paired) {
    cat("Merging paired reads...\n")
    merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                         minOverlap=opt$min_overlap,
                         maxMismatch=1)
    seqtab <- makeSequenceTable(merged)
  } else {
    seqtab <- makeSequenceTable(dadaFs)
  }

  cat("Sequence table:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs\n")

  # ── Step 5: Remove Chimeras ──────────────────────────────────────────────

  cat("Removing chimeras...\n")
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",
                                       multithread=use_mt)

  asv_count <- ncol(seqtab.nochim)
  sample_count <- nrow(seqtab.nochim)
  cat("After chimera removal:", sample_count, "samples,", asv_count, "ASVs\n")

  if (asv_count == 0) stop("No ASVs survived chimera removal")

  # ── Build Read Tracking Table ────────────────────────────────────────────

  # When there's only 1 sample, dada()/mergePairs() return single objects
  # instead of lists.  Wrap them so sapply() works uniformly.
  if (length(sample.names) == 1) {
    dadaFs <- list(dadaFs)
    if (is_paired) {
      dadaRs <- list(dadaRs)
      merged <- list(merged)
    }
  }

  # Get denoised counts
  get_n <- function(x) sum(getUniques(x))

  if (is_paired) {
    track <- cbind(
      filt_out[keep, , drop=FALSE],
      sapply(dadaFs, get_n),
      sapply(dadaRs, get_n),
      sapply(merged, get_n),
      rowSums(seqtab.nochim)
    )
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  } else {
    track <- cbind(
      filt_out[keep, , drop=FALSE],
      sapply(dadaFs, get_n),
      rowSums(seqtab.nochim)
    )
    colnames(track) <- c("input", "filtered", "denoised", "nonchim")
  }

  rownames(track) <- sample.names

  # ── Write Outputs ────────────────────────────────────────────────────────

  # ASV table: rows = ASVs, columns = samples
  # First column is ASV_ID, second is the actual sequence, then one column per sample
  asv_ids <- paste0("ASV_", seq_len(asv_count))
  seqs <- colnames(seqtab.nochim)

  # Transpose: seqtab.nochim is samples x sequences → we want sequences x samples
  asv_table <- as.data.frame(t(seqtab.nochim))
  colnames(asv_table) <- sample.names
  asv_table <- cbind(ASV_ID = asv_ids, sequence = seqs, asv_table)

  asv_path <- file.path(opt$output_dir, "asv_table.tsv")
  write.table(asv_table, asv_path, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Wrote:", asv_path, "\n")

  # Read tracking table
  track_path <- file.path(opt$output_dir, "track_reads.tsv")
  track_df <- cbind(sample = rownames(track), as.data.frame(track))
  write.table(track_df, track_path, sep="\t", row.names=FALSE, quote=FALSE)
  cat("Wrote:", track_path, "\n")

  # ── Success ──────────────────────────────────────────────────────────────

  status <- toJSON(list(
    status       = "success",
    asv_count    = asv_count,
    sample_count = sample_count
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=0)

}, error = function(e) {
  status <- toJSON(list(
    status  = "error",
    message = conditionMessage(e)
  ), auto_unbox=TRUE)
  cat("\n", status, "\n", sep="")
  quit(status=1)
})
