#!/usr/bin/env bash
# =============================================================================
# setup_ubuntu.sh  —  MST-Pipeline conda setup for Ubuntu
# =============================================================================
#
# Creates TWO conda environments:
#
#   mst
#     • Python 3.11  : fastapi, uvicorn, dash, dash-bootstrap-components,
#                      dash-uploader, plotly, sqlalchemy, pandas, numpy,
#                      scipy, biopython, biom-format, h5py
#     • CLI tools    : vsearch, bbmap (bbmap.sh + reformat.sh), cutadapt
#     • R + DADA2    : r-base, bioconductor-dada2
#
#   ST
#     • Python 3.9   : sourcetracker (Gibbs sampler), pandas, numpy
#     • Called via:   conda run --no-capture-output -n ST python3 _run_gibbs.py
#
# Usage:
#   bash setup_ubuntu.sh
#
# Requirements:
#   • Miniconda or Anaconda installed
#     https://docs.conda.io/en/latest/miniconda.html
# =============================================================================

set -euo pipefail

ENV_MAIN="mst"
ENV_ST="ST"

# ── Terminal colours ──────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; RESET='\033[0m'

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[ OK ]${RESET}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
die()     { echo -e "${RED}[ERR ]${RESET}  $*" >&2; exit 1; }
header()  { echo -e "\n${BOLD}$*${RESET}"; }

# ── Detect conda / mamba — install mamba if missing ──────────────────────────
if command -v mamba &>/dev/null; then
    CM="mamba"
    info "Using mamba (fast solver)."
elif command -v conda &>/dev/null; then
    info "mamba not found — installing into base environment…"
    conda install -y -n base -c conda-forge mamba
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    if command -v mamba &>/dev/null; then
        CM="mamba"
        success "mamba installed successfully."
    else
        warn "mamba install succeeded but is not on PATH yet — falling back to conda."
        CM="conda"
    fi
else
    die "conda not found.\nInstall Miniconda: https://docs.conda.io/en/latest/miniconda.html"
fi

# ── Helper: remove env if it already exists ───────────────────────────────────
_env_exists() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

_maybe_remove() {
    local env="$1"
    if _env_exists "$env"; then
        warn "Environment '$env' already exists."
        read -rp "  Remove and recreate? [y/N] " ans
        if [[ "$ans" =~ ^[Yy]$ ]]; then
            conda env remove -y -n "$env"
            success "Removed '$env'."
            return 0
        else
            info "Keeping existing '$env' — skipping creation."
            return 1
        fi
    fi
    return 0
}

# =============================================================================
#  STEP 1 — Build the main MST environment
# =============================================================================
header "══════════════════════════════════════════════════════"
header " Step 1/3 — Building environment: ${ENV_MAIN}"
header "══════════════════════════════════════════════════════"

if _maybe_remove "$ENV_MAIN"; then

    # ── 1a. Create env with ALL conda packages in one solve ───────────────────
    # Solving everything at once avoids the std::bad_alloc OOM that happens
    # when adding R+DADA2 into an existing large environment.
    info "Creating '${ENV_MAIN}' — Python 3.11 + R/DADA2 + CLI tools (one-shot solve)…"
    info "This may take several minutes…"
    $CM create -y -n "$ENV_MAIN" \
        -c conda-forge -c bioconda -c defaults \
        python=3.11 \
        pandas \
        numpy \
        scipy \
        biopython \
        h5py \
        vsearch \
        bbmap \
        cutadapt \
        r-base \
        bioconductor-dada2 \
        r-optparse \
        r-jsonlite \
        "tbb=2020.*"
    success "Conda packages installed (Python + R/DADA2 + CLI tools + TBB pinned)."

    # ── 1b. pip packages (web framework + BIOM) ──────────────────────────────
    info "Installing pip packages (fastapi, dash, biom-format, etc.)…"
    conda run -n "$ENV_MAIN" pip install --no-input \
        fastapi \
        "uvicorn[standard]" \
        python-multipart \
        dash \
        dash-bootstrap-components \
        dash-uploader \
        plotly \
        sqlalchemy \
        biom-format
    success "Pip packages installed."

fi  # end _maybe_remove mst

# =============================================================================
#  STEP 2 — Build the ST environment (sourcetracker, Python 3.9)
# =============================================================================
header "══════════════════════════════════════════════════════"
header " Step 2/3 — Building environment: ${ENV_ST}"
header "══════════════════════════════════════════════════════"

if _maybe_remove "$ENV_ST"; then

    info "Creating '${ENV_ST}' with Python 3.9 + sourcetracker…"
    $CM create -y -n "$ENV_ST" \
        -c conda-forge -c bioconda -c defaults \
        python=3.9 \
        pandas \
        numpy \
        sourcetracker
    success "sourcetracker environment installed."

fi  # end _maybe_remove ST

# =============================================================================
#  STEP 3 — Verify
# =============================================================================
header "══════════════════════════════════════════════════════"
header " Step 3/3 — Verification"
header "══════════════════════════════════════════════════════"

# ── Main env Python packages ─────────────────────────────────────────────────
info "Checking Python packages in '${ENV_MAIN}'…"
conda run -n "$ENV_MAIN" python - <<'PYCHECK'
import importlib, sys
failures = []
checks = [
    ("fastapi",                 "fastapi"),
    ("uvicorn",                 "uvicorn"),
    ("dash",                    "dash"),
    ("dash-bootstrap-components","dash_bootstrap_components"),
    ("plotly",                  "plotly"),
    ("sqlalchemy",              "sqlalchemy"),
    ("pandas",                  "pandas"),
    ("numpy",                   "numpy"),
    ("scipy",                   "scipy"),
    ("biopython",               "Bio"),
    ("biom-format",             "biom"),
    ("h5py",                    "h5py"),
]
for label, module in checks:
    try:
        m = importlib.import_module(module)
        ver = getattr(m, "__version__", "ok")
        print(f"  ✓ {label} {ver}")
    except ImportError:
        print(f"  ✗ {label}  MISSING")
        failures.append(label)
if failures:
    sys.exit(1)
PYCHECK
success "Main env Python packages OK."

# ── CLI tools ─────────────────────────────────────────────────────────────────
info "Checking CLI tools in '${ENV_MAIN}'…"
conda run -n "$ENV_MAIN" bash <<'CLICHECK'
ok=0
check() {
    if command -v "$1" &>/dev/null; then
        echo "  ✓ $1"
    else
        echo "  ✗ $1  NOT FOUND"
        ok=1
    fi
}
check bbmap.sh
check reformat.sh
check cutadapt
check vsearch
check Rscript
exit $ok
CLICHECK
success "CLI tools OK."

# ── R packages ───────────────────────────────────────────────────────────────
info "Checking R packages…"
conda run -n "$ENV_MAIN" Rscript - <<'RCHECK'
pkgs <- c("dada2", "optparse", "jsonlite")
ok <- TRUE
for (pkg in pkgs) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("  ✓ %s %s\n", pkg, as.character(packageVersion(pkg))))
    } else {
        cat(sprintf("  ✗ %s  NOT FOUND\n", pkg))
        ok <- FALSE
    }
}
if (!ok) quit(status = 1)
RCHECK
success "R packages OK."

# ── sourcetracker env ─────────────────────────────────────────────────────────
info "Checking sourcetracker in '${ENV_ST}'…"
conda run -n "$ENV_ST" python - <<'STCHECK'
import sys
try:
    from sourcetracker._sourcetracker import gibbs, subsample_dataframe
    print("  ✓ gibbs + subsample_dataframe importable")
except ImportError as e:
    print(f"  ✗ {e}")
    sys.exit(1)
STCHECK
success "sourcetracker Gibbs OK."

# =============================================================================
#  Done — usage instructions
# =============================================================================
echo
echo -e "${GREEN}${BOLD}╔══════════════════════════════════════════════════════════════╗"
echo -e "║                    Setup complete!                          ║"
echo -e "╚══════════════════════════════════════════════════════════════╝${RESET}"
echo
echo -e "  ${BOLD}Run the MST-Pipeline web app:${RESET}"
echo -e "    conda activate ${ENV_MAIN}"
echo -e "    bash run.sh"
echo
echo -e "  ${BOLD}Or without activating:${RESET}"
echo -e "    conda run -n ${ENV_MAIN} uvicorn app.main:app --host 0.0.0.0 --port 8050"
echo
echo -e "  Then open ${CYAN}http://localhost:8050${RESET} in your browser."
echo
