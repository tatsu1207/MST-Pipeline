#!/usr/bin/env bash
# ============================================================================
# MST-Pipeline — Run Script
# ============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="mst"

# Ensure conda is available
if ! command -v conda &>/dev/null; then
    for CONDA_PREFIX in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "/opt/conda" "$HOME/.conda"; do
        if [ -f "${CONDA_PREFIX}/etc/profile.d/conda.sh" ]; then
            source "${CONDA_PREFIX}/etc/profile.d/conda.sh"
            break
        fi
    done
fi

if ! command -v conda &>/dev/null; then
    echo "Error: conda not found. Please install Miniconda or Anaconda first."
    exit 1
fi

eval "$(conda shell.bash hook 2>/dev/null)"

# Check if environment exists
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Error: Conda environment '${ENV_NAME}' not found."
    echo "Run ./setup_ubuntu.sh first."
    exit 1
fi

# Activate environment
conda activate "${ENV_NAME}"

# Derive port from UID: 7000 + UID (e.g. UID 1000 → port 8000)
PORT=$((8000 + $(id -u)))

# Kill existing process on our port if running
EXISTING_PID=$(lsof -ti:${PORT} 2>/dev/null || true)
if [ -n "${EXISTING_PID}" ]; then
    echo "Stopping existing process on port ${PORT} (PID: ${EXISTING_PID})..."
    kill ${EXISTING_PID} 2>/dev/null || true
    sleep 1
fi

LOGFILE="${SCRIPT_DIR}/run.log"

echo "Starting MST-Pipeline in background..."
echo "   Open http://localhost:${PORT} in your browser"
echo "   Logs: ${LOGFILE}"
echo "   Stop: kill \$(cat ${SCRIPT_DIR}/run.pid)"
echo ""

cd "${SCRIPT_DIR}"
nohup uvicorn app.main:app --reload --reload-exclude data --host 0.0.0.0 --port ${PORT} >> "${LOGFILE}" 2>&1 &
echo $! > "${SCRIPT_DIR}/run.pid"
echo "Running (PID: $!)"
