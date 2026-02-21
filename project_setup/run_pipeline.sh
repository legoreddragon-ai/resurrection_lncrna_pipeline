#!/bin/bash
# =============================================================================
# run_pipeline.sh — Master End-to-End Pipeline
# Resurrection Plant lncRNA Discovery Pipeline
# =============================================================================
# Usage:
#   bash setup/run_pipeline.sh [--skip-download] [--skip-trim] [--skip-align]
#
# Runs all stages in order:
#   Stage 0: Preflight checks
#   Stage 1: Download data
#   Stage 2: Trim reads
#   Stage 3: Align reads
#   Stage 4: StringTie assembly
#   Stage 5: Merge expression matrix
#   Stage 6: Differential expression
#   Stage 7: Visualizations
# =============================================================================

set -uo pipefail

# ── Resolve paths dynamically ─────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

LOGFILE="${LOG_DIR}/pipeline_$(date '+%Y%m%d_%H%M%S').log"

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_DOWNLOAD=false
SKIP_TRIM=false
SKIP_ALIGN=false

for arg in "$@"; do
    case $arg in
        --skip-download) SKIP_DOWNLOAD=true ;;
        --skip-trim)     SKIP_TRIM=true ;;
        --skip-align)    SKIP_ALIGN=true ;;
    esac
done

# ── Logging helper ────────────────────────────────────────────────────────────
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGFILE"
}

log "============================================================"
log "  Resurrection Plant lncRNA Pipeline — Full Run"
log "  Base directory: ${BASE_DIR}"
log "  Log file: ${LOGFILE}"
log "============================================================"

# ── Stage 0: Preflight ───────────────────────────────────────────────────────
log "STAGE 0: Preflight checks..."
bash "${SCRIPT_DIR}/00_preflight.sh" "${BASE_DIR}"
if [ $? -ne 0 ]; then
    log "ERROR: Preflight failed. Fix issues above before continuing."
    exit 1
fi
log "Preflight passed."

# ── Stage 1: Download ─────────────────────────────────────────────────────────
if [ "$SKIP_DOWNLOAD" = false ]; then
    log "STAGE 1: Downloading samples..."
    bash "${SCRIPT_DIR}/01_download.sh" "${BASE_DIR}"
    [ $? -ne 0 ] && log "ERROR: Download stage failed." && exit 1
    log "Download complete."
else
    log "STAGE 1: Skipped (--skip-download)"
fi

# ── Stage 2: Trim ─────────────────────────────────────────────────────────────
if [ "$SKIP_TRIM" = false ]; then
    log "STAGE 2: Trimming reads..."
    bash "${SCRIPT_DIR}/02_trim.sh" "${BASE_DIR}"
    [ $? -ne 0 ] && log "ERROR: Trim stage failed." && exit 1
    log "Trimming complete."
else
    log "STAGE 2: Skipped (--skip-trim)"
fi

# ── Stage 3: Align ────────────────────────────────────────────────────────────
if [ "$SKIP_ALIGN" = false ]; then
    log "STAGE 3: Aligning reads..."
    bash "${SCRIPT_DIR}/03_align.sh" "${BASE_DIR}"
    [ $? -ne 0 ] && log "ERROR: Alignment stage failed." && exit 1
    log "Alignment complete."
else
    log "STAGE 3: Skipped (--skip-align)"
fi

# ── Stage 4: StringTie ────────────────────────────────────────────────────────
log "STAGE 4: StringTie assembly..."
bash "${SCRIPT_DIR}/04_stringtie.sh" "${BASE_DIR}"
[ $? -ne 0 ] && log "ERROR: StringTie stage failed." && exit 1
log "StringTie complete."

# ── Stage 5: Merge matrix ─────────────────────────────────────────────────────
log "STAGE 5: Merging expression matrix..."
python3 "${SCRIPT_DIR}/05_merge_matrix.py" "${BASE_DIR}"
[ $? -ne 0 ] && log "ERROR: Matrix merge failed." && exit 1
log "Matrix merge complete."

# ── Stage 6: DE analysis ──────────────────────────────────────────────────────
log "STAGE 6: Differential expression analysis..."
python3 "${SCRIPT_DIR}/06_de_analysis.py" "${BASE_DIR}"
[ $? -ne 0 ] && log "ERROR: DE analysis failed." && exit 1
log "DE analysis complete."

# ── Stage 7: Visualizations ───────────────────────────────────────────────────
log "STAGE 7: Generating visualizations..."
python3 "${SCRIPT_DIR}/07_visualize.py" "${BASE_DIR}"
[ $? -ne 0 ] && log "ERROR: Visualization failed." && exit 1
log "Visualizations complete."

# ── Final summary ─────────────────────────────────────────────────────────────
log "============================================================"
log "  PIPELINE COMPLETE"
log "  Key outputs:"
log "    ${BASE_DIR}/results/stringtie_tpm_matrix.csv"
log "    ${BASE_DIR}/results/lncrna_differential_expression.csv"
log "    ${BASE_DIR}/results/significant_stress_responsive_lncrnas.csv"
log "    ${BASE_DIR}/results/volcano_plot.png"
log "    ${BASE_DIR}/results/heatmap_significant_lncrnas.png"
log "    ${BASE_DIR}/results/boxplots_top_lncrnas.png"
log "  Full log: ${LOGFILE}"
log "============================================================"
