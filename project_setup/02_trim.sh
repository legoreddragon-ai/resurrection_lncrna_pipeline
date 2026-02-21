#!/bin/bash
# =============================================================================
# 02_trim.sh — Quality Trimming with Trimmomatic
# =============================================================================
# Usage: bash setup/02_trim.sh <BASE_DIR>
# Trims all 10 samples. Skips samples already trimmed.
# Output: data/fastq_trimmed/<sample>_1_paired.fastq.gz
#                             <sample>_2_paired.fastq.gz
# =============================================================================

set -uo pipefail
BASE_DIR="${1:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/02_trim.log"; }

THREADS=8
RAW_DIR="${BASE_DIR}/data/fastq_raw"
TRIM_DIR="${BASE_DIR}/data/fastq_trimmed"
mkdir -p "$TRIM_DIR"

# Find Trimmomatic adapter file
ADAPTERS=$(find ~/miniconda3 /usr -name "TruSeq3-PE.fa" 2>/dev/null | head -1)
if [ -z "$ADAPTERS" ]; then
    ADAPTERS=$(find / -name "TruSeq3-PE.fa" 2>/dev/null | head -1)
fi
if [ -z "$ADAPTERS" ]; then
    log "ERROR: TruSeq3-PE.fa adapter file not found. Specify path manually."
    log "Expected location: ~/miniconda3/share/trimmomatic/adapters/TruSeq3-PE.fa"
    exit 1
fi
log "Using adapters: ${ADAPTERS}"

SAMPLES=(
    "cp_hyd_r1"
    "cp_hyd_r3"
    "cp_rehyd_r1"
    "cp_rehyd_r2"
    "cp_rehyd_r3"
    "cp_wc2_r2"
    "cp_wc2_r3"
    "cp_wc60_r1"
    "cp_wc60_r2"
    "cp_wc60_r3"
)

log "=== Starting trimming for ${#SAMPLES[@]} samples ==="

for NAME in "${SAMPLES[@]}"; do
    R1_IN="${RAW_DIR}/${NAME}_1.fastq.gz"
    R2_IN="${RAW_DIR}/${NAME}_2.fastq.gz"
    R1_OUT="${TRIM_DIR}/${NAME}_1_paired.fastq.gz"
    R2_OUT="${TRIM_DIR}/${NAME}_2_paired.fastq.gz"
    R1_UNPAIRED="${TRIM_DIR}/${NAME}_1_unpaired.fastq.gz"
    R2_UNPAIRED="${TRIM_DIR}/${NAME}_2_unpaired.fastq.gz"

    # Skip if output already exists and non-empty
    if [ -s "$R1_OUT" ] && [ -s "$R2_OUT" ]; then
        log "SKIP: ${NAME} already trimmed"
        continue
    fi

    # Check input exists
    if [ ! -s "$R1_IN" ] || [ ! -s "$R2_IN" ]; then
        log "ERROR: Input files missing or empty for ${NAME}: ${R1_IN}"
        continue
    fi

    log "Trimming ${NAME}..."
    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$R1_IN" "$R2_IN" \
        "$R1_OUT" "$R1_UNPAIRED" \
        "$R2_OUT" "$R2_UNPAIRED" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2>>"${LOG_DIR}/02_trim_${NAME}.log"
    # NOTE: GCLocker warnings from Trimmomatic are harmless Java GC messages.
    # They do NOT indicate failure. We intentionally do not use set -e here.

    if [ -s "$R1_OUT" ] && [ -s "$R2_OUT" ]; then
        R1_SIZE=$(du -sh "$R1_OUT" | cut -f1)
        R2_SIZE=$(du -sh "$R2_OUT" | cut -f1)
        log "Done: ${NAME} (${R1_SIZE} + ${R2_SIZE})"
    else
        log "ERROR: Trimmomatic produced empty output for ${NAME}"
    fi
done

# Summary
log "=== Trim summary ==="
DONE=$(ls "${TRIM_DIR}"/*_paired.fastq.gz 2>/dev/null | wc -l)
log "Paired trimmed files: ${DONE} (expected: 20)"
[ "$DONE" -lt 20 ] && log "WARNING: Expected 20 paired files"

log "=== Stage 2 complete ==="
