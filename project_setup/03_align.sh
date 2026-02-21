#!/bin/bash
# =============================================================================
# 03_align.sh — HISAT2 Alignment + SAMtools Sort & Index
# =============================================================================
# Usage: bash setup/03_align.sh <BASE_DIR>
# Aligns all 10 trimmed samples to transcriptome index.
# Skips samples whose sorted BAM + index already exist.
# Output: results/hisat2_aligned/<sample>.sorted.bam
#                                 <sample>.sorted.bam.bai
# =============================================================================

set -uo pipefail
BASE_DIR="${1:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/03_align.log"; }

THREADS=8
TRIM_DIR="${BASE_DIR}/data/fastq_trimmed"
BAM_DIR="${BASE_DIR}/results/hisat2_aligned"
REF_DIR="${BASE_DIR}/references"
INDEX="${REF_DIR}/Cp_transcriptome_index"
mkdir -p "$BAM_DIR"

# Verify index exists
if ! ls "${INDEX}.1.ht2" &>/dev/null; then
    log "ERROR: HISAT2 index not found at ${INDEX}"
    log "Run 01_download.sh first to build the index."
    exit 1
fi

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

log "=== Starting alignment for ${#SAMPLES[@]} samples ==="

for NAME in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${NAME}.sorted.bam"
    BAI="${BAM}.bai"

    if [ -s "$BAM" ] && [ -s "$BAI" ]; then
        log "SKIP: ${NAME} already aligned ($(du -sh "$BAM" | cut -f1))"
        continue
    fi

    R1="${TRIM_DIR}/${NAME}_1_paired.fastq.gz"
    R2="${TRIM_DIR}/${NAME}_2_paired.fastq.gz"

    if [ ! -s "$R1" ] || [ ! -s "$R2" ]; then
        log "ERROR: Trimmed input missing for ${NAME} — run 02_trim.sh first"
        continue
    fi

    log "Aligning ${NAME}..."
    hisat2 \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        --dta \
        --rna-strandness RF \
        2>"${LOG_DIR}/03_align_${NAME}.log" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$BAM"

    if [ $? -ne 0 ] || [ ! -s "$BAM" ]; then
        log "ERROR: Alignment failed for ${NAME}"
        continue
    fi

    samtools index "$BAM"

    # Report alignment rate
    RATE=$(grep 'overall alignment rate' "${LOG_DIR}/03_align_${NAME}.log" | tail -1)
    log "Done: ${NAME} — ${RATE}"
done

# Summary
log "=== Alignment summary ==="
for NAME in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${NAME}.sorted.bam"
    if [ -s "$BAM" ]; then
        RATE=$(grep 'overall alignment rate' "${LOG_DIR}/03_align_${NAME}.log" 2>/dev/null | tail -1 | awk '{print $1}')
        log "  ${NAME}: $(du -sh "$BAM" | cut -f1) | ${RATE} alignment"
    else
        log "  ${NAME}: MISSING"
    fi
done

log "=== Stage 3 complete ==="
