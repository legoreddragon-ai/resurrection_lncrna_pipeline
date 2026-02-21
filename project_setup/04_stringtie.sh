#!/bin/bash
# =============================================================================
# 04_stringtie.sh — StringTie Transcript Quantification
# =============================================================================
# Usage: bash setup/04_stringtie.sh <BASE_DIR>
# Quantifies transcript expression for all 10 samples.
# CRITICAL: Uses -G (annotation) and -e (expression mode) flags.
#           Without -G, StringTie assembles novel transcripts instead of
#           quantifying known ones — producing ~6 lines not ~48045.
# Output: results/stringtie_output/<sample>.gtf
#                                   <sample>_abundance.txt
# =============================================================================

set -uo pipefail
BASE_DIR="${1:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/04_stringtie.log"; }

THREADS=8
BAM_DIR="${BASE_DIR}/results/hisat2_aligned"
ST_DIR="${BASE_DIR}/results/stringtie_output"
GTF="${BASE_DIR}/references/annotation.gtf"
mkdir -p "$ST_DIR"

# Verify annotation exists
if [ ! -s "$GTF" ]; then
    log "ERROR: annotation.gtf not found at ${GTF}"
    log "Run 01_download.sh to prepare reference files."
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

log "=== Starting StringTie for ${#SAMPLES[@]} samples ==="

for NAME in "${SAMPLES[@]}"; do
    BAM="${BAM_DIR}/${NAME}.sorted.bam"
    GTF_OUT="${ST_DIR}/${NAME}.gtf"
    ABUND="${ST_DIR}/${NAME}_abundance.txt"

    # Skip only if abundance file has correct number of lines
    if [ -s "$ABUND" ]; then
        LINES=$(wc -l < "$ABUND")
        if [ "$LINES" -gt 1000 ]; then
            log "SKIP: ${NAME} already assembled (${LINES} lines in abundance)"
            continue
        else
            log "WARNING: ${NAME} abundance file has only ${LINES} lines — likely missing -G flag. Re-running."
            rm -f "$GTF_OUT" "$ABUND"
        fi
    fi

    if [ ! -s "$BAM" ]; then
        log "ERROR: BAM file missing for ${NAME} — run 03_align.sh first"
        continue
    fi

    log "Assembling ${NAME}..."
    stringtie "$BAM" \
        -G "$GTF" \
        -o "$GTF_OUT" \
        -A "$ABUND" \
        -p "$THREADS" \
        -e \
        2>"${LOG_DIR}/04_stringtie_${NAME}.log"

    if [ $? -ne 0 ] || [ ! -s "$ABUND" ]; then
        log "ERROR: StringTie failed for ${NAME}"
        continue
    fi

    LINES=$(wc -l < "$ABUND")
    if [ "$LINES" -lt 1000 ]; then
        log "ERROR: ${NAME} abundance file has only ${LINES} lines — check annotation GTF"
    else
        log "Done: ${NAME} (${LINES} genes quantified)"
    fi
done

# Summary
log "=== StringTie summary ==="
for NAME in "${SAMPLES[@]}"; do
    ABUND="${ST_DIR}/${NAME}_abundance.txt"
    if [ -s "$ABUND" ]; then
        LINES=$(wc -l < "$ABUND")
        log "  ${NAME}: ${LINES} lines"
    else
        log "  ${NAME}: MISSING"
    fi
done

log "=== Stage 4 complete ==="
