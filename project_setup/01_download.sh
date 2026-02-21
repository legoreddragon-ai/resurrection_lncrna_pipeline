#!/bin/bash
# =============================================================================
# 01_download.sh — Download RNA-Seq Data + Reference Files
# =============================================================================
# Usage: bash setup/01_download.sh <BASE_DIR>
# Downloads:
#   - 10 paired-end FASTQ samples from NCBI SRA
#   - Transcriptome FASTA from GEO
#   - Builds HISAT2 index
#   - Prepares GTF annotation
# =============================================================================

set -uo pipefail
BASE_DIR="${1:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$LOG_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/01_download.log"; }

THREADS=8
RAW_DIR="${BASE_DIR}/data/fastq_raw"
REF_DIR="${BASE_DIR}/references"
mkdir -p "$RAW_DIR" "$REF_DIR"

# ── Sample manifest ───────────────────────────────────────────────────────────
declare -A SAMPLES=(
    ["SRR12542161"]="cp_hyd_r1"
    ["SRR12542166"]="cp_hyd_r3"
    ["SRR12542171"]="cp_rehyd_r1"
    ["SRR12542176"]="cp_rehyd_r2"
    ["SRR12542181"]="cp_rehyd_r3"
    ["SRR12542186"]="cp_wc2_r2"
    ["SRR12542191"]="cp_wc2_r3"
    ["SRR12542196"]="cp_wc60_r1"
    ["SRR12542206"]="cp_wc60_r2"
    ["SRR12542211"]="cp_wc60_r3"
)

# ── 1. Download FASTQ samples ─────────────────────────────────────────────────
log "=== Downloading FASTQ samples ==="
FAILED=()

for ACC in "${!SAMPLES[@]}"; do
    NAME="${SAMPLES[$ACC]}"
    R1="${RAW_DIR}/${NAME}_1.fastq.gz"
    R2="${RAW_DIR}/${NAME}_2.fastq.gz"

    if [ -s "$R1" ] && [ -s "$R2" ]; then
        log "SKIP: ${NAME} already downloaded ($(du -sh "$R1" | cut -f1) + $(du -sh "$R2" | cut -f1))"
        continue
    fi

    log "Downloading ${NAME} (${ACC})..."

    # Primary: fasterq-dump
    fasterq-dump "$ACC" \
        --outdir "${RAW_DIR}/" \
        --outfile "${NAME}.fastq" \
        --split-files \
        --threads "$THREADS" \
        --progress \
        2>>"${LOG_DIR}/01_download.log"

    if [ $? -eq 0 ] && [ -f "${RAW_DIR}/${NAME}_1.fastq" ]; then
        log "Compressing ${NAME}..."
        gzip "${RAW_DIR}/${NAME}_1.fastq"
        gzip "${RAW_DIR}/${NAME}_2.fastq"
        log "Done: ${NAME}"
    else
        log "WARNING: fasterq-dump failed for ${ACC} — trying direct wget..."
        # Fallback: direct NCBI download
        SUBDIR=$(echo "$ACC" | rev | cut -c1-3 | rev)
        NCBI_URL="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/SRR012/${ACC:3:5}/${ACC}/${ACC}.lite.1"
        wget --no-check-certificate -c "$NCBI_URL" \
             -O "${RAW_DIR}/${ACC}.sra" \
             2>>"${LOG_DIR}/01_download.log"

        if [ -s "${RAW_DIR}/${ACC}.sra" ]; then
            fasterq-dump "${RAW_DIR}/${ACC}.sra" \
                --outdir "${RAW_DIR}/" \
                --outfile "${NAME}.fastq" \
                --split-files \
                2>>"${LOG_DIR}/01_download.log"
            gzip "${RAW_DIR}/${NAME}_1.fastq" 2>/dev/null
            gzip "${RAW_DIR}/${NAME}_2.fastq" 2>/dev/null
            rm -f "${RAW_DIR}/${ACC}.sra"
        else
            log "ERROR: All download methods failed for ${ACC}"
            FAILED+=("$ACC")
        fi
    fi
done

# Verify downloads
log "=== Verifying downloads ==="
ZERO_FILES=$(find "$RAW_DIR" -name '*.fastq.gz' -size 0 2>/dev/null)
if [ -n "$ZERO_FILES" ]; then
    log "ERROR: 0-byte files found — delete and re-run:"
    echo "$ZERO_FILES"
    exit 1
fi
COUNT=$(ls "${RAW_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
log "FASTQ files present: ${COUNT} (expected: 20)"
[ "$COUNT" -lt 20 ] && log "WARNING: Expected 20 files, found ${COUNT}"

if [ ${#FAILED[@]} -gt 0 ]; then
    log "Failed accessions: ${FAILED[*]}"
    log "Re-run this script after resolving network issues — already downloaded samples will be skipped."
fi

# ── 2. Download reference FASTA ───────────────────────────────────────────────
log "=== Downloading reference files ==="
FA="${REF_DIR}/GSE157098_Cp_transcriptome_assembly_V2.fa"
if [ -s "$FA" ]; then
    log "SKIP: Transcriptome FASTA already present"
else
    log "Downloading transcriptome FASTA..."
    wget --no-check-certificate -c \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_transcriptome_assembly_V2.fa.gz" \
        -O "${FA}.gz" \
        2>>"${LOG_DIR}/01_download.log"
    gunzip "${FA}.gz"

    CONTIGS=$(grep -c '>' "$FA")
    log "Transcriptome FASTA downloaded: ${CONTIGS} contigs"
    [ "$CONTIGS" -lt 40000 ] && log "WARNING: Expected ~48045 contigs, got ${CONTIGS}"
fi

# ── 3. Download annotation ────────────────────────────────────────────────────
ANNOT="${REF_DIR}/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt"
if [ -s "$ANNOT" ]; then
    log "SKIP: Annotation already present"
else
    log "Downloading annotation..."
    wget --no-check-certificate -c \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz" \
        -O "${ANNOT}.gz" \
        2>>"${LOG_DIR}/01_download.log"
    gunzip "${ANNOT}.gz"
    log "Annotation downloaded"
fi

# ── 4. Build HISAT2 index ─────────────────────────────────────────────────────
INDEX="${REF_DIR}/Cp_transcriptome_index"
if ls "${INDEX}.1.ht2" &>/dev/null; then
    log "SKIP: HISAT2 index already built"
else
    log "Building HISAT2 index (this takes 10-15 minutes)..."
    hisat2-build "$FA" "$INDEX" -p "$THREADS" \
        2>>"${LOG_DIR}/01_download.log"
    COUNT=$(ls "${REF_DIR}"/Cp_transcriptome_index*.ht2 2>/dev/null | wc -l)
    log "HISAT2 index built: ${COUNT} files"
    [ "$COUNT" -ne 8 ] && log "WARNING: Expected 8 index files, got ${COUNT}"
fi

# ── 5. Prepare GTF annotation ─────────────────────────────────────────────────
GTF="${REF_DIR}/annotation.gtf"
GFF="${REF_DIR}/annotation.gff3"
if [ -s "$GTF" ]; then
    log "SKIP: annotation.gtf already present"
elif [ -s "$GFF" ]; then
    log "Converting GFF3 to GTF..."
    gffread "$GFF" -T -o "$GTF" 2>>"${LOG_DIR}/01_download.log"
    LINES=$(wc -l < "$GTF")
    log "annotation.gtf created: ${LINES} lines"
else
    log "WARNING: Neither annotation.gtf nor annotation.gff3 found in ${REF_DIR}"
    log "Please place annotation.gff3 in ${REF_DIR} and rerun, or provide annotation.gtf directly."
fi

log "=== Stage 1 complete ==="
