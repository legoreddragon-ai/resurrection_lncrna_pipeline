#!/bin/bash
# =============================================================================
# 00_preflight.sh — Environment & Prerequisites Check
# =============================================================================
# Usage: bash setup/00_preflight.sh <BASE_DIR>
# Called automatically by run_pipeline.sh
# Can also be run standalone: bash setup/00_preflight.sh
# =============================================================================

BASE_DIR="${1:-$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")}"
PASS=0
FAIL=0
WARN=0

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok()   { echo -e "  ${GREEN}[PASS]${NC} $*"; ((PASS++)); }
fail() { echo -e "  ${RED}[FAIL]${NC} $*"; ((FAIL++)); }
warn() { echo -e "  ${YELLOW}[WARN]${NC} $*"; ((WARN++)); }

echo "============================================================"
echo "  PREFLIGHT CHECK"
echo "  Base directory: ${BASE_DIR}"
echo "============================================================"

# ── OS ────────────────────────────────────────────────────────────────────────
echo ""
echo "[ System ]"
OS=$(uname -s)
ARCH=$(uname -m)
if [ "$OS" = "Linux" ] && [ "$ARCH" = "x86_64" ]; then
    ok "OS: Linux x86_64"
else
    warn "OS: ${OS} ${ARCH} — pipeline designed for Linux x86_64"
fi

# RAM
RAM_GB=$(free -g | awk '/^Mem:/{print $2}')
if [ "$RAM_GB" -ge 8 ]; then
    ok "RAM: ${RAM_GB} GB"
else
    warn "RAM: ${RAM_GB} GB — recommend 8+ GB"
fi

# Disk
DISK_GB=$(df -BG "${BASE_DIR}" 2>/dev/null | awk 'NR==2{gsub("G",""); print $4}')
if [ -n "$DISK_GB" ] && [ "$DISK_GB" -ge 80 ]; then
    ok "Disk available: ${DISK_GB} GB"
else
    warn "Disk available: ${DISK_GB:-unknown} GB — recommend 80+ GB"
fi

# ── Tools ─────────────────────────────────────────────────────────────────────
echo ""
echo "[ Required Tools ]"

check_tool() {
    local tool=$1
    local version_flag=${2:---version}
    if command -v "$tool" &>/dev/null; then
        local ver
        ver=$("$tool" $version_flag 2>&1 | head -1)
        ok "${tool}: ${ver}"
    else
        fail "${tool}: NOT FOUND — install before continuing"
    fi
}

check_tool python3 "--version"
check_tool pip3 "--version"
check_tool git "--version"
check_tool curl "--version"
check_tool wget "--version"
check_tool gzip "--version"

echo ""
echo "[ Bioinformatics Tools ]"
check_tool trimmomatic "-version"
check_tool hisat2 "--version"
check_tool samtools "--version"
check_tool stringtie "--version"

# fasterq-dump version flag is different
if command -v fasterq-dump &>/dev/null; then
    VER=$(fasterq-dump --version 2>&1 | head -1)
    ok "fasterq-dump: ${VER}"
    # Check for old version
    if echo "$VER" | grep -q "2\.9"; then
        warn "fasterq-dump 2.9.x detected — may have TLS issues. Upgrade to 3.x recommended."
    fi
else
    fail "fasterq-dump: NOT FOUND — install SRA Toolkit 3.x"
fi

# ── Python packages ───────────────────────────────────────────────────────────
echo ""
echo "[ Python Packages ]"
python3 -c "import pandas; print('pandas:', pandas.__version__)" 2>/dev/null && ok "pandas" || fail "pandas not installed"
python3 -c "import numpy; print('numpy:', numpy.__version__)" 2>/dev/null && ok "numpy" || fail "numpy not installed"
python3 -c "import scipy; print('scipy:', scipy.__version__)" 2>/dev/null && ok "scipy" || fail "scipy not installed"
python3 -c "import matplotlib; print('matplotlib:', matplotlib.__version__)" 2>/dev/null && ok "matplotlib" || fail "matplotlib not installed"

# ── NCBI connectivity ─────────────────────────────────────────────────────────
echo ""
echo "[ Network / NCBI Connectivity ]"
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" --max-time 10 https://www.ncbi.nlm.nih.gov 2>/dev/null)
if [ "$HTTP_CODE" = "200" ]; then
    ok "NCBI HTTPS reachable (HTTP ${HTTP_CODE})"
else
    warn "NCBI HTTPS returned ${HTTP_CODE} — downloads may fail"
fi

# ── Project directories ───────────────────────────────────────────────────────
echo ""
echo "[ Project Directories ]"
DIRS=(
    "data/fastq_raw"
    "data/fastq_trimmed"
    "results/hisat2_aligned"
    "results/stringtie_output"
    "results/figures"
    "references"
    "logs"
    "scripts"
)
for d in "${DIRS[@]}"; do
    FULL="${BASE_DIR}/${d}"
    if [ -d "$FULL" ]; then
        ok "${d}/"
    else
        warn "${d}/ missing — creating..."
        mkdir -p "$FULL"
    fi
done

# ── Reference files ───────────────────────────────────────────────────────────
echo ""
echo "[ Reference Files ]"
REF="${BASE_DIR}/references"

if [ -f "${REF}/GSE157098_Cp_transcriptome_assembly_V2.fa" ]; then
    CONTIGS=$(grep -c '>' "${REF}/GSE157098_Cp_transcriptome_assembly_V2.fa" 2>/dev/null || echo "?")
    ok "Transcriptome FASTA (${CONTIGS} contigs)"
else
    warn "Transcriptome FASTA not found — run Stage 1 or download manually"
fi

if ls "${REF}/Cp_transcriptome_index.1.ht2" &>/dev/null; then
    COUNT=$(ls "${REF}"/Cp_transcriptome_index*.ht2 2>/dev/null | wc -l)
    ok "HISAT2 index (${COUNT} files)"
else
    warn "HISAT2 index not built — run Stage 1"
fi

if [ -f "${REF}/annotation.gtf" ]; then
    LINES=$(wc -l < "${REF}/annotation.gtf")
    ok "annotation.gtf (${LINES} lines)"
else
    warn "annotation.gtf not found — run Stage 1"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo "  PREFLIGHT SUMMARY: ${PASS} passed | ${WARN} warnings | ${FAIL} failed"
echo "============================================================"

if [ "$FAIL" -gt 0 ]; then
    echo -e "  ${RED}Fix all FAIL items before running the pipeline.${NC}"
    exit 1
elif [ "$WARN" -gt 0 ]; then
    echo -e "  ${YELLOW}Warnings present — review above. Pipeline may still run.${NC}"
    exit 0
else
    echo -e "  ${GREEN}All checks passed. Ready to run.${NC}"
    exit 0
fi
