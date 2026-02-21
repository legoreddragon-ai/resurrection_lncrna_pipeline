#!/bin/bash
# =============================================================================
# Resurrection lncRNA Pipeline - Remaining 9 Samples
# Tasks: Trim → Align → Assemble → Merge → DE Analysis → Visualize
# =============================================================================

set -uo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

THREADS=8
BASE_DIR="/mnt/neil/resurrection_lncrna_pipeline"
DATA_RAW="${BASE_DIR}/data/fastq_raw"
DATA_TRIMMED="${BASE_DIR}/data/fastq_trimmed"
REF_DIR="${BASE_DIR}/references"
RESULTS_DIR="${BASE_DIR}/results"
HISAT2_DIR="${RESULTS_DIR}/hisat2_aligned"
STRINGTIE_DIR="${RESULTS_DIR}/stringtie_output"
LOG_DIR="${BASE_DIR}/logs"

GENOME_INDEX="${REF_DIR}/Cp_transcriptome_index"
ANNOTATION_GTF="${REF_DIR}/annotation.gtf"

# Trimmomatic adapter file (adjust path if needed)
ADAPTERS="$(dirname $(which trimmomatic))/../share/trimmomatic/adapters/TruSeq3-PE.fa"

# Samples to process (the 9 remaining)
declare -A SAMPLES=(
    ["cp_hyd_r3"]="cp_hyd_r3"
    ["cp_rehyd_r1"]="cp_rehyd_r1"
    ["cp_rehyd_r2"]="cp_rehyd_r2"
    ["cp_rehyd_r3"]="cp_rehyd_r3"
    ["cp_wc2_r2"]="cp_wc2_r2"
    ["cp_wc2_r3"]="cp_wc2_r3"
    ["cp_wc60_r1"]="cp_wc60_r1"
    ["cp_wc60_r2"]="cp_wc60_r2"
    ["cp_wc60_r3"]="cp_wc60_r3"
)

# All 10 samples for matrix merge (includes cp_hyd_r1 already done)
ALL_SAMPLES=(
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

# =============================================================================
# HELPERS
# =============================================================================

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_DIR}/pipeline.log"; }
check_file() { [ -s "$1" ] || { log "ERROR: Missing or empty file: $1"; exit 1; }; }

mkdir -p "$DATA_TRIMMED" "$HISAT2_DIR" "$STRINGTIE_DIR" "$LOG_DIR" "$RESULTS_DIR"

# =============================================================================
# TASK 1: TRIMMOMATIC - Trim remaining 9 samples
# =============================================================================

log "====== TASK 1: TRIMMING ======"

for NAME in "${!SAMPLES[@]}"; do
    R1="${DATA_RAW}/${NAME}_1.fastq.gz"
    R2="${DATA_RAW}/${NAME}_2.fastq.gz"
    OUT1="${DATA_TRIMMED}/${NAME}_1_paired.fastq.gz"
    OUT2="${DATA_TRIMMED}/${NAME}_2_paired.fastq.gz"
    UNPAIRED1="${DATA_TRIMMED}/${NAME}_1_unpaired.fastq.gz"
    UNPAIRED2="${DATA_TRIMMED}/${NAME}_2_unpaired.fastq.gz"

    # Skip if already done
    if [ -s "$OUT1" ] && [ -s "$OUT2" ]; then
        log "SKIP: $NAME already trimmed"
        continue
    fi

    check_file "$R1"
    check_file "$R2"

    log "Trimming $NAME..."
    trimmomatic PE \
        -threads "$THREADS" \
        -phred33 \
        "$R1" "$R2" \
        "$OUT1" "$UNPAIRED1" \
        "$OUT2" "$UNPAIRED2" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10:2:keepBothReads \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2>> "${LOG_DIR}/${NAME}_trimmomatic.log"

    check_file "$OUT1"
    log "Done trimming: $NAME"
done

log "====== TASK 1 COMPLETE: All samples trimmed ======"

# =============================================================================
# TASK 2: HISAT2 ALIGNMENT
# =============================================================================

log "====== TASK 2: ALIGNMENT ======"

# Verify index exists
if [ ! -f "${GENOME_INDEX}.1.ht2" ]; then
    log "ERROR: HISAT2 index not found at ${GENOME_INDEX}"
    log "Build it with: hisat2-build ${REF_DIR}/genome.fa ${GENOME_INDEX}"
    exit 1
fi

for NAME in "${!SAMPLES[@]}"; do
    BAM="${HISAT2_DIR}/${NAME}.sorted.bam"

    if [ -s "$BAM" ] && [ -s "${BAM}.bai" ]; then
        log "SKIP: $NAME already aligned"
        continue
    fi

    R1="${DATA_TRIMMED}/${NAME}_1_paired.fastq.gz"
    R2="${DATA_TRIMMED}/${NAME}_2_paired.fastq.gz"
    check_file "$R1"
    check_file "$R2"

    log "Aligning $NAME..."
    hisat2 \
        -x "$GENOME_INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        --dta \
        --rna-strandness RF \
        2>> "${LOG_DIR}/${NAME}_hisat2.log" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$BAM"

    samtools index "$BAM"
    check_file "$BAM"
    log "Done aligning: $NAME"
done

log "====== TASK 2 COMPLETE: All samples aligned ======"

# =============================================================================
# TASK 3: STRINGTIE TRANSCRIPT ASSEMBLY
# =============================================================================

log "====== TASK 3: STRINGTIE ASSEMBLY ======"

check_file "$ANNOTATION_GTF"

for NAME in "${!SAMPLES[@]}"; do
    BAM="${HISAT2_DIR}/${NAME}.sorted.bam"
    GTF_OUT="${STRINGTIE_DIR}/${NAME}.gtf"
    ABUND_OUT="${STRINGTIE_DIR}/${NAME}_abundance.txt"

    if [ -s "$GTF_OUT" ] && [ -s "$ABUND_OUT" ]; then
        log "SKIP: $NAME already assembled"
        continue
    fi

    check_file "$BAM"

    log "Assembling $NAME..."
    stringtie "$BAM" \
        -G "$ANNOTATION_GTF" \
        -o "$GTF_OUT" \
        -A "$ABUND_OUT" \
        -p "$THREADS" \
        -e \
        2>> "${LOG_DIR}/${NAME}_stringtie.log"

    check_file "$GTF_OUT"
    check_file "$ABUND_OUT"
    log "Done assembling: $NAME"
done

log "====== TASK 3 COMPLETE: All samples assembled ======"

# =============================================================================
# TASK 4: MERGE EXPRESSION MATRICES
# =============================================================================

log "====== TASK 4: MERGING EXPRESSION MATRICES ======"

TPM_MATRIX="${RESULTS_DIR}/stringtie_tpm_matrix.csv"

python3 << PYEOF
import pandas as pd
import os

stringtie_dir = "${STRINGTIE_DIR}"
results_dir = "${RESULTS_DIR}"
all_samples = ${ALL_SAMPLES[@]/#/[\"} 
# Rebuild as proper python list
all_samples = """${ALL_SAMPLES[*]}""".split()

dfs = []
missing = []

for sample in all_samples:
    abund_file = os.path.join(stringtie_dir, f"{sample}_abundance.txt")
    if not os.path.exists(abund_file) or os.path.getsize(abund_file) == 0:
        print(f"WARNING: Missing abundance file for {sample}: {abund_file}")
        missing.append(sample)
        continue
    df = pd.read_csv(abund_file, sep='\t')
    df = df[['Gene ID', 'Gene Name', 'TPM']].copy()
    df = df.rename(columns={'TPM': sample})
    dfs.append(df)
    print(f"Loaded: {sample} ({len(df)} genes)")

if not dfs:
    print("ERROR: No abundance files found!")
    exit(1)

if missing:
    print(f"WARNING: Missing samples: {missing}")

# Merge all on Gene ID
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df.drop(columns=['Gene Name']), on='Gene ID', how='outer')

merged = merged.fillna(0)
out_path = os.path.join(results_dir, "stringtie_tpm_matrix.csv")
merged.to_csv(out_path, index=False)
print(f"Saved TPM matrix: {out_path}")
print(f"Shape: {merged.shape[0]} genes x {merged.shape[1]} columns")
PYEOF

log "====== TASK 4 COMPLETE: TPM matrix merged ======"

# =============================================================================
# TASK 5: DIFFERENTIAL EXPRESSION ANALYSIS
# =============================================================================

log "====== TASK 5: DIFFERENTIAL EXPRESSION ANALYSIS ======"

python3 << PYEOF
import pandas as pd
import numpy as np
from scipy import stats

# Load lncRNA DE results and TPM matrix
results_dir = "${RESULTS_DIR}"
lncrna_de_path = f"{results_dir}/lncrna_differential_expression.csv"
tpm_path = f"{results_dir}/stringtie_tpm_matrix.csv"

print("Loading TPM matrix...")
tpm = pd.read_csv(tpm_path)
print(f"TPM matrix shape: {tpm.shape}")

# Define sample groups
hydrated_cols    = [c for c in tpm.columns if 'hyd_r' in c and 'rehyd' not in c and 'wc' not in c]
dehydrated_cols  = [c for c in tpm.columns if 'wc' in c]
rehydrated_cols  = [c for c in tpm.columns if 'rehyd' in c]

print(f"Hydrated samples:   {hydrated_cols}")
print(f"Dehydrated samples: {dehydrated_cols}")
print(f"Rehydrated samples: {rehydrated_cols}")

# Load coding potential to filter for lncRNAs
coding_path = f"{results_dir}/coding_potential.txt"
try:
    coding = pd.read_csv(coding_path, sep='\t', header=None, names=['Gene ID', 'coding_score'])
    lncrna_ids = set(coding[coding['coding_score'] < 0.5]['Gene ID'])
    print(f"lncRNA candidates: {len(lncrna_ids)}")
    tpm_lnc = tpm[tpm['Gene ID'].isin(lncrna_ids)].copy()
except Exception as e:
    print(f"Could not load coding potential ({e}), using all genes")
    tpm_lnc = tpm.copy()

# Require at least 1 TPM in hydrated
if hydrated_cols:
    tpm_lnc = tpm_lnc[tpm_lnc[hydrated_cols].mean(axis=1) >= 1]

print(f"After TPM filter: {len(tpm_lnc)} lncRNAs")

# Differential expression: hydrated vs dehydrated
results = []
for _, row in tpm_lnc.iterrows():
    hyd_vals  = row[hydrated_cols].values.astype(float)   if hydrated_cols   else np.array([])
    dehy_vals = row[dehydrated_cols].values.astype(float) if dehydrated_cols else np.array([])

    if len(hyd_vals) < 2 or len(dehy_vals) < 2:
        continue

    mean_hyd  = np.mean(hyd_vals)
    mean_dehy = np.mean(dehy_vals)

    # Avoid log2(0)
    log2fc = np.log2((mean_dehy + 0.1) / (mean_hyd + 0.1))

    _, pval = stats.ttest_ind(hyd_vals, dehy_vals)

    results.append({
        'gene_id':       row['Gene ID'],
        'gene_name':     row.get('Gene Name', ''),
        'mean_hydrated': round(mean_hyd, 3),
        'mean_dehydrated': round(mean_dehy, 3),
        'log2FC':        round(log2fc, 4),
        'pvalue':        pval,
        **{s: round(row[s], 3) for s in hydrated_cols + dehydrated_cols + rehydrated_cols if s in row}
    })

de_df = pd.DataFrame(results)

# BH FDR correction
from scipy.stats import rankdata
pvals = de_df['pvalue'].values
n = len(pvals)
ranked = rankdata(pvals)
fdr = np.minimum(1, pvals * n / ranked)
# Ensure monotonicity
for i in range(n - 2, -1, -1):
    fdr[i] = min(fdr[i], fdr[i + 1])
de_df['padj'] = fdr

# Save full results
full_out = f"{results_dir}/lncrna_differential_expression.csv"
de_df.to_csv(full_out, index=False)
print(f"Saved full DE results: {full_out} ({len(de_df)} lncRNAs)")

# Filter significant
sig = de_df[(de_df['padj'] < 0.05) & (de_df['log2FC'].abs() > 1)].copy()
sig = sig.sort_values('log2FC')

sig_out = f"{results_dir}/significant_stress_responsive_lncrnas.csv"
sig.to_csv(sig_out, index=False)
print(f"\nSignificant stress-responsive lncRNAs: {len(sig)}")
print(sig[['gene_id', 'gene_name', 'mean_hydrated', 'mean_dehydrated', 'log2FC', 'padj']].to_string())
PYEOF

log "====== TASK 5 COMPLETE: Differential expression done ======"

# =============================================================================
# TASK 6: VISUALIZATIONS
# =============================================================================

log "====== TASK 6: VISUALIZATIONS ======"

python3 << PYEOF
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster import hierarchy

results_dir = "${RESULTS_DIR}"

de_df  = pd.read_csv(f"{results_dir}/lncrna_differential_expression.csv")
sig_df = pd.read_csv(f"{results_dir}/significant_stress_responsive_lncrnas.csv")
tpm_df = pd.read_csv(f"{results_dir}/stringtie_tpm_matrix.csv")

# ── 1. VOLCANO PLOT ──────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 7))

de_df['neg_log10_padj'] = -np.log10(de_df['padj'].clip(lower=1e-300))
nonsig = de_df[(de_df['padj'] >= 0.05) | (de_df['log2FC'].abs() <= 1)]
sig_up = de_df[(de_df['padj'] < 0.05) & (de_df['log2FC'] > 1)]
sig_dn = de_df[(de_df['padj'] < 0.05) & (de_df['log2FC'] < -1)]

ax.scatter(nonsig['log2FC'], nonsig['neg_log10_padj'], s=6, color='#cccccc', alpha=0.5, label='Not significant')
ax.scatter(sig_up['log2FC'], sig_up['neg_log10_padj'], s=20, color='#e74c3c', alpha=0.8, label=f'Up in dehydrated (n={len(sig_up)})')
ax.scatter(sig_dn['log2FC'], sig_dn['neg_log10_padj'], s=20, color='#3498db', alpha=0.8, label=f'Down in dehydrated (n={len(sig_dn)})')

# Label top candidates
top = de_df.nsmallest(14, 'padj')
for _, row in top.iterrows():
    label = row['gene_name'] if str(row.get('gene_name', '')) not in ['', 'nan'] else row['gene_id']
    ax.annotate(label, (row['log2FC'], row['neg_log10_padj']),
                fontsize=6, ha='left', va='bottom',
                xytext=(3, 3), textcoords='offset points')

ax.axvline(x=1,  color='gray', linestyle='--', linewidth=0.8)
ax.axvline(x=-1, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(y=-np.log10(0.05), color='gray', linestyle='--', linewidth=0.8)
ax.set_xlabel('log₂ Fold Change (Dehydrated / Hydrated)', fontsize=12)
ax.set_ylabel('-log₁₀ adjusted p-value', fontsize=12)
ax.set_title('Volcano Plot: lncRNA Differential Expression\n(Hydrated vs. Dehydrated)', fontsize=13)
ax.legend(frameon=True, fontsize=10)
plt.tight_layout()
plt.savefig(f"{results_dir}/volcano_plot.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: volcano_plot.png")

# ── 2. HEATMAP ───────────────────────────────────────────────────────────────
sample_cols = [c for c in tpm_df.columns if c.startswith('cp_')]
if len(sig_df) > 0 and sample_cols:
    heat_data = tpm_df[tpm_df['Gene ID'].isin(sig_df['gene_id'])].copy()
    heat_data = heat_data.set_index('Gene ID')[sample_cols]

    # Log-transform
    heat_log = np.log2(heat_data + 1)

    # Z-score per gene
    heat_z = heat_log.subtract(heat_log.mean(axis=1), axis=0).divide(
        heat_log.std(axis=1).replace(0, 1), axis=0)

    # Sort columns by condition
    col_order = sorted(sample_cols, key=lambda x: (
        0 if 'hyd' in x and 'rehyd' not in x and 'wc' not in x else
        1 if 'wc' in x else 2))
    heat_z = heat_z[col_order]

    fig, ax = plt.subplots(figsize=(12, max(6, len(heat_z) * 0.5)))
    cmap = plt.cm.RdBu_r
    im = ax.imshow(heat_z.values, aspect='auto', cmap=cmap, vmin=-2, vmax=2)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(len(heat_z)))
    ax.set_yticklabels(heat_z.index, fontsize=8)
    plt.colorbar(im, ax=ax, label='Z-score (log₂ TPM)')
    ax.set_title(f'Expression Heatmap: {len(sig_df)} Significant Stress-Responsive lncRNAs', fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{results_dir}/heatmap_significant_lncrnas.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: heatmap_significant_lncrnas.png")

# ── 3. BOXPLOTS BY CONDITION ─────────────────────────────────────────────────
if len(sig_df) > 0 and sample_cols:
    hyd_cols   = [c for c in sample_cols if 'hyd' in c and 'rehyd' not in c and 'wc' not in c]
    dehy_cols  = [c for c in sample_cols if 'wc' in c]
    rehyd_cols = [c for c in sample_cols if 'rehyd' in c]

    top14 = sig_df.head(14)
    n = len(top14)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 3.5))
    axes = axes.flatten()

    for i, (_, row) in enumerate(top14.iterrows()):
        gene_tpm = tpm_df[tpm_df['Gene ID'] == row['gene_id']]
        if gene_tpm.empty:
            continue
        gene_tpm = gene_tpm.iloc[0]

        vals_hyd   = gene_tpm[hyd_cols].values.astype(float)   if hyd_cols   else np.array([0])
        vals_dehy  = gene_tpm[dehy_cols].values.astype(float)  if dehy_cols  else np.array([0])
        vals_rehyd = gene_tpm[rehyd_cols].values.astype(float) if rehyd_cols else np.array([0])

        ax = axes[i]
        bp = ax.boxplot([vals_hyd, vals_dehy, vals_rehyd],
                        labels=['Hydrated', 'Dehydrated', 'Rehydrated'],
                        patch_artist=True,
                        medianprops=dict(color='black', linewidth=2))

        colors = ['#2ecc71', '#e74c3c', '#3498db']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        label = str(row.get('gene_name', ''))
        if label in ['', 'nan']:
            label = row['gene_id']
        ax.set_title(f"{label}\nlog₂FC={row['log2FC']:.2f}, padj={row['padj']:.3f}", fontsize=8)
        ax.set_ylabel('TPM', fontsize=8)
        ax.tick_params(labelsize=7)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle('Expression by Condition: Top Stress-Responsive lncRNAs', fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(f"{results_dir}/boxplots_top_lncrnas.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: boxplots_top_lncrnas.png")

print("\nAll visualizations complete.")
PYEOF

log "====== TASK 6 COMPLETE: All visualizations saved ======"
log ""
log "=============================================="
log "PIPELINE COMPLETE"
log "Key outputs:"
log "  ${RESULTS_DIR}/stringtie_tpm_matrix.csv"
log "  ${RESULTS_DIR}/lncrna_differential_expression.csv"
log "  ${RESULTS_DIR}/significant_stress_responsive_lncrnas.csv"
log "  ${RESULTS_DIR}/volcano_plot.png"
log "  ${RESULTS_DIR}/heatmap_significant_lncrnas.png"
log "  ${RESULTS_DIR}/boxplots_top_lncrnas.png"
log "=============================================="
