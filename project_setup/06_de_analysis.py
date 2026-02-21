#!/usr/bin/env python3
# =============================================================================
# 06_de_analysis.py — Differential Expression Analysis
# =============================================================================
# Usage: python3 setup/06_de_analysis.py <BASE_DIR>
#
# Identifies stress-responsive lncRNAs by comparing hydrated vs dehydrated.
#
# Statistical approach:
#   - Welch's t-test (unequal variance, appropriate for unequal group sizes)
#   - Benjamini-Hochberg FDR correction
#   - NOTE: With only n=2 hydrated replicates, FDR yields 0 significant genes.
#     This is mathematically correct. Top 14 candidates are selected using
#     raw p-value + |log2FC| + expression magnitude as a combined score.
#     This limitation must be stated in any write-up.
#
# Outputs:
#   results/lncrna_differential_expression.csv  — all 17,574 expressed lncRNAs
#   results/significant_stress_responsive_lncrnas.csv — top 14 candidates
# =============================================================================

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from scipy.stats import rankdata

BASE_DIR    = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.parent.resolve()
RESULTS_DIR = BASE_DIR / "results"

TPM_FILE     = RESULTS_DIR / "stringtie_tpm_matrix.csv"
CODING_FILE  = RESULTS_DIR / "coding_potential.txt"
DE_FILE      = RESULTS_DIR / "lncrna_differential_expression.csv"
SIG_FILE     = RESULTS_DIR / "significant_stress_responsive_lncrnas.csv"

# Sample groups
HYDRATED_COLS   = ["cp_hyd_r1", "cp_hyd_r3"]
DEHYDRATED_COLS = ["cp_wc2_r2", "cp_wc2_r3", "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3"]
REHYDRATED_COLS = ["cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3"]
ALL_SAMPLE_COLS = HYDRATED_COLS + DEHYDRATED_COLS + REHYDRATED_COLS

print("=" * 60)
print("  Stage 6: Differential Expression Analysis")
print("=" * 60)

# ── Load TPM matrix ────────────────────────────────────────────────────────────
print(f"\nLoading TPM matrix: {TPM_FILE}")
tpm = pd.read_csv(TPM_FILE)
print(f"Shape: {tpm.shape}")

# Verify sample columns present
missing_cols = [c for c in ALL_SAMPLE_COLS if c not in tpm.columns]
if missing_cols:
    print(f"ERROR: Missing sample columns: {missing_cols}")
    sys.exit(1)

# ── Load lncRNA classification ─────────────────────────────────────────────────
print(f"\nLoading coding potential: {CODING_FILE}")
coding = pd.read_csv(CODING_FILE, sep="\t")
lncrna_ids = set(coding[coding["Classification"] == "lncRNA"]["Transcript_ID"])
print(f"lncRNA candidates: {len(lncrna_ids):,}")

tpm_lnc = tpm[tpm["Gene ID"].isin(lncrna_ids)].copy()
print(f"lncRNAs in TPM matrix: {len(tpm_lnc):,}")

# ── Expression filter ──────────────────────────────────────────────────────────
mask = (
    (tpm_lnc[HYDRATED_COLS].mean(axis=1)   >= 1) |
    (tpm_lnc[DEHYDRATED_COLS].mean(axis=1) >= 1) |
    (tpm_lnc[REHYDRATED_COLS].mean(axis=1) >= 1)
)
tpm_lnc = tpm_lnc[mask]
print(f"After expression filter (TPM>=1 in any condition): {len(tpm_lnc):,}")

# ── DE analysis ───────────────────────────────────────────────────────────────
print("\nRunning differential expression analysis...")
results = []

for _, row in tpm_lnc.iterrows():
    hyd_vals   = row[HYDRATED_COLS].values.astype(float)
    dehy_vals  = row[DEHYDRATED_COLS].values.astype(float)
    rehyd_vals = row[REHYDRATED_COLS].values.astype(float)

    mean_hyd   = np.mean(hyd_vals)
    mean_dehy  = np.mean(dehy_vals)
    mean_rehyd = np.mean(rehyd_vals)

    # Pseudocount of 0.1 to avoid log(0)
    log2fc = np.log2((mean_dehy + 0.1) / (mean_hyd + 0.1))

    # Welch's t-test (does not assume equal variance)
    _, pval = stats.ttest_ind(hyd_vals, dehy_vals, equal_var=False)

    results.append({
        "gene_id":         row["Gene ID"],
        "mean_hydrated":   round(mean_hyd,   3),
        "mean_dehydrated": round(mean_dehy,  3),
        "mean_rehydrated": round(mean_rehyd, 3),
        "log2FC":          round(log2fc,     4),
        "pvalue":          pval,
        **{s: round(row[s], 3) for s in ALL_SAMPLE_COLS},
    })

de_df = pd.DataFrame(results)

# ── BH FDR correction ──────────────────────────────────────────────────────────
pvals = de_df["pvalue"].fillna(1).values
n = len(pvals)
ranked = rankdata(pvals)
fdr = np.minimum(1, pvals * n / ranked)
idx = np.argsort(pvals)[::-1]
for i in range(len(idx) - 1):
    fdr[idx[i + 1]] = min(fdr[idx[i + 1]], fdr[idx[i]])
de_df["padj"] = fdr

# ── Save full DE results ───────────────────────────────────────────────────────
de_df.to_csv(DE_FILE, index=False)
print(f"Saved full DE results: {DE_FILE} ({len(de_df):,} lncRNAs)")

# ── Diagnostics ───────────────────────────────────────────────────────────────
print(f"\nP-value distribution:")
for thresh in [0.001, 0.01, 0.05, 0.1]:
    n_sig = (de_df["pvalue"] < thresh).sum()
    print(f"  pval < {thresh}: {n_sig:,}")

print(f"\nAfter FDR correction:")
for thresh in [0.05, 0.1, 0.2]:
    n_sig = (de_df["padj"] < thresh).sum()
    print(f"  padj < {thresh}: {n_sig:,}")

print(f"\n|log2FC| distribution:")
for thresh in [1, 2, 3]:
    n_sig = (de_df["log2FC"].abs() > thresh).sum()
    print(f"  |log2FC| > {thresh}: {n_sig:,}")

# ── Select top 14 candidates ───────────────────────────────────────────────────
print("\n--- Candidate Selection ---")
print("NOTE: FDR correction yields 0 significant genes with n=2 hydrated replicates.")
print("Selecting top 14 by: pvalue<0.05, |log2FC|>2, max_TPM>=5, ranked by score.")

candidates = de_df[
    (de_df["pvalue"] < 0.05) &
    (de_df["log2FC"].abs() > 2) &
    (de_df[["mean_hydrated", "mean_dehydrated", "mean_rehydrated"]].max(axis=1) >= 5)
].copy()
print(f"Candidates meeting filters: {len(candidates):,}")

if len(candidates) == 0:
    print("WARNING: No candidates found. Relaxing filters to pvalue<0.1, |log2FC|>1...")
    candidates = de_df[
        (de_df["pvalue"] < 0.1) &
        (de_df["log2FC"].abs() > 1) &
        (de_df[["mean_hydrated", "mean_dehydrated", "mean_rehydrated"]].max(axis=1) >= 1)
    ].copy()

# Score: fold change magnitude × log2(dehydrated expression)
candidates["score"] = (
    candidates["log2FC"].abs() *
    np.log2(candidates["mean_dehydrated"] + 1)
)
top14 = candidates.sort_values("score", ascending=False).head(14)
top14.to_csv(SIG_FILE, index=False)

print(f"\nTop 14 stress-responsive lncRNAs:")
print(top14[["gene_id", "mean_hydrated", "mean_dehydrated", "mean_rehydrated",
             "log2FC", "pvalue", "padj"]].to_string(index=False))

print(f"\nSaved: {SIG_FILE}")
print("\n[Stage 6 complete]")
