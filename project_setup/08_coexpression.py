#!/usr/bin/env python3
"""
08_coexpression.py — Co-expression Network Analysis
Finds protein-coding genes that co-express with the 14 stress-responsive lncRNA candidates.

Usage:
    python3 project_setup/08_coexpression.py

Outputs:
    results/coexpression_partners.csv   — all gene pairs with r > threshold
    results/coexpression_summary.csv    — per-lncRNA partner count & top hits
    results/coexpression_heatmap.png    — correlation heatmap of 14 lncRNAs vs top partners
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import os
import sys

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MATRIX      = os.path.join(BASE_DIR, "results", "stringtie_tpm_matrix.csv")
CANDIDATES  = os.path.join(BASE_DIR, "results", "significant_stress_responsive_lncrnas.csv")
OUT_PAIRS   = os.path.join(BASE_DIR, "results", "coexpression_partners.csv")
OUT_SUMMARY = os.path.join(BASE_DIR, "results", "coexpression_summary.csv")
OUT_HEATMAP = os.path.join(BASE_DIR, "results", "coexpression_heatmap.png")

# Pearson r threshold — 0.9 is stringent; lower to 0.85 if too few partners
R_THRESHOLD = 0.90
MIN_EXPR    = 1.0   # at least 1 TPM in one sample to be tested

# Sample columns in the order they appear in the matrix
SAMPLE_COLS = [
    "cp_hyd_r1", "cp_hyd_r3",
    "cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3",
    "cp_wc2_r2", "cp_wc2_r3",
    "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3",
]

# Same samples in the order used by the DE results (for lncRNA vector extraction)
CAND_SAMPLE_COLS = [
    "cp_hyd_r1", "cp_hyd_r3",
    "cp_wc2_r2", "cp_wc2_r3", "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3",
    "cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3",
]

print("=" * 60)
print("Step 4A — Co-expression Network Analysis")
print("=" * 60)

# ── Load data ─────────────────────────────────────────────────────────────────
print("\n[1/5] Loading TPM matrix...")
matrix = pd.read_csv(MATRIX)
matrix = matrix.rename(columns={"Gene ID": "gene_id", "Gene Name": "gene_name"})
print(f"      {len(matrix):,} genes × {len(SAMPLE_COLS)} samples")

print("[2/5] Loading lncRNA candidates...")
cands = pd.read_csv(CANDIDATES)
print(f"      {len(cands)} candidates loaded")

# ── Build lncRNA expression vectors ──────────────────────────────────────────
# Map candidate sample cols → matrix sample cols (same names, different order)
lncrna_vectors = {}
for _, row in cands.iterrows():
    gid = row["gene_id"]
    # pull TPM in matrix column order
    vec = np.array([row[c] for c in SAMPLE_COLS if c in row.index], dtype=float)
    if len(vec) == len(SAMPLE_COLS):
        lncrna_vectors[gid] = vec
    else:
        # fallback: look up in matrix
        match = matrix[matrix["gene_id"] == gid]
        if not match.empty:
            lncrna_vectors[gid] = match[SAMPLE_COLS].values[0].astype(float)

print(f"      {len(lncrna_vectors)} lncRNA vectors ready")

# ── Filter matrix to expressed genes ─────────────────────────────────────────
print("[3/5] Filtering to expressed genes (TPM ≥ 1 in any sample)...")
expr_mask = (matrix[SAMPLE_COLS] >= MIN_EXPR).any(axis=1)
expressed = matrix[expr_mask].copy()
# Exclude the 14 lncRNA candidates themselves from the partner search
expressed = expressed[~expressed["gene_id"].isin(lncrna_vectors.keys())]
print(f"      {len(expressed):,} genes to test")

# ── Pearson correlation scan ──────────────────────────────────────────────────
print(f"[4/5] Computing Pearson correlations (r threshold = {R_THRESHOLD})...")
print(f"      Scanning {len(expressed):,} genes × {len(lncrna_vectors)} lncRNAs...")

records = []
gene_matrix = expressed[SAMPLE_COLS].values.astype(float)
gene_ids    = expressed["gene_id"].values
gene_names  = expressed["gene_name"].values if "gene_name" in expressed.columns else gene_ids

for lnc_id, lnc_vec in lncrna_vectors.items():
    # Skip lncRNAs with zero variance (all-zero vector)
    if lnc_vec.std() == 0:
        continue
    for i in range(len(gene_matrix)):
        g_vec = gene_matrix[i]
        if g_vec.std() == 0:
            continue
        r, pval = pearsonr(lnc_vec, g_vec)
        if r >= R_THRESHOLD:
            records.append({
                "lncrna_id":   lnc_id,
                "partner_id":  gene_ids[i],
                "partner_name": gene_names[i],
                "pearson_r":   round(r, 4),
                "pvalue":      round(pval, 6),
            })

pairs = pd.DataFrame(records)
print(f"      Found {len(pairs):,} co-expression pairs above threshold")

if len(pairs) == 0:
    print("\n  ⚠  No pairs found at r > 0.90. Retrying at r > 0.85...")
    R_THRESHOLD = 0.85
    records = []
    for lnc_id, lnc_vec in lncrna_vectors.items():
        if lnc_vec.std() == 0:
            continue
        for i in range(len(gene_matrix)):
            g_vec = gene_matrix[i]
            if g_vec.std() == 0:
                continue
            r, pval = pearsonr(lnc_vec, g_vec)
            if r >= R_THRESHOLD:
                records.append({
                    "lncrna_id":   lnc_id,
                    "partner_id":  gene_ids[i],
                    "partner_name": gene_names[i],
                    "pearson_r":   round(r, 4),
                    "pvalue":      round(pval, 6),
                })
    pairs = pd.DataFrame(records)
    print(f"      Found {len(pairs):,} co-expression pairs at r > 0.85")

pairs = pairs.sort_values(["lncrna_id", "pearson_r"], ascending=[True, False])
pairs.to_csv(OUT_PAIRS, index=False)
print(f"      Saved → {OUT_PAIRS}")

# ── Summary per lncRNA ────────────────────────────────────────────────────────
summary_rows = []
for lnc_id in lncrna_vectors.keys():
    sub = pairs[pairs["lncrna_id"] == lnc_id]
    top = sub.head(3)["partner_name"].tolist() if not sub.empty else []
    summary_rows.append({
        "lncrna_id":     lnc_id,
        "n_partners":    len(sub),
        "mean_r":        round(sub["pearson_r"].mean(), 4) if not sub.empty else 0,
        "max_r":         round(sub["pearson_r"].max(), 4)  if not sub.empty else 0,
        "top_partners":  "; ".join(top),
    })

summary = pd.DataFrame(summary_rows).sort_values("n_partners", ascending=False)
summary.to_csv(OUT_SUMMARY, index=False)
print(f"      Saved → {OUT_SUMMARY}")

# ── Print summary table ───────────────────────────────────────────────────────
print("\n── Co-expression Summary ──────────────────────────────────")
print(f"{'lncRNA':<28} {'Partners':>8} {'Mean r':>8} {'Top partner'}")
print("-" * 70)
for _, row in summary.iterrows():
    top = row["top_partners"].split(";")[0].strip() if row["top_partners"] else "—"
    print(f"{row['lncrna_id']:<28} {int(row['n_partners']):>8} {row['mean_r']:>8.4f}  {top}")

# ── Heatmap of top partners ───────────────────────────────────────────────────
print(f"\n[5/5] Generating co-expression heatmap...")

if len(pairs) > 0:
    # Take top 30 partner genes by max r across all lncRNAs
    top_partners = (pairs.groupby("partner_id")["pearson_r"]
                    .max()
                    .sort_values(ascending=False)
                    .head(30)
                    .index.tolist())

    # Build r matrix: lncRNAs × top partners
    lnc_ids = list(lncrna_vectors.keys())
    r_matrix = pd.DataFrame(index=lnc_ids, columns=top_partners, dtype=float)

    for lnc_id in lnc_ids:
        for pid in top_partners:
            row = pairs[(pairs["lncrna_id"] == lnc_id) & (pairs["partner_id"] == pid)]
            r_matrix.loc[lnc_id, pid] = row["pearson_r"].values[0] if not row.empty else np.nan

    # Shorten gene IDs for display
    r_matrix.index   = [i.replace("Cp_V2_contig_", "contig_") for i in r_matrix.index]
    r_matrix.columns = [c.replace("Cp_V2_contig_", "contig_") for c in r_matrix.columns]

    fig, ax = plt.subplots(figsize=(max(14, len(top_partners) * 0.45), 7))
    sns.heatmap(
        r_matrix.astype(float),
        ax=ax,
        cmap="RdYlGn",
        vmin=0.8, vmax=1.0,
        linewidths=0.3,
        linecolor="#cccccc",
        annot=False,
        cbar_kws={"label": "Pearson r", "shrink": 0.7},
        mask=r_matrix.isnull(),
    )
    ax.set_title(
        f"Co-expression Network: 14 lncRNAs × Top {len(top_partners)} Partners\n"
        f"Craterostigma plantagineum  (Pearson r ≥ {R_THRESHOLD})",
        fontsize=13, fontweight="bold", pad=15,
    )
    ax.set_xlabel("Co-expressed Partner Genes", fontsize=10)
    ax.set_ylabel("lncRNA Candidates", fontsize=10)
    ax.tick_params(axis="x", labelsize=7, rotation=90)
    ax.tick_params(axis="y", labelsize=8, rotation=0)
    plt.tight_layout()
    plt.savefig(OUT_HEATMAP, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"      Saved → {OUT_HEATMAP}")
else:
    print("      Skipping heatmap — no partners found")

print("\n" + "=" * 60)
print("Step 4A Complete")
print(f"  Co-expression pairs : {len(pairs):,}")
print(f"  lncRNAs with partners: {(summary['n_partners'] > 0).sum()} / {len(summary)}")
print(f"  Next step: python3 project_setup/09_go_enrichment.py")
print("=" * 60)
