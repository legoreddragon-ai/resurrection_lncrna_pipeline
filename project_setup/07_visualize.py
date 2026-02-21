#!/usr/bin/env python3
# =============================================================================
# 07_visualize.py — Generate All Result Visualizations
# =============================================================================
# Usage: python3 setup/07_visualize.py <BASE_DIR>
# Generates:
#   results/volcano_plot.png
#   results/heatmap_significant_lncrnas.png
#   results/boxplots_top_lncrnas.png
# =============================================================================

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE_DIR    = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.parent.resolve()
RESULTS_DIR = BASE_DIR / "results"
FIG_DIR     = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

DE_FILE  = RESULTS_DIR / "lncrna_differential_expression.csv"
SIG_FILE = RESULTS_DIR / "significant_stress_responsive_lncrnas.csv"
TPM_FILE = RESULTS_DIR / "stringtie_tpm_matrix.csv"

HYDRATED_COLS   = ["cp_hyd_r1", "cp_hyd_r3"]
DEHYDRATED_COLS = ["cp_wc2_r2", "cp_wc2_r3", "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3"]
REHYDRATED_COLS = ["cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3"]

print("=" * 60)
print("  Stage 7: Generating Visualizations")
print("=" * 60)

# ── Load data ─────────────────────────────────────────────────────────────────
print("\nLoading results...")
de_df  = pd.read_csv(DE_FILE)
sig_df = pd.read_csv(SIG_FILE)
tpm_df = pd.read_csv(TPM_FILE)
print(f"DE results: {len(de_df):,} | Significant: {len(sig_df)} | TPM matrix: {tpm_df.shape}")

# ── 1. VOLCANO PLOT ───────────────────────────────────────────────────────────
print("\nGenerating volcano plot...")
fig, ax = plt.subplots(figsize=(11, 8))

de_df["neg_log10_padj"] = -np.log10(de_df["pvalue"].clip(lower=1e-300))

nonsig  = de_df[(de_df["pvalue"] >= 0.05) | (de_df["log2FC"].abs() <= 1)]
sig_up  = de_df[(de_df["pvalue"] < 0.05)  & (de_df["log2FC"] > 1)]
sig_dn  = de_df[(de_df["pvalue"] < 0.05)  & (de_df["log2FC"] < -1)]

ax.scatter(nonsig["log2FC"],  nonsig["neg_log10_padj"],  s=5,  color="#CCCCCC", alpha=0.4, label="Not significant")
ax.scatter(sig_up["log2FC"],  sig_up["neg_log10_padj"],  s=18, color="#E74C3C", alpha=0.8, label=f"Up in dehydrated (n={len(sig_up)})")
ax.scatter(sig_dn["log2FC"],  sig_dn["neg_log10_padj"],  s=18, color="#2E86C1", alpha=0.8, label=f"Down in dehydrated (n={len(sig_dn)})")

# Label top 14 candidates
if len(sig_df) > 0:
    for _, row in sig_df.iterrows():
        gene_row = de_df[de_df["gene_id"] == row["gene_id"]]
        if gene_row.empty:
            continue
        x = gene_row["log2FC"].values[0]
        y = -np.log10(gene_row["pvalue"].values[0] + 1e-300)
        label = str(row["gene_id"]).replace("Cp_V2_contig_", "contig_")
        ax.annotate(label, (x, y), fontsize=6, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points", color="#1A1A1A")

ax.axvline(x=1,  color="#888888", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(x=-1, color="#888888", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(y=-np.log10(0.05), color="#888888", linestyle="--", linewidth=0.8, alpha=0.7)

ax.set_xlabel("log₂ Fold Change (Dehydrated / Hydrated)", fontsize=12)
ax.set_ylabel("-log₁₀ p-value", fontsize=12)
ax.set_title("Volcano Plot: lncRNA Differential Expression\nHydrated vs. Dehydrated — Craterostigma plantagineum", fontsize=13)
ax.legend(frameon=True, fontsize=10)
ax.grid(True, alpha=0.2)
plt.tight_layout()

out = RESULTS_DIR / "volcano_plot.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {out}")

# ── 2. HEATMAP ────────────────────────────────────────────────────────────────
print("Generating heatmap...")
sample_cols = [c for c in tpm_df.columns if c.startswith("cp_")]

if len(sig_df) > 0 and sample_cols:
    heat_data = tpm_df[tpm_df["Gene ID"].isin(sig_df["gene_id"])].copy()
    heat_data = heat_data.set_index("Gene ID")[sample_cols]

    # Log2 transform then z-score per gene
    heat_log = np.log2(heat_data + 1)
    heat_z = heat_log.subtract(heat_log.mean(axis=1), axis=0)
    std = heat_log.std(axis=1).replace(0, 1)
    heat_z = heat_z.divide(std, axis=0)

    # Order columns: hydrated → dehydrated → rehydrated
    col_order = sorted(sample_cols, key=lambda x: (
        0 if ("hyd" in x and "rehyd" not in x and "wc" not in x) else
        1 if "wc" in x else 2
    ))
    heat_z = heat_z[col_order]

    # Shorten row labels
    heat_z.index = heat_z.index.str.replace("Cp_V2_contig_", "contig_")

    fig, ax = plt.subplots(figsize=(13, max(6, len(heat_z) * 0.55 + 2)))
    im = ax.imshow(heat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)

    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_order, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(heat_z)))
    ax.set_yticklabels(heat_z.index, fontsize=8)

    # Column group labels
    n_hyd  = sum(1 for c in col_order if "hyd" in c and "rehyd" not in c and "wc" not in c)
    n_dehy = sum(1 for c in col_order if "wc" in c)
    n_rehyd= sum(1 for c in col_order if "rehyd" in c)
    for start, count, label, color in [
        (0,               n_hyd,  "Hydrated",    "#2ECC71"),
        (n_hyd,           n_dehy, "Dehydrated",  "#E74C3C"),
        (n_hyd + n_dehy,  n_rehyd,"Rehydrated",  "#3498DB"),
    ]:
        ax.axvline(x=start - 0.5, color=color, linewidth=2, alpha=0.6)
        ax.text(start + count / 2 - 0.5, -1.2, label,
                ha="center", va="top", fontsize=9, color=color, fontweight="bold",
                transform=ax.get_xaxis_transform())

    cb = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("Z-score (log₂ TPM)", fontsize=10)
    ax.set_title(f"Expression Heatmap: {len(sig_df)} Stress-Responsive lncRNAs\nCraterostigma plantagineum", fontsize=12)
    plt.tight_layout()

    out = RESULTS_DIR / "heatmap_significant_lncrnas.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

# ── 3. BOXPLOTS ───────────────────────────────────────────────────────────────
print("Generating boxplots...")
if len(sig_df) > 0:
    top14 = sig_df.head(14)
    n = len(top14)
    ncols = 3
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(14, nrows * 3.8))
    axes = axes.flatten()

    for i, (_, row) in enumerate(top14.iterrows()):
        gene_row = tpm_df[tpm_df["Gene ID"] == row["gene_id"]]
        if gene_row.empty:
            continue
        gene_row = gene_row.iloc[0]

        hyd_vals   = gene_row[HYDRATED_COLS].values.astype(float)
        dehy_vals  = gene_row[DEHYDRATED_COLS].values.astype(float)
        rehyd_vals = gene_row[REHYDRATED_COLS].values.astype(float)

        ax = axes[i]
        bp = ax.boxplot(
            [hyd_vals, dehy_vals, rehyd_vals],
            labels=["Hydrated", "Dehydrated", "Rehydrated"],
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            flierprops=dict(marker="o", markersize=4),
        )
        for patch, color in zip(bp["boxes"], ["#2ECC71", "#E74C3C", "#3498DB"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)

        label = str(row["gene_id"]).replace("Cp_V2_contig_", "contig_")
        ax.set_title(f"{label}\nlog₂FC={row['log2FC']:.2f}  p={row['pvalue']:.3f}",
                     fontsize=8)
        ax.set_ylabel("TPM", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3, axis="y")

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("Expression by Condition: Top 14 Stress-Responsive lncRNAs\nCraterostigma plantagineum",
                 fontsize=12, y=1.01)
    plt.tight_layout()

    out = RESULTS_DIR / "boxplots_top_lncrnas.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out}")

print("\n[Stage 7 complete]")
print("\nAll outputs:")
for f in ["volcano_plot.png", "heatmap_significant_lncrnas.png", "boxplots_top_lncrnas.png"]:
    p = RESULTS_DIR / f
    if p.exists():
        print(f"  {p}  ({p.stat().st_size // 1024} KB)")
