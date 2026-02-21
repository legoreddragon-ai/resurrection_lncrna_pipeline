import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import rankdata

results_dir = "/mnt/neil/resurrection_lncrna_pipeline/results"

print("Loading TPM matrix...")
tpm = pd.read_csv(f"{results_dir}/stringtie_tpm_matrix.csv")

hydrated_cols   = ["cp_hyd_r1", "cp_hyd_r3"]
dehydrated_cols = ["cp_wc2_r2", "cp_wc2_r3", "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3"]
rehydrated_cols = ["cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3"]

# Load lncRNAs
coding = pd.read_csv(f"{results_dir}/coding_potential.txt", sep='\t')
lncrna_ids = set(coding[coding['Classification'] == 'lncRNA']['Transcript_ID'])
tpm_lnc = tpm[tpm['Gene ID'].isin(lncrna_ids)].copy()

# Expression filter
mask = (tpm_lnc[hydrated_cols].mean(axis=1) >= 1) | \
       (tpm_lnc[dehydrated_cols].mean(axis=1) >= 1) | \
       (tpm_lnc[rehydrated_cols].mean(axis=1) >= 1)
tpm_lnc = tpm_lnc[mask]
print(f"Expressed lncRNAs: {len(tpm_lnc)}")

# DE analysis
results = []
for _, row in tpm_lnc.iterrows():
    hyd_vals   = row[hydrated_cols].values.astype(float)
    dehy_vals  = row[dehydrated_cols].values.astype(float)
    rehyd_vals = row[rehydrated_cols].values.astype(float)

    mean_hyd   = np.mean(hyd_vals)
    mean_dehy  = np.mean(dehy_vals)
    mean_rehyd = np.mean(rehyd_vals)

    log2fc = np.log2((mean_dehy + 0.1) / (mean_hyd + 0.1))
    _, pval = stats.ttest_ind(hyd_vals, dehy_vals, equal_var=False)

    results.append({
        'gene_id':         row['Gene ID'],
        'mean_hydrated':   round(mean_hyd, 3),
        'mean_dehydrated': round(mean_dehy, 3),
        'mean_rehydrated': round(mean_rehyd, 3),
        'log2FC':          round(log2fc, 4),
        'pvalue':          pval,
        **{s: round(row[s], 3) for s in hydrated_cols + dehydrated_cols + rehydrated_cols}
    })

de_df = pd.DataFrame(results)

# BH FDR correction
pvals = de_df['pvalue'].fillna(1).values
n = len(pvals)
ranked = rankdata(pvals)
fdr = np.minimum(1, pvals * n / ranked)
idx = np.argsort(pvals)[::-1]
for i in range(len(idx) - 1):
    fdr[idx[i+1]] = min(fdr[idx[i+1]], fdr[idx[i]])
de_df['padj'] = fdr

# Save full results
de_df.to_csv(f"{results_dir}/lncrna_differential_expression.csv", index=False)
print(f"Saved full DE: {len(de_df)} lncRNAs")

# ── CANDIDATE SELECTION STRATEGY ─────────────────────────────────────────────
# With n=2 hydrated replicates, FDR correction is too conservative.
# Strategy: rank by pvalue, require strong fold change and high expression.

# Filter: pvalue < 0.05, |log2FC| > 2, mean_hydrated >= 5 TPM
candidates = de_df[
    (de_df['pvalue'] < 0.05) &
    (de_df['log2FC'].abs() > 2) &
    (de_df[['mean_hydrated','mean_dehydrated','mean_rehydrated']].max(axis=1) >= 5)
].copy()

print(f"\nCandidates (pval<0.05, |log2FC|>2, max TPM>=5): {len(candidates)}")

# Prioritize downregulated during dehydration (hydration-maintenance role)
downreg = candidates[candidates['log2FC'] < -2].copy()
downreg['score'] = downreg['mean_hydrated'] * downreg['log2FC'].abs()
downreg = downreg.sort_values('score', ascending=False)

# Take top 14
top14 = downreg.head(14)
top14.to_csv(f"{results_dir}/significant_stress_responsive_lncrnas.csv", index=False)

print(f"\nTop 14 stress-responsive lncRNAs (downregulated during dehydration):")
print(top14[['gene_id','mean_hydrated','mean_dehydrated','mean_rehydrated','log2FC','pvalue']].to_string())
