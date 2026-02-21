import pandas as pd
import os

stringtie_dir = "/mnt/neil/resurrection_lncrna_pipeline/results/stringtie_output"
results_dir = "/mnt/neil/resurrection_lncrna_pipeline/results"

all_samples = [
    "cp_hyd_r1", "cp_hyd_r3",
    "cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3",
    "cp_wc2_r2", "cp_wc2_r3",
    "cp_wc60_r1", "cp_wc60_r2", "cp_wc60_r3"
]

dfs = []
missing = []

for sample in all_samples:
    abund_file = os.path.join(stringtie_dir, f"{sample}_abundance.txt")
    if not os.path.exists(abund_file) or os.path.getsize(abund_file) == 0:
        print(f"WARNING: Missing {sample}")
        missing.append(sample)
        continue
    df = pd.read_csv(abund_file, sep='\t')
    df = df[['Gene ID', 'Gene Name', 'TPM']].copy()
    df = df.rename(columns={'TPM': sample})
    dfs.append(df)
    print(f"Loaded: {sample} ({len(df)} genes)")

if missing:
    print(f"Missing: {missing}")

merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df.drop(columns=['Gene Name']), on='Gene ID', how='outer')
merged = merged.fillna(0)

out_path = os.path.join(results_dir, "stringtie_tpm_matrix.csv")
merged.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
print(f"Shape: {merged.shape}")
print(merged.head())
