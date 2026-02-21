#!/usr/bin/env python3
# =============================================================================
# 05_merge_matrix.py — Merge StringTie Abundance Files into TPM Matrix
# =============================================================================
# Usage: python3 setup/05_merge_matrix.py <BASE_DIR>
# Combines all 10 abundance.txt files into a single expression matrix.
# Output: results/stringtie_tpm_matrix.csv  (48045 genes × 12 columns)
# =============================================================================

import sys
import os
import pandas as pd
from pathlib import Path

BASE_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.parent.resolve()
ST_DIR   = BASE_DIR / "results" / "stringtie_output"
OUT_FILE = BASE_DIR / "results" / "stringtie_tpm_matrix.csv"

ALL_SAMPLES = [
    "cp_hyd_r1",   "cp_hyd_r3",
    "cp_rehyd_r1", "cp_rehyd_r2", "cp_rehyd_r3",
    "cp_wc2_r2",   "cp_wc2_r3",
    "cp_wc60_r1",  "cp_wc60_r2",  "cp_wc60_r3",
]

print("=" * 60)
print("  Stage 5: Merging TPM Expression Matrix")
print("=" * 60)

dfs = []
missing = []

for sample in ALL_SAMPLES:
    abund_file = ST_DIR / f"{sample}_abundance.txt"

    if not abund_file.exists() or abund_file.stat().st_size == 0:
        print(f"  MISSING: {sample}")
        missing.append(sample)
        continue

    df = pd.read_csv(abund_file, sep="\t")

    # Sanity check: must have standard StringTie columns
    required = {"Gene ID", "Gene Name", "TPM"}
    if not required.issubset(df.columns):
        print(f"  ERROR: {sample} has unexpected columns: {list(df.columns)}")
        sys.exit(1)

    # Sanity check: must have full gene count
    if len(df) < 1000:
        print(f"  ERROR: {sample} has only {len(df)} rows — "
              f"StringTie likely ran without -G flag. Rerun 04_stringtie.sh.")
        sys.exit(1)

    df = df[["Gene ID", "Gene Name", "TPM"]].copy()
    df = df.rename(columns={"TPM": sample})
    dfs.append(df)
    print(f"  Loaded: {sample:20s} ({len(df):,} genes)")

if missing:
    print(f"\nMissing samples: {missing}")
    print("Proceeding with available samples...")

if not dfs:
    print("ERROR: No abundance files found. Cannot build matrix.")
    sys.exit(1)

# Merge all samples on Gene ID
print(f"\nMerging {len(dfs)} samples...")
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df.drop(columns=["Gene Name"]), on="Gene ID", how="outer")

merged = merged.fillna(0)
merged.to_csv(OUT_FILE, index=False)

print(f"\nSaved: {OUT_FILE}")
print(f"Shape: {merged.shape[0]:,} genes × {merged.shape[1]} columns")
print(f"Columns: {list(merged.columns)}")
print(f"\nSample statistics:")
sample_cols = [c for c in merged.columns if c.startswith("cp_")]
print(merged[sample_cols].describe().round(2).to_string())
print("\n[Stage 5 complete]")
