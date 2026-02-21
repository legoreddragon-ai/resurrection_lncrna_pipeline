#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

PROJECT_DIR = Path.home() / "resurrection_lncrna_pipeline"

def load_expression_matrix(filename="stringtie_tpm_matrix.csv"):
    """Load TPM expression matrix"""
    filepath = PROJECT_DIR / "results" / filename
    if not filepath.exists():
        raise FileNotFoundError(f"Expression matrix not found: {filepath}")
    return pd.read_csv(filepath, index_col=0)

def load_coding_potential(filename="coding_potential.txt"):
    """Load coding potential classifications"""
    filepath = PROJECT_DIR / "results" / filename
    if not filepath.exists():
        raise FileNotFoundError(f"Coding potential file not found: {filepath}")
    return pd.read_csv(filepath, sep='\t')

def load_annotation(filename="GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt"):
    """Load annotation file"""
    filepath = PROJECT_DIR / "references" / filename
    if not filepath.exists():
        raise FileNotFoundError(f"Annotation file not found: {filepath}")
    return pd.read_csv(filepath, sep='\t')

def identify_lncrnas(tpm_matrix, coding_potential):
    """Filter lncRNAs from expression matrix"""
    lncrna_ids = set(coding_potential[
        coding_potential['Classification'] == 'lncRNA'
    ]['Transcript_ID'])
    
    # Filter matrix for lncRNAs
    lncrna_tpm = tpm_matrix[tpm_matrix.index.isin(lncrna_ids)]
    return lncrna_tpm

def differential_expression_analysis(lncrna_tpm):
    """Perform differential expression analysis"""
    
    print("Performing differential expression analysis...")
    
    # Define sample groups
    hydrated = [col for col in lncrna_tpm.columns 
                if 'hyd' in col and 'rehyd' not in col]
    rehydrated = [col for col in lncrna_tpm.columns if 'rehyd' in col]
    dehydrated = [col for col in lncrna_tpm.columns 
                  if 'wc2' in col or 'wc60' in col]
    
    print(f"  Hydrated samples: {hydrated}")
    print(f"  Rehydrated samples: {rehydrated}")
    print(f"  Dehydrated samples: {dehydrated}")
    
    # Log-transform
    log_tpm = np.log2(lncrna_tpm + 1)
    
    # Calculate means
    hydrated_mean = log_tpm[hydrated].mean(axis=1)
    rehydrated_mean = log_tpm[rehydrated].mean(axis=1)
    dehydrated_mean = log_tpm[dehydrated].mean(axis=1)
    
    # T-tests
    pvalues = []
    log2fcs = []
    
    for idx in log_tpm.index:
        dehyd_vals = log_tpm.loc[idx, dehydrated].values
        hyd_vals = log_tpm.loc[idx, hydrated].values
        
        mean_dehyd = dehyd_vals.mean() + 1
        mean_hyd = hyd_vals.mean() + 1
        
        t_stat, p_val = stats.ttest_ind(dehyd_vals, hyd_vals)
        log2fc = np.log2(mean_dehyd / mean_hyd)
        
        pvalues.append(p_val)
        log2fcs.append(log2fc)
    
    # Create results dataframe
    de_results = pd.DataFrame({
        'Gene_ID': log_tpm.index,
        'Hydrated_Mean_TPM': 2**hydrated_mean.values - 1,
        'Dehydrated_Mean_TPM': 2**dehydrated_mean.values - 1,
        'Rehydrated_Mean_TPM': 2**rehydrated_mean.values - 1,
        'Log2FC_Dehyd_vs_Hyd': log2fcs,
        'Pvalue': pvalues
    })
    
    # Multiple testing correction
    from scipy.stats import rankdata
    m = len(pvalues)
    ranks = rankdata(pvalues)
    de_results['AdjustedPvalue'] = de_results['Pvalue'] * m / ranks
    
    # Filter significant
    sig_lncrnas = de_results[
        (de_results['AdjustedPvalue'] < 0.05) & 
        (np.abs(de_results['Log2FC_Dehyd_vs_Hyd']) > 1)
    ]
    
    print(f"\nSignificant stress-responsive lncRNAs: {len(sig_lncrnas)}")
    print(f"  Upregulated in dehydration: {len(sig_lncrnas[sig_lncrnas['Log2FC_Dehyd_vs_Hyd'] > 1])}")
    print(f"  Downregulated in dehydration: {len(sig_lncrnas[sig_lncrnas['Log2FC_Dehyd_vs_Hyd'] < -1])}")
    
    return de_results, sig_lncrnas

def save_results(de_results, sig_lncrnas):
    """Save analysis results to CSV files"""
    
    results_dir = PROJECT_DIR / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    de_results.to_csv(results_dir / "lncrna_differential_expression.csv", index=False)
    sig_lncrnas.to_csv(results_dir / "significant_stress_responsive_lncrnas.csv", index=False)
    
    print(f"\nResults saved to:")
    print(f"  - {results_dir / 'lncrna_differential_expression.csv'}")
    print(f"  - {results_dir / 'significant_stress_responsive_lncrnas.csv'}")

def main():
    """Run complete analysis pipeline"""
    
    try:
        print("="*80)
        print("lncRNA DIFFERENTIAL EXPRESSION ANALYSIS")
        print("="*80)
        print()
        
        # Load data
        print("Loading data...")
        tpm_matrix = load_expression_matrix()
        coding_potential = load_coding_potential()
        print(f"  Expression matrix: {tpm_matrix.shape}")
        print(f"  Coding classifications: {len(coding_potential)}")
        print()
        
        # Identify lncRNAs
        print("Identifying lncRNAs...")
        lncrna_tpm = identify_lncrnas(tpm_matrix, coding_potential)
        print(f"  lncRNAs found: {len(lncrna_tpm)}")
        print()
        
        # Differential expression
        de_results, sig_lncrnas = differential_expression_analysis(lncrna_tpm)
        print()
        
        # Save results
        save_results(de_results, sig_lncrnas)
        print()
        
        print("="*80)
        print("ANALYSIS COMPLETE")
        print("="*80)
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
