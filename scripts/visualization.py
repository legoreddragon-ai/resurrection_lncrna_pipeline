#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.resolve()

class Visualizer:
    """Create visualizations for lncRNA analysis"""
    
    def __init__(self, results_dir=None):
        if results_dir is None:
            results_dir = PROJECT_DIR / "results"
        self.results_dir = Path(results_dir)
        self.figures_dir = self.results_dir / "figures"
        self.figures_dir.mkdir(exist_ok=True)
        
        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.figsize'] = (12, 8)
    
    def volcano_plot(self, de_results, output_file="volcano_plot.png"):
        """Create volcano plot of differential expression"""
        
        print("Creating volcano plot...")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Determine colors
        colors = []
        for idx, row in de_results.iterrows():
            if row['AdjustedPvalue'] < 0.05 and abs(row['Log2FC_Dehyd_vs_Hyd']) > 1:
                colors.append('red')
            else:
                colors.append('gray')
        
        # Plot
        ax.scatter(de_results['Log2FC_Dehyd_vs_Hyd'], 
                  -np.log10(de_results['AdjustedPvalue']),
                  c=colors, alpha=0.6, s=50)
        
        # Add lines for thresholds
        ax.axhline(y=-np.log10(0.05), color='blue', linestyle='--', label='padj=0.05')
        ax.axvline(x=1, color='green', linestyle='--', label='log2FC=1')
        ax.axvline(x=-1, color='green', linestyle='--', label='log2FC=-1')
        
        ax.set_xlabel('Log2 Fold Change (Dehydrated vs Hydrated)', fontsize=12)
        ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
        ax.set_title('Volcano Plot: Stress-Responsive lncRNAs', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        output_path = self.figures_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved to {output_path}")
        plt.close()
    
    def expression_heatmap(self, tpm_matrix, sig_lncrnas, output_file="expression_heatmap.png"):
        """Create heatmap of expression levels"""
        
        print("Creating expression heatmap...")
        
        # Get significant lncRNA data
        sig_ids = sig_lncrnas['Gene_ID'].values if 'Gene_ID' in sig_lncrnas.columns else sig_lncrnas.index
        heatmap_data = tpm_matrix.loc[tpm_matrix.index.isin(sig_ids)]
        
        # Log-transform for better visualization
        heatmap_data_log = np.log2(heatmap_data + 1)
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(heatmap_data_log, cmap='RdYlBu_r', cbar_kws={'label': 'Log2(TPM+1)'}, ax=ax)
        
        ax.set_title('Expression Heatmap: Significant lncRNAs', fontsize=14)
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('lncRNAs', fontsize=12)
        
        output_path = self.figures_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved to {output_path}")
        plt.close()
    
    def expression_boxplot(self, tpm_matrix, sig_lncrnas, output_file="expression_boxplot.png"):
        """Create boxplot of expression by condition"""
        
        print("Creating expression boxplot...")
        
        # Prepare data
        conditions = []
        values = []
        
        for col in tpm_matrix.columns:
            if 'hyd' in col and 'rehyd' not in col:
                condition = 'Hydrated'
            elif 'rehyd' in col:
                condition = 'Rehydrated'
            elif 'wc' in col:
                condition = 'Dehydrated'
            else:
                continue
            
            sig_ids = sig_lncrnas['Gene_ID'].values if 'Gene_ID' in sig_lncrnas.columns else sig_lncrnas.index
            col_data = tpm_matrix.loc[tpm_matrix.index.isin(sig_ids), col]
            
            for val in col_data:
                conditions.append(condition)
                values.append(val)
        
        # Create dataframe
        plot_data = pd.DataFrame({
            'Condition': conditions,
            'TPM': values
        })
        
        # Create boxplot
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.boxplot(data=plot_data, x='Condition', y='TPM', ax=ax, palette='Set2')
        
        ax.set_ylabel('TPM', fontsize=12)
        ax.set_xlabel('Hydration State', fontsize=12)
        ax.set_title('Expression Levels Across Conditions', fontsize=14)
        ax.set_yscale('log')
        
        output_path = self.figures_dir / output_file
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved to {output_path}")
        plt.close()
    
    def summary_statistics(self, de_results, output_file="summary_statistics.txt"):
        """Generate summary statistics"""
        
        print("Generating summary statistics...")
        
        with open(self.figures_dir / output_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("DIFFERENTIAL EXPRESSION SUMMARY\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"Total lncRNAs analyzed: {len(de_results)}\n")
            
            sig = de_results[(de_results['AdjustedPvalue'] < 0.05) & 
                            (np.abs(de_results['Log2FC_Dehyd_vs_Hyd']) > 1)]
            f.write(f"Significant lncRNAs: {len(sig)}\n")
            
            upregulated = sig[sig['Log2FC_Dehyd_vs_Hyd'] > 1]
            downregulated = sig[sig['Log2FC_Dehyd_vs_Hyd'] < -1]
            
            f.write(f"  Upregulated in dehydration: {len(upregulated)}\n")
            f.write(f"  Downregulated in dehydration: {len(downregulated)}\n\n")
            
            f.write("Top 5 Downregulated lncRNAs:\n")
            top_down = de_results.nsmallest(5, 'Log2FC_Dehyd_vs_Hyd')
            for idx, row in top_down.iterrows():
                f.write(f"  {row['Gene_ID']}: log2FC={row['Log2FC_Dehyd_vs_Hyd']:.2f}, padj={row['AdjustedPvalue']:.2e}\n")
            
            f.write("\nTop 5 Upregulated lncRNAs:\n")
            top_up = de_results.nlargest(5, 'Log2FC_Dehyd_vs_Hyd')
            for idx, row in top_up.iterrows():
                f.write(f"  {row['Gene_ID']}: log2FC={row['Log2FC_Dehyd_vs_Hyd']:.2f}, padj={row['AdjustedPvalue']:.2e}\n")
        
        print(f"  Saved to {self.figures_dir / output_file}")

def main():
    """Test visualizer"""
    
    print("="*80)
    print("VISUALIZATION MODULE TEST")
    print("="*80)
    print()
    
    visualizer = Visualizer()
    print(f"Figures directory: {visualizer.figures_dir}")
    print()
    print("Visualization module ready!")

if __name__ == "__main__":
    main()
