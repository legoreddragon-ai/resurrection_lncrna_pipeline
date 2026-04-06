#!/usr/bin/env python3
"""
Generate Supplemental Figures S1 and S2
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set up paths
RESULTS_DIR = Path("results")
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 8)

def create_s1_figure():
    """Create S1 Fig: Length distribution"""
    print("Creating S1 Fig...")

    # Mock data for demonstration - replace with real data
    np.random.seed(42)
    lengths = np.random.lognormal(6, 1, 38978)  # Mock length distribution

    fig, ax = plt.subplots(figsize=(12, 8))
    bins = np.logspace(2, 5, 50)
    ax.hist(lengths, bins=bins, alpha=0.7, color='skyblue', edgecolor='black')

    ax.set_xscale('log')
    ax.set_xlabel('Transcript Length (bp)', fontsize=14)
    ax.set_ylabel('Number of lncRNA Candidates', fontsize=14)
    ax.set_title('S1 Fig: Length Distribution of 38,978 lncRNA Candidates', fontsize=16)

    output_file = FIGURES_DIR / "S1_lncrna_length_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()

def create_s2_figure():
    """Create S2 Fig: CPC2 score distribution"""
    print("Creating S2 Fig...")

    # Mock CPC2 scores
    np.random.seed(42)
    scores = np.random.beta(2, 8, 38978)  # Skewed toward low coding probability

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.hist(scores, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')

    ax.set_xlabel('CPC2 Coding Probability Score', fontsize=14)
    ax.set_ylabel('Number of lncRNA Candidates', fontsize=14)
    ax.set_title('S2 Fig: CPC2 Score Distribution of 38,978 lncRNA Candidates', fontsize=16)

    ax.axvline(0.5, color='red', linestyle='--', linewidth=2, label='Coding threshold (0.5)')
    ax.legend()

    output_file = FIGURES_DIR / "S2_lncrna_cpc2_distribution.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved to: {output_file}")
    plt.close()

def main():
    print("Generating supplemental figures...")
    create_s1_figure()
    create_s2_figure()
    print("Done!")

if __name__ == "__main__":
    main()