#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.resolve()

class DataLoader:
    """Load and manage project data"""
    
    def __init__(self):
        self.project_dir = PROJECT_DIR
        self.data_dir = self.project_dir / "data"
        self.results_dir = self.project_dir / "results"
        self.references_dir = self.project_dir / "references"
    
    def load_tpm_matrix(self, filename="stringtie_tpm_matrix.csv"):
        """Load TPM expression matrix"""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"Warning: TPM matrix not found at {filepath}")
            return None
        
        print(f"Loading TPM matrix from {filename}...")
        df = pd.read_csv(filepath, index_col=0)
        print(f"  Genes: {df.shape[0]}, Samples: {df.shape[1]}")
        return df
    
    def load_coding_potential(self, filename="coding_potential.txt"):
        """Load coding potential classifications"""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"Warning: Coding potential file not found at {filepath}")
            return None
        
        print(f"Loading coding potential from {filename}...")
        df = pd.read_csv(filepath, sep='\t')
        print(f"  Total transcripts: {len(df)}")
        lncrna_count = len(df[df['Classification'] == 'lncRNA'])
        coding_count = len(df[df['Classification'] == 'coding'])
        print(f"  lncRNAs: {lncrna_count}, Coding: {coding_count}")
        return df
    
    def load_annotation(self, filename="GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt"):
        """Load annotation file"""
        filepath = self.references_dir / filename
        if not filepath.exists():
            print(f"Warning: Annotation file not found at {filepath}")
            return None
        
        print(f"Loading annotation from {filename}...")
        df = pd.read_csv(filepath, sep='\t')
        print(f"  Annotated transcripts: {len(df)}")
        return df
    
    def load_differential_expression(self, filename="lncrna_differential_expression.csv"):
        """Load differential expression results"""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"Warning: DE results not found at {filepath}")
            return None
        
        print(f"Loading DE results from {filename}...")
        df = pd.read_csv(filepath)
        print(f"  Total lncRNAs analyzed: {len(df)}")
        sig_count = len(df[(df['AdjustedPvalue'] < 0.05) & (np.abs(df['Log2FC_Dehyd_vs_Hyd']) > 1)])
        print(f"  Significant lncRNAs: {sig_count}")
        return df
    
    def load_significant_lncrnas(self, filename="significant_stress_responsive_lncrnas.csv"):
        """Load significant stress-responsive lncRNAs"""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"Warning: Significant lncRNAs not found at {filepath}")
            return None
        
        print(f"Loading significant lncRNAs from {filename}...")
        df = pd.read_csv(filepath)
        print(f"  Significant lncRNAs found: {len(df)}")
        return df
    
    def load_id_mapping(self, filename="transcript_id_mapping.csv"):
        """Load transcript ID mapping"""
        filepath = self.results_dir / filename
        if not filepath.exists():
            print(f"Warning: ID mapping not found at {filepath}")
            return None
        
        print(f"Loading ID mapping from {filename}...")
        df = pd.read_csv(filepath)
        print(f"  Transcript mappings: {len(df)}")
        return df
    
    def check_project_structure(self):
        """Check if all required directories exist"""
        print("\nChecking project structure...")
        
        required_dirs = [
            self.data_dir / "fastq_trimmed",
            self.results_dir / "stringtie_output",
            self.references_dir
        ]
        
        all_exist = True
        for directory in required_dirs:
            if directory.exists():
                print(f"  [OK] {directory.relative_to(self.project_dir)}")
            else:
                print(f"  [MISSING] {directory.relative_to(self.project_dir)}")
                all_exist = False
        
        return all_exist
    
    def list_results(self):
        """List all available result files"""
        print("\nAvailable result files:")
        
        if not self.results_dir.exists():
            print("  Results directory not found")
            return
        
        csv_files = list(self.results_dir.glob("*.csv"))
        txt_files = list(self.results_dir.glob("*.txt"))
        fa_files = list(self.results_dir.glob("*.fa"))
        
        for file in sorted(csv_files + txt_files + fa_files):
            size_mb = file.stat().st_size / (1024*1024)
            print(f"  - {file.name} ({size_mb:.2f} MB)")

def main():
    """Test data loader"""
    
    print("="*80)
    print("DATA LOADER TEST")
    print("="*80)
    print()
    
    loader = DataLoader()
    
    loader.check_project_structure()
    print()
    loader.list_results()
    print()
    
    print("="*80)
    print("DATA LOADER READY")
    print("="*80)

if __name__ == "__main__":
    main()
