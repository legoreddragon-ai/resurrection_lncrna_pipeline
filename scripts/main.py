#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path.home() / "resurrection_lncrna_pipeline"

def main():
    print("="*80)
    print("RESURRECTION PLANT lncRNA DISCOVERY PIPELINE")
    print("="*80)
    print()
    
    print("Project directory:", PROJECT_DIR)
    print()
    
    # Check if project structure exists
    required_dirs = [
        PROJECT_DIR / "data" / "fastq_trimmed",
        PROJECT_DIR / "references",
        PROJECT_DIR / "results" / "stringtie_output",
        PROJECT_DIR / "scripts"
    ]
    
    print("Checking project structure...")
    for directory in required_dirs:
        if directory.exists():
            print(f"  [OK] {directory.relative_to(PROJECT_DIR)}")
        else:
            print(f"  [MISSING] {directory.relative_to(PROJECT_DIR)}")
    print()
    
    # Check for key files
    print("Checking key files...")
    key_files = [
        PROJECT_DIR / "README.md",
        PROJECT_DIR / "pipeline.sh",
        PROJECT_DIR / "LICENSE",
        PROJECT_DIR / ".gitignore"
    ]
    
    for file in key_files:
        if file.exists():
            size = file.stat().st_size
            print(f"  [OK] {file.name} ({size} bytes)")
        else:
            print(f"  [MISSING] {file.name}")
    print()
    
    # Check for results
    results_dir = PROJECT_DIR / "results"
    if results_dir.exists():
        csv_files = list(results_dir.glob("*.csv"))
        txt_files = list(results_dir.glob("*.txt"))
        
        if csv_files or txt_files:
            print("Found result files:")
            for file in csv_files + txt_files:
                size_mb = file.stat().st_size / (1024*1024)
                print(f"  - {file.name} ({size_mb:.2f} MB)")
        else:
            print("No result files found yet (analysis not run)")
    print()
    
    print("="*80)
    print("PIPELINE STATUS")
    print("="*80)
    print()
    print("Project initialized successfully!")
    print()
    print("Next steps:")
    print("1. Download FASTQ files from SRA (PRJNA660052)")
    print("2. Place trimmed FASTQ in data/fastq_trimmed/")
    print("3. Run: bash pipeline.sh")
    print()

if __name__ == "__main__":
    main()
