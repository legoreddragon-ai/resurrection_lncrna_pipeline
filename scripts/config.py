#!/usr/bin/env python3

from pathlib import Path

# Project directories
PROJECT_DIR = Path(__file__).parent.parent.resolve()
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
REFERENCES_DIR = PROJECT_DIR / "references"
SCRIPTS_DIR = PROJECT_DIR / "scripts"
DOCS_DIR = PROJECT_DIR / "docs"

# Data directories
FASTQ_RAW_DIR = DATA_DIR / "fastq_raw"
FASTQ_TRIMMED_DIR = DATA_DIR / "fastq_trimmed"

# Results directories
STRINGTIE_OUTPUT_DIR = RESULTS_DIR / "stringtie_output"
HISAT2_ALIGNED_DIR = RESULTS_DIR / "hisat2_aligned"
FIGURES_DIR = RESULTS_DIR / "figures"

# Input files
ANNOTATION_FILE = REFERENCES_DIR / "GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt"
TRANSCRIPTOME_FASTA = REFERENCES_DIR / "GSE157098_Cp_transcriptome_assembly_V2.fa"
ANNOTATION_GFF3 = REFERENCES_DIR / "annotation.gff3"
ANNOTATION_GTF = REFERENCES_DIR / "annotation.gtf"

# Result files
TPM_MATRIX_FILE = RESULTS_DIR / "stringtie_tpm_matrix.csv"
CODING_POTENTIAL_FILE = RESULTS_DIR / "coding_potential.txt"
ID_MAPPING_FILE = RESULTS_DIR / "transcript_id_mapping.csv"
DE_RESULTS_FILE = RESULTS_DIR / "lncrna_differential_expression.csv"
SIGNIFICANT_LNCRNAS_FILE = RESULTS_DIR / "significant_stress_responsive_lncrnas.csv"
COEXPRESSION_NETWORK_FILE = RESULTS_DIR / "lncrna_coexpression_network.csv"

# Analysis parameters
TRANSCRIPT_MIN_LENGTH = 200  # Minimum transcript length in nucleotides
CODING_SCORE_THRESHOLD = 0.5  # Threshold for coding vs lncRNA classification
CORRELATION_THRESHOLD = 0.7  # Correlation threshold for co-expression
P_VALUE_THRESHOLD = 0.05  # P-value threshold for significance
LOG2FC_THRESHOLD = 1  # Log2 fold-change threshold
FDR_METHOD = "benjamini_hochberg"  # Multiple testing correction method

# Sample grouping
SAMPLE_GROUPS = {
    'hydrated': ['hyd'],
    'rehydrated': ['rehyd'],
    'dehydrated': ['wc2', 'wc60']
}

# Database URLs
GEO_ACCESSION = "GSE157098"
SRA_PROJECT = "PRJNA660052"
GEO_FTP = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/"

# Visualization parameters
FIGURE_DPI = 300
FIGURE_STYLE = "whitegrid"
FIGURE_SIZE_SMALL = (10, 8)
FIGURE_SIZE_LARGE = (14, 10)

# Logging
LOG_LEVEL = "INFO"
LOG_FILE = RESULTS_DIR / "pipeline.log"

def create_directories():
    """Create all required directories if they don't exist"""
    for directory in [DATA_DIR, RESULTS_DIR, REFERENCES_DIR, SCRIPTS_DIR, 
                     DOCS_DIR, FASTQ_RAW_DIR, FASTQ_TRIMMED_DIR,
                     STRINGTIE_OUTPUT_DIR, HISAT2_ALIGNED_DIR, FIGURES_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

def print_config():
    """Print current configuration"""
    print("="*80)
    print("PROJECT CONFIGURATION")
    print("="*80)
    print()
    
    print("Directories:")
    print(f"  Project: {PROJECT_DIR}")
    print(f"  Data: {DATA_DIR}")
    print(f"  Results: {RESULTS_DIR}")
    print(f"  References: {REFERENCES_DIR}")
    print()
    
    print("Analysis Parameters:")
    print(f"  Min transcript length: {TRANSCRIPT_MIN_LENGTH} bp")
    print(f"  Coding threshold: {CODING_SCORE_THRESHOLD}")
    print(f"  P-value threshold: {P_VALUE_THRESHOLD}")
    print(f"  Log2FC threshold: {LOG2FC_THRESHOLD}")
    print(f"  Correlation threshold: {CORRELATION_THRESHOLD}")
    print()
    
    print("Data Sources:")
    print(f"  GEO Accession: {GEO_ACCESSION}")
    print(f"  SRA Project: {SRA_PROJECT}")
    print()

def main():
    """Test configuration"""
    print_config()
    create_directories()
    print("Directories created/verified.")

if __name__ == "__main__":
    main()
