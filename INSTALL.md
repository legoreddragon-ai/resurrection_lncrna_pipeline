# Installation and Usage Guide

## Project Overview

Resurrection Plant lncRNA Discovery Pipeline - A complete computational pipeline for identifying and analyzing desiccation-induced long non-coding RNAs in resurrection plants.

## System Requirements

- Python 3.7+
- Ubuntu/Linux (or WSL on Windows)
- 4GB+ RAM
- 50GB+ disk space (for data)

## Installation

### Step 1: Clone or Download the Project
```bash
cd ~
# If cloning from GitHub:
git clone https://github.com/[YOUR_USERNAME]/resurrection_lncrna_pipeline.git

# Or if you've already created it:
cd resurrection_lncrna_pipeline
```

### Step 2: Install Python Dependencies
```bash
# Install required packages
pip install -r requirements.txt

# Optional: Install conda packages for bioinformatics tools
conda install -c bioconda trimmomatic hisat2 samtools stringtie
```

### Step 3: Verify Installation
```bash
# Run configuration test
python3 scripts/config.py

# Run data loader test
python3 scripts/data_loader.py

# Run main status check
python3 scripts/main.py
```

All tests should pass without errors.

## Project Structure
```
resurrection_lncrna_pipeline/
├── README.md                 # Project overview
├── INSTALL.md               # This file
├── LICENSE                  # MIT License
├── .gitignore               # Git ignore rules
├── requirements.txt         # Python dependencies
├── setup.py                 # Package setup
├── run.py                   # Main entry point
├── pipeline.sh              # Bash pipeline script
│
├── data/
│   ├── fastq_raw/          # Raw FASTQ files (from SRA)
│   └── fastq_trimmed/      # Quality-trimmed FASTQ files
│
├── references/
│   ├── annotation.gff3      # Gene annotation (GFF3)
│   ├── annotation.gtf       # Gene annotation (GTF)
│   ├── GSE157098_Cp_transcriptome_assembly_V2.fa
│   └── GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt
│
├── results/
│   ├── stringtie_output/    # StringTie assembly results
│   ├── figures/             # Generated visualizations
│   ├── stringtie_tpm_matrix.csv
│   ├── lncrna_differential_expression.csv
│   ├── significant_stress_responsive_lncrnas.csv
│   ├── coding_potential.txt
│   └── transcript_id_mapping.csv
│
├── scripts/
│   ├── __init__.py          # Python package init
│   ├── main.py              # Main status script
│   ├── config.py            # Configuration
│   ├── data_loader.py       # Data loading utilities
│   ├── analysis.py          # Analysis functions
│   └── visualization.py      # Visualization functions
│
└── docs/
    └── (Documentation files)
```

## Usage

### Option 1: Quick Status Check
```bash
python3 run.py
```

Shows project structure and available data.

### Option 2: Check Project Structure
```bash
python3 run.py --check-structure
```

Verifies all required directories exist.

### Option 3: Load and Display Data
```bash
python3 run.py --load-data
```

Loads available data files and displays statistics.

### Option 4: Run Analysis
```bash
python3 run.py --analyze
```

Performs differential expression analysis on available data.

### Option 5: Create Visualizations
```bash
python3 run.py --visualize
```

Creates plots from analysis results.

### Option 6: Run Complete Pipeline
```bash
python3 run.py --full
```

Runs all steps: check structure, load data, analyze, visualize.

### Option 7: Run Bash Pipeline
```bash
bash pipeline.sh
```

Runs the bash pipeline script.

## Downloading Data

### Download Annotation Files

The pipeline automatically downloads annotation files from GEO when needed.

### Download FASTQ Files from SRA

Download files from PRJNA660052:
```bash
# Using SRA Toolkit
prefetch PRJNA660052
fasterq-dump PRJNA660052

# Or download specific accessions
# Visit: https://www.ncbi.nlm.nih.gov/sra/PRJNA660052
```

Place downloaded FASTQ files in `data/fastq_raw/`

### Trim FASTQ Files
```bash
# Using Trimmomatic
trimmomatic PE -phred33 \
  data/fastq_raw/sample_1.fastq.gz data/fastq_raw/sample_2.fastq.gz \
  data/fastq_trimmed/sample_1.fastq.gz data/fastq_trimmed/sample_1_unpaired.fastq.gz \
  data/fastq_trimmed/sample_2.fastq.gz data/fastq_trimmed/sample_2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Place trimmed files in `data/fastq_trimmed/`

## Configuration

Edit `scripts/config.py` to customize:

- Minimum transcript length
- Coding potential threshold
- P-value threshold
- Log2 fold-change threshold
- Correlation threshold
- Analysis parameters

## Troubleshooting

### Import Errors

If you get "ModuleNotFoundError: No module named 'pandas'":
```bash
pip install pandas numpy scipy biopython matplotlib seaborn --break-system-packages
```

### Missing Data Files

Check that files are in correct locations:
```bash
ls -lh ~/resurrection_lncrna_pipeline/references/
ls -lh ~/resurrection_lncrna_pipeline/data/fastq_trimmed/
ls -lh ~/resurrection_lncrna_pipeline/results/
```

### Permission Denied

Make scripts executable:
```bash
chmod +x ~/resurrection_lncrna_pipeline/run.py
chmod +x ~/resurrection_lncrna_pipeline/pipeline.sh
chmod +x ~/resurrection_lncrna_pipeline/scripts/*.py
```

## Next Steps

1. Download FASTQ files from SRA (PRJNA660052)
2. Trim FASTQ files using Trimmomatic
3. Place trimmed files in `data/fastq_trimmed/`
4. Run `python3 run.py --full` to execute pipeline
5. Check `results/significant_stress_responsive_lncrnas.csv` for key findings

## Citation

If you use this pipeline, please cite:
```
Neil, Resurrection Plant lncRNA Discovery Pipeline, 2026
https://github.com/[YOUR_USERNAME]/resurrection_lncrna_pipeline
```

## License

This project is licensed under the MIT License - see LICENSE file for details.

## Support

For issues or questions, please create a GitHub issue or contact the author.
