# lncRNA Discovery Pipeline
Computational pipeline for identifying desiccation-associated long non-coding RNAs (lncRNAs) in *Craterostigma plantagineum* and comparing response patterns with *Selaginella lepidophylla* using public RNA-seq data.

This repository includes the full analysis workflow from raw data download through alignment, transcript quantification, differential expression analysis, co-expression analysis, WGCNA, BLAST-based annotation, GO enrichment, and downstream validation steps.

## Associated manuscript

**Discovery of Desiccation-Induced Long Non-Coding RNAs in the Resurrection Plant *Craterostigma plantagineum*: A Transcriptome-Wide Computational Analysis with Cross-Species Validation in *Selaginella lepidophylla***  
**Author:** Neil Sumanth

## Datasets

- *Craterostigma plantagineum*: GEO **GSE157098** / SRA **PRJNA660052**
- *Selaginella lepidophylla*: SRA **PRJNA420971**

## Overview

The pipeline was developed to identify lncRNAs associated with desiccation response in *C. plantagineum* and to test whether related functional patterns also appear in *S. lepidophylla*.

Main analysis components:

- RNA-seq download and preprocessing
- read trimming
- HISAT2 alignment
- StringTie transcript quantification
- TPM matrix generation
- differential expression analysis
- visualization
- Pearson co-expression analysis
- WGCNA module analysis
- BLAST annotation of co-expressed partners
- GO enrichment
- Rfam contamination screening
- CPC2 coding potential assessment
- miRNA target / sponge prediction
- cross-species comparison

## Key results

### *Craterostigma plantagineum*

- Transcripts quantified: **48,045**
- lncRNA candidates identified: **38,978** (**81.1%**)
- lncRNAs tested for differential expression: **17,574**
- Core DESeq2-supported candidates: **14**
- DESeq2 adjusted p-value range for core candidates: **4.5 × 10^-7 to 1.7 × 10^-6**
- Exploratory upregulated candidates: **661**
- Exploratory downregulated candidates: **0**
- Co-expression pairs: **1,062** at **r ≥ 0.90**
- Annotated co-expressed partner genes: **356**
- Rfam contamination hits: **0**
- Predicted miRNA interactions: **106**
- Arabidopsis BLAST hits: **0**

### *Selaginella lepidophylla*

- Upregulated candidates: **303**
- Downregulated candidates: **98**
- Chloroplast/photosynthesis enrichment: **1.99×**, **FDR = 3 × 10^-22**

### Shared biological pattern

Both species showed enrichment of chloroplast/photosynthesis-related functions among stress-associated transcripts:

- *C. plantagineum*: **3.73×**, **FDR = 6 × 10^-14**
- *S. lepidophylla*: **1.99×**, **FDR = 3 × 10^-22**

### Candidate highlight

One of the strongest candidates in the *C. plantagineum* analysis was **contig_2815**, which emerged as:

- a WGCNA hub in the **MEred** module
- a predicted **miR164a** ceRNA target candidate
- associated with **ERD7-like** annotation
- co-expressed with stress-related partners including **PP2C-60**, **WRKY3**, **MYB73**, and **ABI5**

## Repository structure

```text
project_setup/
├── 00_preflight.sh
├── 01_download.sh
├── 02_trim.sh
├── 03_align.sh
├── 04_stringtie.sh
├── 05_merge_matrix.py
├── 06_de_analysis.py
├── 07_visualize.py
├── 08_coexpression.py
├── 09_blast_annotate.py
├── 10_go_enrichment.py
├── 11_final_report.py
└── 12_wgcna.R

Pipeline steps
Download RNA-seq data and reference files
Trim reads with Trimmomatic
Align reads with HISAT2
Quantify transcripts with StringTie
Merge abundance tables into TPM matrices
Run differential expression analysis
Generate volcano plots, heatmaps, and boxplots
Build co-expression networks
Run WGCNA module analysis
Annotate co-expressed genes with BLAST
Perform GO enrichment
Screen candidates against Rfam
Assess coding potential with CPC2
Predict miRNA interactions
Compare patterns across species
Quick start
git clone https://github.com/legoreddragon-ai/resurrection_lncrna_pipeline.git
cd resurrection_lncrna_pipeline

bash project_setup/00_preflight.sh
bash project_setup/01_download.sh
bash project_setup/02_trim.sh
bash project_setup/03_align.sh
bash project_setup/04_stringtie.sh

python3 project_setup/05_merge_matrix.py
python3 project_setup/06_de_analysis.py
python3 project_setup/07_visualize.py
python3 project_setup/08_coexpression.py
python3 project_setup/09_blast_annotate.py
python3 project_setup/10_go_enrichment.py
python3 project_setup/11_final_report.py
Rscript project_setup/12_wgcna.R

System requirements
Recommended environment:
Ubuntu 20.04+ or WSL2
16-32 GB RAM
150+ GB free disk
8+ CPU cores
Main dependencies
Python / command-line tools
Python 3.8+
Trimmomatic
HISAT2
SAMtools
StringTie
BLAST+
Infernal
R packages
DESeq2
WGCNA
Expected outputs
Representative output files include:
results/stringtie_tpm_matrix.csv
results/deseq2_desiccated_vs_hydrated.csv
results/coexpression_partners.csv
results/blast_annotated_partners.csv
results/go_enrichment_results.csv
results/wgcna_module_trait_correlation.csv
results/wgcna_candidate_modules.csv
results/final_report.md
Notes on analysis design
The C. plantagineum dataset has limited hydrated replication, which reduces power for genome-wide FDR-based discovery.
The 14 core candidates are the strongest DESeq2-supported set in this study.
The 661 exploratory candidates should be treated as hypothesis-generating rather than equally strong evidence.
For S. lepidophylla, no public GFF annotation was available in this workflow, so transcript discovery relied on de novo assembly and merge steps.
CPC2 classified 9 of 14 core candidates as confidently non-coding and flagged 5 of 14 as transcripts of interest that may include sORFs.
miRNA predictions were generated using proxy resources and should be treated as supportive computational evidence.
This repository contains a computational analysis pipeline only; experimental validation is still needed.
Species-specific notes
Craterostigma plantagineum
Primary desiccation-response lncRNA discovery was performed in C. plantagineum using public RNA-seq data from hydrated, desiccated, and recovery-associated states.
Selaginella lepidophylla
Cross-species comparison was performed in S. lepidophylla, including full time-course analysis. The time-course results suggested both early-responding and later-responding transcript subsets during desiccation progression.
Limitations
Limited replication in the source dataset affects statistical power.
S. lepidophylla analysis was performed without a public reference annotation.
CPC2 and miRNA predictions are computational and not definitive functional validation.
No wet-lab follow-up is included in this repository.
Citation
If you use this repository, please cite:
Sumanth N.
Discovery of Desiccation-Induced Long Non-Coding RNAs in the Resurrection Plant Craterostigma plantagineum: A Transcriptome-Wide Computational Analysis with Cross-Species Validation in Selaginella lepidophylla.
2026. Preprint.
Data sources
GEO: GSE157098
SRA: PRJNA660052
SRA: PRJNA420971



