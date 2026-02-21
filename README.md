# Discovery and Characterization of Desiccation-Induced lncRNAs in Resurrection Plants
# This is work in progress.  

## Overview

This project reanalyzes public RNA-seq data from *Craterostigma plantagineum* (a resurrection plant) to identify and characterize stress-responsive long non-coding RNAs (lncRNAs) during desiccation cycles.

### Key Findings

14 stress-responsive lncRNAs identified, all downregulated during dehydration:
- Top candidate: Cp_V2_contig_159 (STRG.167) - 7.8-fold reduction during desiccation
- Expression levels: Hydrated (1,320 TPM) → Dehydrated (5 TPM) → Rehydrated (13 TPM)
- Suggests "hydration-maintenance" regulatory role

## Project Status

Current Progress: 75% Complete

### Completed Steps
1. Literature Review - Hypothesis, research questions, methods documented
2. Data Acquisition - 60 SRA accessions (PRJNA660052) downloaded and organized
3. RNA-Seq Processing - Quality control, trimming, alignment (HISAT2), assembly (StringTie)
4. lncRNA Identification - 48,045 transcripts filtered; 38,978 lncRNAs classified
5. Differential Expression - 14 significant stress-responsive lncRNAs identified (padj<0.05, |log2FC|>1)

### In Progress / TODO
6. Co-expression Network - Identify genes clustering with the 14 lncRNAs
7. Functional Enrichment - GO/KEGG analysis (LEA, ABA, antioxidant pathways)
8. Visualization & Validation - Heatmaps, networks, BLAST novelty assessment
9. Documentation - Final results tables, pipeline documentation

## Data & Results

### Input Data
- Species: *Craterostigma plantagineum* (resurrection plant)
- Samples: 10 biological replicates across 3 conditions
  - Hydrated: 2 replicates
  - Rehydrated: 3 replicates
  - Dehydrated (water content 2-60%): 5 replicates
- Sequencing: Paired-end RNA-seq

### Key Results Files

File | Size | Description
-----|------|------------
results/significant_stress_responsive_lncrnas.csv | 1.2K | 14 key lncRNAs with DE statistics
results/lncrna_differential_expression.csv | 3.9M | All 31,267 lncRNAs with expression stats
results/stringtie_tpm_matrix.csv | 4.0M | Expression matrix (42,191 genes x 10 samples)
results/coding_potential.txt | 1.9M | ORF-based coding classification
references/annotation.gff3 | 8.3M | Gene annotation (GFF3 format)

## Methods

### RNA-Seq Pipeline
1. Quality Control: Trimmomatic adapter & quality trimming
2. Alignment: HISAT2 to transcriptome reference
3. Assembly: StringTie transcript assembly with abundance quantification
4. lncRNA Classification: ORF-based coding potential scoring (threshold: 50%)
5. Differential Expression: t-tests (Benjamini-Hochberg FDR correction)
   - Filters: adjusted p-value < 0.05, |log2 fold-change| > 1

## Installation & Usage

### Requirements
```bash
conda install -c bioconda trimmomatic hisat2 samtools stringtie
pip install pandas numpy scipy biopython
```

### Quick Start
```bash
# Run full pipeline
bash pipeline.sh
```

## Project Structure
```
resurrection_lncrna_pipeline/
├── README.md
├── LICENSE
├── pipeline.sh
├── data/
│   └── fastq_trimmed/     (Quality-trimmed FASTQ files)
├── references/
│   ├── annotation.gff3    (Gene annotation)
│   ├── annotation.gtf     (Gene annotation, GTF format)
│   └── GSE157098_Cp_transcriptome_assembly_V2.fa
├── results/
│   ├── significant_stress_responsive_lncrnas.csv  (KEY RESULTS)
│   ├── lncrna_differential_expression.csv
│   ├── stringtie_tpm_matrix.csv
│   ├── coding_potential.txt
│   └── stringtie_output/
│       └── *_abundance.txt
└── scripts/
    └── (Analysis scripts)
```

## Key Statistics

- Total transcripts: 48,045
- Classified as lncRNA: 38,978 (81%)
- Classified as coding: 9,067 (19%)
- Significant stress-responsive lncRNAs: 14 (padj<0.05, |log2FC|>1)
- Expression range: Up to 1,320 TPM (hydrated) to <1 TPM (dehydrated)

## Future Directions

1. Network Analysis - Co-expression modules with coding genes
2. Functional Enrichment - GO/KEGG pathways (LEA, ABA signaling, antioxidants)
3. Target Prediction - miRNA sponge & cis-regulatory roles
4. Cross-species Comparison - Selaginella spp. & Arabidopsis lncRNA databases
5. Wet-lab Validation - qPCR, CRISPR knockdown (future work)

## License

This project is open source under the MIT License.

## Contact

For questions or collaboration, please reach out via GitHub Issues.
