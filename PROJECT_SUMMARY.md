# Project Summary

## Resurrection Plant lncRNA Discovery Pipeline

**Status: COMPLETE**
Last updated: March 2026
Manuscript: v6 final — bioRxiv submission ready

---

## Overview

First systematic transcriptome-wide characterisation of desiccation-induced long non-coding RNAs (lncRNAs) in any resurrection plant species. Primary analysis in *Craterostigma plantagineum* with cross-species validation in *Selaginella lepidophylla*.

**Author:** Neil Sumanth, Independent Researcher
**Repository:** github.com/legoreddragon-ai/resurrection_lncrna_pipeline
**Data:** GEO GSE157098 (Cp) / SRA PRJNA420971 (Sl)

---

## Core Findings

### 1. Unidirectional Activation Signature

In *C. plantagineum*, 661 lncRNAs were nominally upregulated during desiccation and none downregulated — a strictly one-sided activation programme that persists through rehydration (199 up, 0 down). This pattern is absent in *S. lepidophylla* (303 up / 98 down), making it a species-specific feature of the angiosperm resurrection response.

### 2. Fourteen High-Confidence Candidates

| Candidate | log2FC | DESeq2 padj | Des. TPM | CPC2 prob | Notes |
|-----------|--------|-------------|----------|-----------|-------|
| contig_551 | 13.38 | 4.5x10^-7 | 1,957 | 0.9995 | ABCC14-like in Sl; transcript of interest |
| contig_4950 | 12.35 | 1.1x10^-6 | 957 | 0.776 | TRAPPC6B-like in Sl; transcript of interest |
| contig_43263 | 12.30 | 5.9x10^-7 | 925 | 0.252 | Co-expresses with PP2C-60 |
| contig_2815 | 12.29 | 1.2x10^-6 | 922 | 0.700 | Top hub; miR164 sponge; ERD7; transcript of interest |
| contig_6094 | 12.17 | 1.2x10^-6 | 850 | 0.121 | Confirmed non-coding |
| contig_670 | 12.13 | 1.6x10^-6 | 823 | 0.023 | MYB73, NF-YC2 co-expression |
| contig_47679 | 12.02 | 1.2x10^-6 | 764 | 0.008 | Confirmed non-coding |
| contig_45133 | 12.02 | 1.6x10^-6 | 762 | 0.028 | NF-YC2 co-expression |
| contig_692 | 12.00 | 1.2x10^-6 | 755 | 0.661 | Transcript of interest |
| contig_3232 | 11.88 | 1.4x10^-6 | 693 | 0.093 | Confirmed non-coding |
| contig_24180 | 11.76 | 1.7x10^-6 | 638 | 0.017 | MYB73, NF-YC2 co-expression |
| contig_12238 | 11.74 | 1.4x10^-6 | 630 | 0.665 | 229 co-expression partners; transcript of interest |
| contig_15059 | 11.65 | 1.6x10^-6 | 591 | 0.059 | Confirmed non-coding |
| contig_8442 | 11.65 | 1.0x10^-6 | 590 | 0.037 | miR830-5p interaction |

All 14 show binary expression: 0 TPM hydrated, 590-1,957 TPM desiccated.
9/14 confirmed non-coding by CPC2 (prob < 0.5).
5/14 flagged as transcripts of interest (CPC2 prob > 0.5; may contain sORFs).
All 14 confirmed free of structural RNA contamination (Rfam, 0 hits).

### 3. Chloroplast Protection as Central Theme

Chloroplast and photosynthesis genes are the dominant functional category among co-expression partners in both species:
- *C. plantagineum*: 3.73-fold enriched (FDR = 6x10^-14)
- *S. lepidophylla*: 1.99-fold enriched (FDR = 3x10^-22)

This functional convergence across two phylogenetically distant lineages (angiosperm vs lycophyte, ~400 million years apart) strongly supports the hypothesis that lncRNA-mediated chloroplast protection is a convergently evolved feature of vegetative desiccation tolerance.

### 4. Mechanistic Divergence Between Species

While chloroplast protection is shared, the upstream regulatory architecture differs:
- *C. plantagineum*: ABA-independent; strictly unidirectional activation
- *S. lepidophylla*: ABA-dominant (67.44-fold enrichment); bidirectional response

This pattern of functional convergence with mechanistic divergence is consistent with independent evolutionary origins of desiccation tolerance in the two lineages.

### 5. contig_2815 as Top Mechanistic Candidate

Multiple converging lines of evidence identify contig_2815 as the highest-priority candidate for experimental follow-up:
- Top WGCNA hub gene in the MEred desiccation module (intramodular connectivity 231.76)
- Strongest miRNA sponge prediction: miR164a/b/c, expectation 2.5 (miR164 targets NAC TFs)
- BLAST annotation matches ERD7, a chloroplast-targeted early dehydration response protein
- Co-expresses with PP2C-60 (ABA signalling), WRKY3, MYB73, ABI5-like protein 3, TIC32

### 6. Lineage Specificity

Zero BLASTn hits against *Arabidopsis thaliana* TAIR10 cDNA at any threshold confirms that all 14 candidates are novel, resurrection-plant specific regulatory innovations with no conservation in drought-sensitive model plants.

### 7. Cross-Species Conservation (Partial)

3/14 candidates show partial sequence conservation in *S. lepidophylla* at relaxed BLAST threshold (E < 1e-3, 71-80% identity):
- contig_551 — matches ABC transporter ABCC14 (membrane transport)
- contig_22145 — matches Mitochondrial Rho GTPase MIRO1 (organelle dynamics)
- contig_4950 — matches TRAPPC6B (vesicle trafficking and membrane integrity)

All three conserved candidates are associated with membrane or organelle protection functions.

### 8. Temporal Structure in S. lepidophylla

Time-course analysis (0hr, 1hr, 6hr, 24hr, 120hr desiccation, rehydration) revealed two kinetically distinct subsets:
- Early responder (MSTRG.22821): active at 1hr, peaks at 24hr, sustained in rehydration
- Late responders (MSTRG.15112, 15115, 15111, 6238): silent at 1hr/6hr, activate at 24hr, sustained through 120hr and rehydration

---

## Analysis Summary

```
C. plantagineum transcripts:     48,045
lncRNA candidates identified:    38,978 (81.1%)
Expressed lncRNAs:               17,574
Core candidates (DESeq2):        14 (padj 4.5x10^-7 to 1.7x10^-6)
Confirmed non-coding (CPC2):     9/14
Transcripts of interest:         5/14
Exploratory set (Welch):         661 up, 0 down
Co-expression pairs (r>=0.90):   1,062
Annotated partner genes:         356
Rfam contamination hits:         0
miRNA interactions predicted:    106
Arabidopsis BLAST hits:          0

S. lepidophylla transcripts:     75,098 (de novo)
Sl expressed genes:              26,192
Sl desiccation upregulated:      303
Sl desiccation downregulated:    98
Sl chloroplast enrichment:       1.99x (FDR = 3x10^-22)
Sl ABA enrichment:               67.44x (FDR < 0.001)
Sl time-course samples:          17 (0hr, 1hr, 6hr, 24hr, 120hr, rehyd)
```

---

## Pipeline Components

### Data

- *C. plantagineum*: 10 samples (2 hydrated, 5 desiccated, 3 rehydrated) from GEO GSE157098
- *S. lepidophylla*: 17 samples (full time-course) from SRA PRJNA420971

### Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| Trimmomatic | 0.40 | Quality trimming |
| HISAT2 | 2.2.x | Alignment |
| SAMtools | 1.x | BAM processing |
| StringTie | 2.2.x | Transcript quantification |
| DESeq2 | 1.40 | Primary differential expression |
| WGCNA | 1.72 | Co-expression network analysis |
| BLASTx/n | 2.12+ | Functional annotation / conservation |
| Infernal cmscan | 1.1.x | Rfam structural RNA screening |
| CPC2 | 1.0.1 | Coding potential assessment |
| psRNATarget | 2017 v2 | miRNA sponge prediction |

### Key Results Files

| File | Description |
|------|-------------|
| `results/real14_candidates_stats.csv` | Core 14 candidates with expression stats |
| `results/deseq2_desiccated_vs_hydrated.csv` | Full DESeq2 results |
| `results/stringtie_tpm_matrix.csv` | C. plantagineum TPM matrix (48,045 x 10) |
| `results/coexpression_partners.csv` | All co-expression pairs (r >= 0.90) |
| `results/blast_annotated_partners.csv` | BLAST-annotated co-expression partners |
| `results/go_enrichment_results.csv` | Cp GO enrichment results |
| `results/wgcna_module_trait_correlation.csv` | WGCNA module-trait correlations |
| `results/wgcna_hub_genes_grey60.csv` | Grey60 module hub genes |
| `results/wgcna_hub_genes_red.csv` | Red module hub genes |
| `results/mirna_sponge/psrnatarget_results.txt` | psRNATarget predictions |
| `results/coding_potential/real14_cpc2_summary.txt` | CPC2 results for all 14 candidates |
| `results/selaginella_tpm_matrix.csv` | Sl 9-sample TPM matrix |
| `results/selaginella_timecourse_tpm_matrix.csv` | Sl 17-sample time-course matrix |
| `results/selaginella_go_enrichment.csv` | Sl GO enrichment results |
| `results/selaginella_comparison_summary.txt` | Cross-species comparison summary |
| `results/conserved_candidates_summary.txt` | Annotation of 3 conserved candidates |
| `results/aba_pathway_partners.csv` | PP2C/ABA pathway co-expression |
| `results/transcription_factor_coexpression.csv` | TF co-expression partners |
| `results/figures/selaginella_heatmap.png` | Sl top candidates heatmap |
| `results/figures/selaginella_timecourse.png` | Sl time-course kinetics plot |

---

## Manuscript

**Title:** Discovery of Desiccation-Induced Long Non-Coding RNAs in the Resurrection Plant *Craterostigma plantagineum*: A Transcriptome-Wide Computational Analysis with Cross-Species Validation in *Selaginella lepidophylla*

**Structure:** Abstract (Background/Methods/Results/Conclusions) / Introduction / Methods (9 subsections) / Results (10 subsections) / Discussion (7 subsections) / Conclusions / References (22)

**Target:** bioRxiv preprint (immediate) then PLOS ONE / BMC Plant Biology / Frontiers in Plant Science

**Version history:**
- v1: Initial manuscript, Welch t-test, wrong 14 candidates
- v2: Added Selaginella, WGCNA, miRNA sponge, cross-species
- v3: Added Selaginella time-course kinetics
- v4: DESeq2 as primary DE method, corrected 14 candidates, fixed CPC2 classification, corrected references
- v5: All reviewer fixes — language softening, functional convergence framing, ref 22 corrected
- v6 final: sORF note, CPC2 numeric probabilities in Table 1, chloroplast RBP mechanism sentence

---

## Known Limitations

- n=2 hydrated replicates (Cp) limits genome-wide statistical power
- De novo assembly (Sl) may introduce transcript boundary artefacts
- No wet-lab validation (qPCR, RIP-seq, degradome sequencing)
- Arabidopsis miRNA proxy used for psRNATarget — resurrection plant miRNA databases do not exist
- 5/14 candidates flagged by CPC2 — experimental confirmation of lncRNA status required

---

## Priority Future Directions

1. qPCR validation of contig_8442, contig_15059, contig_43263 (strongest DESeq2 p-values)
2. RNA pull-down or RIP-seq for contig_2815 to test miR164 sequestration
3. Chromatin immunoprecipitation to assess lncRNA recruitment to chloroplast gene loci
4. Realignment to VanBuren 2023 high-quality octoploid *C. plantagineum* genome (doi:10.1111/tpj.16165)
5. Resurrection plant-specific miRNA database construction via small RNA-seq
6. Single-molecule long-read sequencing to resolve full-length lncRNA isoforms

---

*Project: Resurrection Plant lncRNA Discovery*
*Last updated: March 2026*
*Status: Complete — manuscript v6 final, bioRxiv submission pending*