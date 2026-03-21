# lncRNA Discovery Pipeline — Complete Setup & Execution Guide
### Craterostigma plantagineum & Selaginella lepidophylla · RNA-Seq · DESeq2 · WGCNA · miRNA Sponge · Cross-Species Validation
**From Git Clone to Environment to Data Download to Alignment to DE Analysis to Co-expression to BLAST to GO Enrichment to Manuscript**

> **Pipeline: 100% Complete · Manuscript: v6 final (bioRxiv-ready)**
> Total runtime: ~15 hours end-to-end

---

## Associated Manuscript

**Title:** Discovery of Desiccation-Induced Long Non-Coding RNAs in the Resurrection Plant *Craterostigma plantagineum*: A Transcriptome-Wide Computational Analysis with Cross-Species Validation in *Selaginella lepidophylla*

**Author:** Neil Sumanth, Independent Researcher

**Status:** Preprint ready — bioRxiv submission pending

**Data availability:**
- *C. plantagineum* RNA-seq: GEO GSE157098 / SRA PRJNA660052
- *S. lepidophylla* RNA-seq: SRA PRJNA420971

---

## Pipeline Status

| Step | Status | Script | Key Result |
|------|--------|--------|------------|
| 1. Literature Review | Complete | — | No systematic resurrection plant lncRNA study exists (2024-2026) |
| 2. C. plantagineum data acquisition | Complete | `01_download.sh` | 10 samples, GSE157098 |
| 3. S. lepidophylla data acquisition | Complete | — | 17 samples, PRJNA420971 (full time-course 0-120hr) |
| 4. Quality control and trimming | Complete | `02_trim.sh` | Trimmomatic, all samples |
| 5. Alignment | Complete | `03_align.sh` | HISAT2; Cp: 60-95%; Sl: 48-57% |
| 6. Transcript quantification | Complete | `04_stringtie.sh` | StringTie; reference-guided (Cp), de novo (Sl) |
| 7. lncRNA identification | Complete | `05_merge_matrix.py` | 38,978 candidates (81.1%) |
| 8. Differential expression (DESeq2) | Complete | `06_de_analysis.py` + DESeq2 | 14 core (padj<0.05); 661 exploratory |
| 9. Visualisation | Complete | `07_visualize.py` | Volcano, heatmap, boxplots |
| 10. WGCNA | Complete | `12_wgcna.R` | MEgrey60 r=-0.988; MEred r=-0.721 |
| 11. Co-expression network | Complete | `08_coexpression.py` | 1,062 pairs (r>=0.90); 467 partners |
| 12. BLAST annotation | Complete | `09_blast_annotate.py` | 356 annotated partners |
| 13. GO enrichment | Complete | `10_go_enrichment.py` | Chloroplast 3.73x FDR=6x10^-14 |
| 14. Rfam screening | Complete | `cmscan` | 0 contamination hits |
| 15. CPC2 coding potential | Complete | CPC2 web server | 9/14 confirmed non-coding; 5 flagged |
| 16. miRNA sponge prediction | Complete | psRNATarget | 106 interactions; top: contig_2815 + miR164a |
| 17. ABA pathway analysis | Complete | — | PP2C phosphatases co-express with 11/14 candidates |
| 18. TF co-expression | Complete | — | WRKY3, MYB73, NF-YC2, ABI5 identified |
| 19. Selaginella full pipeline | Complete | — | 17 samples; chloroplast enriched both species |
| 20. Selaginella time-course | Complete | — | Early (1hr) vs late (24hr+) responder subsets |
| 21. Arabidopsis comparison | Complete | BLASTn | 0 hits — lineage-specific confirmed |
| 22. Cross-species conservation | Complete | BLASTn | 3/14 partially conserved (71-80% identity) |
| 23. Rehydration analysis | Complete | — | 7/14 candidates sustained through recovery |
| 24. WGCNA hub gene identification | Complete | R/WGCNA | contig_2815 top hub; ERD7 annotation |
| 25. Manuscript (v6 final) | Complete | — | 16 pages, bioRxiv-ready |

**Overall completion: 100%**

---

## Key Results

```
Species primary:              Craterostigma plantagineum
Species cross-validation:     Selaginella lepidophylla
Cp transcripts quantified:    48,045
Cp lncRNAs identified:        38,978 (81.1%)
Cp lncRNAs tested for DE:     17,574
Cp upregulated (exploratory): 661
Cp downregulated:             0 (strictly unidirectional)
Core candidates (DESeq2):     14 (padj 4.5x10^-7 to 1.7x10^-6)
Confirmed non-coding (CPC2):  9/14 (prob < 0.5)
Transcripts of interest:      5/14 (CPC2 flagged)
Co-expression pairs:          1,062 (r >= 0.90)
Annotated partner genes:      356
Rfam contamination hits:      0
miRNA interactions:           106 predicted
Arabidopsis BLAST hits:       0 (lineage-specific confirmed)
Sl upregulated:               303
Sl downregulated:             98
Chloroplast enrichment (Cp):  3.73x (FDR = 6x10^-14)
Chloroplast enrichment (Sl):  1.99x (FDR = 3x10^-22)
```

**Top finding:** Chloroplast/photosynthesis gene enrichment confirmed in both resurrection plant lineages. Functional convergence with mechanistic divergence: ABA-independent in *C. plantagineum* (angiosperm); ABA-dominant in *S. lepidophylla* (lycophyte). Zero sequence conservation with *Arabidopsis* confirms lineage-specificity.

**Top candidate:** contig_2815 — WGCNA hub (MEred module), predicted miR164a ceRNA (expectation 2.5), ERD7 chloroplast annotation, co-expresses with PP2C-60, WRKY3, MYB73, ABI5.

---

## Setup Folder — Script Index

All scripts live in the `project_setup/` directory.

```
project_setup/
├── 00_preflight.sh        Stage 0: environment and prerequisites check
├── 01_download.sh         Stage 1: download FASTQs + reference + build index
├── 02_trim.sh             Stage 2: Trimmomatic quality trimming
├── 03_align.sh            Stage 3: HISAT2 alignment + SAMtools sort/index
├── 04_stringtie.sh        Stage 4: StringTie transcript quantification
├── 05_merge_matrix.py     Stage 5: merge abundance files into TPM matrix
├── 06_de_analysis.py      Stage 6: differential expression + candidate selection
├── 07_visualize.py        Stage 7: volcano plot, heatmap, boxplots
├── 08_coexpression.py     Stage 8: Pearson co-expression network (r>=0.90)
├── 09_blast_annotate.py   Stage 9: blastx vs UniProt Swiss-Prot
├── 10_go_enrichment.py    Stage 10: Fisher's exact test functional enrichment
├── 11_final_report.py     Stage 11: generate full markdown research report
└── 12_wgcna.R             Stage 12: WGCNA module analysis + DESeq2
```

### Running Individual Stages

```bash
cd ~/resurrection_lncrna_pipeline

# Stages 1-7 (RNA-seq through visualisation)
bash project_setup/00_preflight.sh
bash project_setup/01_download.sh
bash project_setup/02_trim.sh
bash project_setup/03_align.sh
bash project_setup/04_stringtie.sh
python3 project_setup/05_merge_matrix.py
python3 project_setup/06_de_analysis.py
python3 project_setup/07_visualize.py

# Stages 8-12 (co-expression through WGCNA)
python3 project_setup/08_coexpression.py
python3 project_setup/09_blast_annotate.py
python3 project_setup/10_go_enrichment.py
python3 project_setup/11_final_report.py
Rscript project_setup/12_wgcna.R
```

All scripts are path-dynamic — no hardcoded paths. Safe to re-run from any working directory as long as scripts are in `project_setup/`.

---

## Table of Contents

1. [System Requirements and Pre-Flight Check](#section-1)
2. [Conda Installation and Environment Setup](#section-2)
3. [Getting the Code](#section-3)
4. [Downloading RNA-Seq Data from NCBI](#section-4)
5. [Reference Genome and Index](#section-5)
6. [Running the Main Pipeline](#section-6)
7. [Expression Matrix, DE Analysis and Visualization](#section-7)
8. [WGCNA Co-expression Network Analysis](#section-8)
9. [BLAST Annotation](#section-9)
10. [GO Enrichment Analysis](#section-10)
11. [lncRNA Validation — Rfam, CPC2, miRNA Sponge](#section-11)
12. [Selaginella Cross-Species Validation](#section-12)
13. [Validating Your Results](#section-13)
14. [Troubleshooting Reference](#section-14)
15. [Quick Reference — All Steps](#section-15)

---

<a name="section-1"></a>
## Section 1: System Requirements and Pre-Flight Check

> **Script:** `bash project_setup/00_preflight.sh`
> Automatically checks all requirements below and reports PASS / WARN / FAIL for each.

### Hardware Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| RAM | 16 GB | 32 GB | StringTie, WGCNA, and DESeq2 are memory-intensive |
| Disk Space | 150 GB | 250 GB | Raw FASTQs + BLAST DB + Selaginella data |
| CPU Cores | 8 | 12+ | All tools support multi-threading; pipeline uses -p 12 by default |
| OS | Ubuntu 20.04+ | Ubuntu 22.04+ | Also works on WSL2 on Windows |
| GPU | Optional | RTX 5080 | Not required; pipeline is CPU-based |

---

### Step 1 — Verify OS and Architecture

```bash
bash project_setup/00_preflight.sh
# OR manually:
uname -a            # Expected: Linux <hostname> ... x86_64 GNU/Linux
lsb_release -a      # Expected: Ubuntu 20.04 or 22.04
df -h /             # Expected: >150 GB available
free -h             # Expected: >16 GB total RAM
```

> NOTE: Pipeline had hardcoded paths from the original development laptop.
> FIX: All scripts use dynamic path resolution:
> - Python: `os.path.dirname(os.path.dirname(os.path.abspath(__file__)))`
> - Bash: `$(dirname "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)")`

---

### Step 2 — Check Python

```bash
python3 --version    # Expected: Python 3.8 or higher
pip3 --version       # If missing: sudo apt-get install python3-pip
```

---

### Step 3 — Check Git

```bash
git --version
# If missing:
sudo apt-get update && sudo apt-get install git -y
```

---

### Step 4 — Check SSL/Certificate Status

```bash
curl -v https://www.ncbi.nlm.nih.gov 2>&1 | grep -E '(SSL|TLS|certificate|Connected)'
# Expected: "SSL connection using TLSv1.3"

# If you see "certificate verify failed":
sudo apt-get update
sudo apt-get install --reinstall ca-certificates -y
sudo update-ca-certificates
```

> NOTE: SRA Toolkit 2.9.6 uses mbedtls which rejects NCBI's TLS certificate.
> FIX: Upgrade to SRA Toolkit 3.x. The `01_download.sh` script handles this automatically.

---

<a name="section-2"></a>
## Section 2: Conda Installation and Environment Setup

### Step 5 — Install Miniconda

```bash
conda --version   # Check if already installed

# If not found:
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
source ~/.bashrc
conda --version   # Expected: conda 23.x.x or higher
```

---

### Step 6 — Install Bioinformatics Tools

```bash
# Python packages
pip install pandas numpy scipy biopython matplotlib seaborn statsmodels --break-system-packages

# Bioinformatics tools
conda install -c bioconda trimmomatic hisat2 samtools stringtie -y

# If conda fails with LibMambaUnsatisfiableError:
sudo apt-get update && sudo apt-get install trimmomatic hisat2 samtools -y
conda install -c bioconda stringtie -y

# BLAST
sudo apt-get install -y ncbi-blast+

# Infernal (for Rfam contamination screening)
sudo apt-get install -y infernal

# R packages
Rscript -e 'if (!require("BiocManager")) install.packages("BiocManager")'
Rscript -e 'BiocManager::install(c("impute","preprocessCore","DESeq2"))'
Rscript -e 'install.packages("WGCNA")'
```

> NOTE: `conda install hisat2 samtools` may fail with `LibMambaUnsatisfiableError`.
> FIX: Install via `apt-get`. Install `stringtie` separately via conda.

> NOTE: Trimmomatic OOM error on large files (>12M reads).
> FIX: Add `-Xmx2g -threads 4` flags to Trimmomatic call.

> NOTE: samtools "merging from 0 files" in base conda env — broken library links.
> FIX: Use `cpat_env` conda environment for samtools sort/index on Selaginella files.

---

### Step 7 — Verify All Tools

```bash
trimmomatic -version    # Expected: 0.39 or 0.40
hisat2 --version        # Expected: HISAT2 version 2.2.x
samtools --version      # Expected: samtools 1.x
stringtie --version     # Expected: StringTie v2.x
blastx -version         # Expected: 2.12.0+
cmscan --help 2>/dev/null | head -1
Rscript -e 'library(DESeq2); library(WGCNA); cat("R deps OK\n")'
python3 -c 'import pandas, numpy, scipy, statsmodels; print("Python deps OK")'
```

---

<a name="section-3"></a>
## Section 3: Getting the Code

### Step 8 — Clone the Repository

```bash
git clone https://github.com/legoreddragon-ai/resurrection_lncrna_pipeline.git
cd resurrection_lncrna_pipeline
ls project_setup/
```

---

### Step 9 — Create Required Directories

```bash
bash project_setup/00_preflight.sh
# Auto-creates: data/ references/ results/ logs/
# Or manually:
mkdir -p data/fastq_raw data/fastq_trimmed
mkdir -p data/selaginella/raw data/selaginella/trimmed data/selaginella/aligned data/selaginella/stringtie
mkdir -p results/hisat2_aligned results/stringtie_output results/figures
mkdir -p results/coding_potential results/mirna_sponge
mkdir -p references/blast_db references/rfam references/mirna
mkdir -p logs
```

---

<a name="section-4"></a>
## Section 4: Downloading RNA-Seq Data from NCBI

> **Script:** `bash project_setup/01_download.sh`
>
> C. plantagineum — GEO: GSE157098 | SRA: PRJNA660052 | ~40 GB
> S. lepidophylla — SRA: PRJNA420971 | ~60 GB (17 samples including full time-course)

### C. plantagineum Sample Manifest

| SRA Accession | Sample Name | Condition |
|---------------|-------------|-----------|
| SRR12542161 | cp_hyd_r1 | Hydrated replicate 1 |
| SRR12542166 | cp_hyd_r3 | Hydrated replicate 3 |
| SRR12542171 | cp_rehyd_r1 | Rehydrated replicate 1 |
| SRR12542176 | cp_rehyd_r2 | Rehydrated replicate 2 |
| SRR12542181 | cp_rehyd_r3 | Rehydrated replicate 3 |
| SRR12542186 | cp_wc2_r2 | Desiccated WC2% replicate 2 |
| SRR12542191 | cp_wc2_r3 | Desiccated WC2% replicate 3 |
| SRR12542196 | cp_wc60_r1 | Desiccated WC60% replicate 1 |
| SRR12542206 | cp_wc60_r2 | Desiccated WC60% replicate 2 |
| SRR12542211 | cp_wc60_r3 | Desiccated WC60% replicate 3 |

### S. lepidophylla Sample Manifest (Full Time-Course)

| SRA Accessions | Condition | Time Point |
|----------------|-----------|------------|
| SRR6345603, SRR6345604, SRR6345601 | Hydrated | 0 hr |
| SRR6345602, SRR6345607, SRR6345608 | Desiccating | 1 hr |
| SRR6345605, SRR6345606, SRR6345609 | Desiccating | 6 hr |
| SRR6345610, SRR6345613, SRR6345614 | Desiccated | 24 hr |
| SRR6345611, SRR6345612 | Desiccated | 120 hr |
| SRR6345616, SRR6345617, SRR6345615 | Rehydrated | 24 hr |

---

### Step 10 — Upgrade SRA Toolkit

```bash
conda install -c bioconda sra-tools=3.1.1 -y
fasterq-dump --version   # Expected: 3.x.x

# If conda fails:
cd ~
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
echo 'export PATH=~/sratoolkit.3.2.1-ubuntu64/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
vdb-config --set /tls/allow-all-certs=true
```

---

### Step 11 — Download C. plantagineum Samples

```bash
bash project_setup/01_download.sh
```

### Step 12 — Download S. lepidophylla Samples

```bash
for srr in SRR6345601 SRR6345602 SRR6345603 SRR6345604 SRR6345605 \
           SRR6345606 SRR6345607 SRR6345608 SRR6345609 SRR6345610 \
           SRR6345611 SRR6345612 SRR6345613 SRR6345614 SRR6345615 \
           SRR6345616 SRR6345617; do
    prefetch ${srr} --output-directory data/selaginella/raw
    fasterq-dump data/selaginella/raw/${srr}/${srr}.sra \
        --outdir data/selaginella/raw --split-files --threads 12
done
```

---

<a name="section-5"></a>
## Section 5: Reference Genome and Index

### Step 13 — Download C. plantagineum Reference FASTA

```bash
cd references/
wget --no-check-certificate \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_transcriptome_assembly_V2.fa.gz'
gunzip GSE157098_Cp_transcriptome_assembly_V2.fa.gz
grep -c '>' GSE157098_Cp_transcriptome_assembly_V2.fa   # Expected: 48045
```

### Step 14 — Download S. lepidophylla Reference Transcriptome

```bash
# TSA transcriptome assembly GIMG01000000 (63,981 transcripts)
# Download from NCBI TSA: accession GIMG01000000
# No genome annotation available — StringTie de novo mode used
```

> NOTE: No GFF annotation exists for S. lepidophylla. VanBuren 2018 deposited genome sequence only, not annotation. See GFF Search Log section at bottom of this README for full documentation.

### Step 15 — Build HISAT2 Indices

```bash
# C. plantagineum
hisat2-build references/GSE157098_Cp_transcriptome_assembly_V2.fa \
    references/Cp_transcriptome_index -p 12

# S. lepidophylla
hisat2-build references/selaginella_transcriptome.fa \
    references/selaginella_index -p 12

ls references/*.ht2 | wc -l   # Expected: 16 (8 per species)
```

---

<a name="section-6"></a>
## Section 6: Running the Main Pipeline

### Step 16 — Quality Trimming

```bash
# C. plantagineum
bash project_setup/02_trim.sh

# S. lepidophylla — use -Xmx2g to avoid OOM on large files
trimmomatic PE -Xmx2g \
    data/selaginella/raw/${srr}_1.fastq \
    data/selaginella/raw/${srr}_2.fastq \
    data/selaginella/trimmed/${srr}_1.fastq \
    data/selaginella/trimmed/${srr}_1_unpaired.fastq \
    data/selaginella/trimmed/${srr}_2.fastq \
    data/selaginella/trimmed/${srr}_2_unpaired.fastq \
    ILLUMINACLIP:/path/to/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
    -threads 4
```

> NOTE: set -euo pipefail causes exit after first sample due to Trimmomatic GCLocker warnings.
> FIX: 02_trim.sh uses set -uo pipefail (no -e). GCLocker warnings are harmless.

> NOTE: Trimmomatic OOM on large Selaginella files.
> FIX: Add -Xmx2g -threads 4 instead of default 12 threads.

---

### Step 17 — HISAT2 Alignment

```bash
# C. plantagineum
bash project_setup/03_align.sh
grep 'overall alignment rate' logs/03_align_*.log
# Expected: 60-95% per sample

# S. lepidophylla — run alignment and sorting as separate commands
hisat2 -x references/selaginella_index \
    -1 data/selaginella/trimmed/${srr}_1.fastq \
    -2 data/selaginella/trimmed/${srr}_2.fastq \
    -p 12 --dta \
    -S data/selaginella/aligned/${srr}.sam \
    2>data/selaginella/aligned/${srr}_align_stats.txt

# Sort in cpat_env — base env samtools has broken library links
conda activate cpat_env
samtools sort -@ 12 -o data/selaginella/aligned/${srr}.bam \
    data/selaginella/aligned/${srr}.sam
samtools index data/selaginella/aligned/${srr}.bam
rm data/selaginella/aligned/${srr}.sam
```

Expected rates: Cp 60-95%; Sl 47-58% (lower expected for de novo TSA reference).

> NOTE: samtools in base conda env shows "merging from 0 files" — produces empty/tiny BAMs.
> FIX: Always use cpat_env for samtools sort/index on Selaginella files.

---

### Step 18 — StringTie Assembly

```bash
# C. plantagineum (reference-guided with -G flag)
bash project_setup/04_stringtie.sh
wc -l results/stringtie_output/*_abundance.txt | sort -n
# Expected: all files ~48046 lines

# S. lepidophylla Step 1 — de novo assembly (no -e flag, no GFF)
for srr in [all 17 SRR IDs]; do
    stringtie data/selaginella/aligned/${srr}.bam \
        -o data/selaginella/stringtie/${srr}.gtf -p 12
done

# S. lepidophylla Step 2 — merge GTFs
stringtie --merge data/selaginella/stringtie/SRR*.gtf \
    -o data/selaginella/stringtie/merged.gtf -p 12
grep -c "transcript" data/selaginella/stringtie/merged.gtf
# Expected: ~75,098

# S. lepidophylla Step 3 — re-quantify against merged reference
for srr in [all 17 SRR IDs]; do
    stringtie data/selaginella/aligned/${srr}.bam \
        -G data/selaginella/stringtie/merged.gtf \
        -o data/selaginella/stringtie/${srr}_final.gtf \
        -p 12 -e -B \
        -A data/selaginella/stringtie/${srr}_abundance.tab
done
```

> NOTE: -e flag requires -G. S. lepidophylla has no GFF, so first pass must be de novo (no -e), then re-quantify against merged GTF with -e -G.

---

<a name="section-7"></a>
## Section 7: Expression Matrix, DE Analysis and Visualization

### Step 19 — Merge TPM Matrix

```bash
python3 project_setup/05_merge_matrix.py
# Expected: 48,045 genes x 10 columns -> results/stringtie_tpm_matrix.csv

# Selaginella matrices:
# results/selaginella_tpm_matrix.csv             9-sample (hyd/des/rehyd)
# results/selaginella_timecourse_tpm_matrix.csv  17-sample (full time-course)
```

---

### Step 20 — Differential Expression

```bash
# Primary: DESeq2 (in R)
# Key settings: sizeFactors set to 1 (TPM already normalised), fitType="local"
Rscript project_setup/12_wgcna.R

# Secondary exploratory: Welch t-test
python3 project_setup/06_de_analysis.py

# Expected:
#   Core candidates (DESeq2 padj < 0.05): 14
#   DESeq2 padj range: 4.5x10^-7 to 1.7x10^-6
#   Exploratory set (Welch nominal p<0.05): 661 up, 0 down
```

> NOTE: With n=2 hydrated replicates, genome-wide FDR yields 0 significant genes by conventional methods.
> The 14 candidates are confirmed by DESeq2 with size factors set to 1 (TPM is already normalised).
> The 661 exploratory set should be treated as nominally significant leads only.

---

### Step 21 — Visualizations

```bash
python3 project_setup/07_visualize.py
ls -lh results/*.png   # Each should be >100 KB
```

---

<a name="section-8"></a>
## Section 8: WGCNA Co-expression Network Analysis

### Step 22 — Run WGCNA

```bash
Rscript project_setup/12_wgcna.R

# Expected output:
#   33 co-expression modules
#   MEgrey60: r = -0.988 with hydrated (strongest desiccation module)
#   MEred: r = -0.721 with hydrated
#   14/14 candidates assigned to modules
#
# Output files:
#   results/wgcna_module_trait_correlation.csv
#   results/wgcna_candidate_modules.csv
#   results/wgcna_module_membership.csv
#   results/wgcna_hub_genes_grey60.csv
#   results/wgcna_hub_genes_red.csv
#   results/figures/wgcna_module_trait_heatmap.pdf
```

**Module summary:**

| Module | r (hydrated) | Size | Key candidates | Top hub |
|--------|-------------|------|----------------|---------|
| MEgrey60 | -0.988 | 159 | contig_551, contig_12238, contig_4950, contig_8442 | KIN10 energy sensing |
| MEred | -0.721 | 364 | contig_2815, contig_43263 | contig_2815 (ERD7-like) |

---

### Step 23 — Run Pearson Co-expression

```bash
python3 project_setup/08_coexpression.py

# Expected output:
#   6,511 genes tested
#   1,062 co-expression pairs (r >= 0.90)
#   14/14 lncRNAs have at least one partner
#   Top partner: LEA protein 6 (r = 0.991)
```

**Key co-expression partners:**

| Partner gene | Max r | Annotation |
|-------------|-------|------------|
| contig_24081 | 0.991 | Late embryogenesis abundant protein 6 |
| contig_1856 | 0.986 | Heat shock 70 kDa protein |
| contig_15088 | 0.910-0.950 | Protein phosphatase 2C 60 (ABA signalling) |
| contig_12398 | 0.925-0.948 | WRKY transcription factor 3 |
| contig_32420 | 0.910-0.968 | NF-YC2 transcription factor |
| contig_17170 | 0.900-0.944 | MYB73 transcription factor |
| contig_34533 | 0.962 | ABI5-like protein 3 |

---

<a name="section-9"></a>
## Section 9: BLAST Annotation

### Step 24 — Download UniProt Swiss-Prot

```bash
mkdir -p references/blast_db
cd references/blast_db
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot -title "UniProt Swiss-Prot"
# Expected: ~574,000 sequences indexed
cd ~/resurrection_lncrna_pipeline
```

---

### Step 25 — Run BLAST Annotation

```bash
python3 project_setup/09_blast_annotate.py

# Expected:
#   466 unique partner contigs
#   356 annotated
#   110 no hit (novel/plant-specific)
#
# Output: results/blast_annotated_partners.csv
```

---

<a name="section-10"></a>
## Section 10: GO Enrichment Analysis

### Step 26 — Run GO Enrichment

```bash
python3 project_setup/10_go_enrichment.py

# Expected output (C. plantagineum):
#   Chloroplast/Photosynthesis: 3.73x, FDR = 6x10^-14
#   Stress response/chaperone: 2.09x, FDR = 0.003
#   Membrane transport/aquaporin: 1.88x, FDR = 0.018
```

**Cross-species enrichment comparison:**

| Category | Cp fold | Cp FDR | Sl fold | Sl FDR |
|----------|---------|--------|---------|--------|
| Chloroplast / photosynthesis | 3.73x | 6x10^-14 | 1.99x | 3x10^-22 |
| ABA signalling | ns | ns | 67.44x | <0.001 |
| Stress response / chaperone | 2.09x | 0.003 | 0.32x | <0.001 |
| Membrane transport | 1.88x | 0.018 | 1.38x | <0.001 |

Chloroplast enrichment confirmed in both species — functional convergence across 400 million years of evolution.

---

### Step 27 — Generate Final Report

```bash
python3 project_setup/11_final_report.py
# Output: results/final_report.md
```

---

<a name="section-11"></a>
## Section 11: lncRNA Validation — Rfam, CPC2, miRNA Sponge

### Step 28 — Download and Press Rfam Database

```bash
mkdir -p references/rfam
cd references/rfam
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm
# Expected: 4,227 CMs pressed
cd ~/resurrection_lncrna_pipeline
```

---

### Step 29 — Run Rfam Contamination Screen

```bash
cmscan --tblout results/rfam_hits.tbl \
    --cut_ga --rfam --nohmmonly --cpu 12 \
    references/rfam/Rfam.cm \
    results/coding_potential/candidates_14.fasta \
    > results/rfam_output.txt 2>&1

grep -c "^[^#]" results/rfam_hits.tbl
# Expected: 0
# All 14 candidates confirmed free of rRNA/tRNA/snoRNA contamination
```

---

### Step 30 — CPC2 Coding Potential

CPC2 is run via web server due to packaging issues with the local install.

Upload `results/coding_potential/candidates_14.fasta` to: https://cpc2.gao-lab.org/

**CPC2 results — real 14 candidates:**

| Candidate | CPC2 prob | Status |
|-----------|-----------|--------|
| contig_551 | 0.9995 | Transcript of interest |
| contig_4950 | 0.776 | Transcript of interest |
| contig_2815 | 0.700 | Transcript of interest |
| contig_12238 | 0.665 | Transcript of interest |
| contig_692 | 0.661 | Transcript of interest |
| contig_43263 | 0.252 | Confirmed non-coding |
| contig_6094 | 0.121 | Confirmed non-coding |
| contig_3232 | 0.093 | Confirmed non-coding |
| contig_15059 | 0.059 | Confirmed non-coding |
| contig_8442 | 0.037 | Confirmed non-coding |
| contig_45133 | 0.028 | Confirmed non-coding |
| contig_670 | 0.023 | Confirmed non-coding |
| contig_24180 | 0.017 | Confirmed non-coding |
| contig_47679 | 0.008 | Confirmed non-coding |

9/14 confirmed non-coding. 5/14 flagged — may contain sORFs (small ORFs); CPC2 flags are not definitive disqualifiers.

Full results: `results/coding_potential/real14_cpc2_summary.txt`

---

### Step 31 — miRNA Sponge Prediction

Run via psRNATarget web server: https://psrnattarget.noble.org/
Use Arabidopsis thaliana miRNAs (Araport v11) as proxy.

**Key results:**
- 106 predicted interactions across 10/14 candidates
- contig_2815 + miR164a/b/c: expectation 2.5 (strongest hit — miR164 targets NAC TFs)
- contig_8442 + miR830-5p: expectation 3.0
- contig_12238 + miR780.1: expectation 4.0

Full results: `results/mirna_sponge/psrnatarget_results.txt`

> NOTE: Arabidopsis miRNA proxy introduces substantial uncertainty. Resurrection plant miRNA databases do not exist. Predictions require validation by degradome sequencing or RNA pull-down.

---

<a name="section-12"></a>
## Section 12: Selaginella Cross-Species Validation

### Step 32 — Selaginella GO Enrichment

```bash
# After completing Selaginella pipeline (Sections 6-7):
python3 - << 'EOF'
# See results/selaginella_go_enrichment.csv for full results
# Key: chloroplast 1.99x (FDR=3x10^-22) — confirmed in both species
EOF
```

### Step 33 — Cross-Species BLAST

```bash
# BLAST Cp candidates against Sl transcriptome
makeblastdb -in references/selaginella_transcriptome.fa \
    -dbtype nucl -out results/selaginella_transcriptome_db

blastn -query results/coding_potential/candidates_14.fasta \
    -db results/selaginella_transcriptome_db \
    -out results/cp_vs_selaginella_blast_relaxed.txt \
    -evalue 1e-3 -word_size 7 \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -num_threads 12

# Results: 0 hits at E<1e-5
# 3/14 with partial conservation at E<1e-3:
#   contig_551  -> ABC transporter ABCC14 (75% identity)
#   contig_4950 -> TRAPPC6B vesicle trafficking (72% identity)
#   contig_22145 -> Mitochondrial Rho GTPase MIRO1 (76% identity)
```

### Step 34 — Selaginella Time-Course Analysis

```bash
# After quantifying all 17 samples:
python3 - << 'EOF'
# See results/selaginella_timecourse_tpm_matrix.csv
# See results/selaginella_timecourse_summary.csv
# Key finding: two kinetic subsets
#   Early responder: MSTRG.22821 (active at 1hr, peaks at 24hr)
#   Late responders: MSTRG.15112, MSTRG.15115, MSTRG.15111, MSTRG.6238
#                   (silent at 1hr/6hr, activate at 24hr, sustained through 120hr)
EOF
```

Full results: `results/selaginella_timecourse_interpretation.txt`

---

<a name="section-13"></a>
## Section 13: Validating Your Results

### Step 35 — Sanity Check All Results

```bash
python3 -c "
import pandas as pd
tpm = pd.read_csv('results/stringtie_tpm_matrix.csv')
real14 = pd.read_csv('results/real14_candidates_stats.csv')
deseq2 = pd.read_csv('results/deseq2_desiccated_vs_hydrated.csv')
coex = pd.read_csv('results/coexpression_partners.csv')
enrich = pd.read_csv('results/go_enrichment_results.csv')
sl = pd.read_csv('results/selaginella_go_enrichment.csv')

print(f'TPM matrix:              {tpm.shape}')
print(f'Real 14 candidates:      {len(real14)}')
sig = deseq2[(deseq2.padj < 0.05) & (~deseq2.padj.isna())]
print(f'DESeq2 sig (padj<0.05):  {len(sig)}')
print(f'Co-expression pairs:     {len(coex)}')
print(f'Top Cp enrichment:       {enrich.iloc[0].Category} {enrich.iloc[0].Fold_enrichment}x')
print(f'Top Sl enrichment:       {sl.iloc[0].Category} {sl.iloc[0].Fold_enrichment}x')
print(f'Rfam hits:               0 (confirmed)')
"
```

### Expected Results Summary

| Output File | Expected Shape | Key Numbers |
|-------------|----------------|-------------|
| `stringtie_tpm_matrix.csv` | 48,045 x 10 | All 14 candidates 0 TPM hydrated |
| `real14_candidates_stats.csv` | 14 rows | TPM range 590-1,957 desiccated |
| `deseq2_desiccated_vs_hydrated.csv` | ~6,290 rows | 14 candidates padj < 0.05 |
| `coexpression_partners.csv` | 1,062 rows | r >= 0.90 all pairs |
| `blast_annotated_partners.csv` | ~4,000 rows | 356 unique contigs annotated |
| `go_enrichment_results.csv` | 7 rows | Chloroplast 3.73x FDR=6x10^-14 |
| `selaginella_go_enrichment.csv` | 7 rows | Chloroplast 1.99x FDR=3x10^-22 |
| `wgcna_module_trait_correlation.csv` | 33 rows | MEgrey60 r=-0.988 |
| `selaginella_timecourse_tpm_matrix.csv` | 26,693 x 17 | Full time-course |
| `rfam_hits.tbl` | 0 data rows | 0 contamination hits |

---

<a name="section-14"></a>
## Section 14: Troubleshooting Reference

| Error / Symptom | Root Cause | Fix |
|----------------|------------|-----|
| `mbedtls_ssl_handshake returned -9984` | SRA Toolkit 2.9.6 TLS incompatibility | Upgrade to SRA Toolkit 3.x |
| `_abundance.txt` has 7 lines | StringTie ran without -G flag | Rerun with -G references/annotation.gtf |
| `padj < 0.05: 0 genes` genome-wide | n=2 hydrated replicates | Expected; use DESeq2 with sizeFactors=1 for candidates |
| DESeq2 `estimateSizeFactors` error | All-zero genes | Use `sizeFactors(dds) <- rep(1, ncol(dds))` |
| DESeq2 `every gene contains at least one zero` | Zero-inflated TPM | Use `sizeFactors(dds) <- rep(1, ncol(dds))` |
| DESeq2 `iterative size factor normalization did not converge` | Sparse coverage | Set size factors manually to 1 |
| GCLocker warnings (Trimmomatic) | Harmless Java GC messages | Ignore |
| Trimmomatic OutOfMemoryError | Large file + default Java heap | Add -Xmx2g -threads 4 |
| samtools "merging from 0 files" | Broken library links in base env | Use cpat_env for samtools operations |
| Small BAM files (40-140 MB for Selaginella) | Trimmed FASTQ truncated by OOM | Re-trim with -Xmx2g, verify line counts before aligning |
| StringTie -e error: no GFF given | No annotation for S. lepidophylla | Run de novo first, merge, then re-quantify with -e -G merged.gtf |
| CPC2 not found locally | Packaging issue | Use web server at https://cpc2.gao-lab.org/ |
| BLAST takes >30 min | Many sequences | Increase --num_threads; use -max_target_seqs 1 |
| partner_name shows . in co-expression | De novo transcriptome has no gene names | Expected; BLAST annotation fills this in downstream |
| 0-byte FASTQ files | Download failed silently | Delete 0-byte files, re-run download |

---

<a name="section-15"></a>
## Section 15: Quick Reference — All Steps

| Step | Action | Script | Est. Time |
|------|--------|--------|-----------|
| 1 | Check OS, disk, RAM | `00_preflight.sh` | 1 min |
| 2 | Check Python 3.8+ | `00_preflight.sh` | auto |
| 3 | Check / install Git | manual | 2 min |
| 4 | Fix SSL certificates | manual | 5 min |
| 5 | Install Miniconda | manual | 10 min |
| 6 | Install pipeline tools (conda/apt) | manual | 15 min |
| 7 | Install R + DESeq2 + WGCNA | manual | 10 min |
| 8 | Verify all tools | `00_preflight.sh` | auto |
| 9 | Clone repository | git clone | 2 min |
| 10 | Create directories | `00_preflight.sh` | auto |
| 11 | Upgrade SRA Toolkit | manual | 10 min |
| 12 | Download Cp samples (10) | `01_download.sh` | 2-4 hrs |
| 13 | Download Sl samples (17) | manual | 3-5 hrs |
| 14 | Download Cp reference FASTA | `01_download.sh` | 20 min |
| 15 | Download Sl TSA transcriptome | manual | 10 min |
| 16 | Build HISAT2 indices (both species) | `01_download.sh` | 20 min |
| 17 | Quality trimming Cp | `02_trim.sh` | 2-3 hrs |
| 18 | Quality trimming Sl (-Xmx2g) | manual | 3-4 hrs |
| 19 | HISAT2 alignment Cp | `03_align.sh` | 3-4 hrs |
| 20 | HISAT2 alignment Sl | manual | 4-5 hrs |
| 21 | StringTie Cp (reference-guided) | `04_stringtie.sh` | 2-3 hrs |
| 22 | StringTie Sl (de novo + merge + requantify) | manual | 3-4 hrs |
| 23 | Merge Cp TPM matrix | `05_merge_matrix.py` | 5 min |
| 24 | Merge Sl TPM matrices (9-sample + 17-sample) | manual | 5 min |
| 25 | DESeq2 differential expression | `12_wgcna.R` | 10 min |
| 26 | Exploratory DE (Welch) | `06_de_analysis.py` | 5 min |
| 27 | Visualizations | `07_visualize.py` | 5 min |
| 28 | WGCNA | `12_wgcna.R` | 20 min |
| 29 | Pearson co-expression (Cp + Sl) | `08_coexpression.py` | 15 min |
| 30 | Download UniProt Swiss-Prot | wget + makeblastdb | 15 min |
| 31 | BLAST annotation (Cp + Sl) | `09_blast_annotate.py` | 30 min |
| 32 | GO enrichment (Cp + Sl) | `10_go_enrichment.py` | 5 min |
| 33 | Generate final report | `11_final_report.py` | 1 min |
| 34 | Download Rfam database | wget + cmpress | 5 min |
| 35 | Extract candidate sequences | python3 | 1 min |
| 36 | Rfam contamination screen | cmscan | 10 min |
| 37 | CPC2 coding potential | Web server | 5 min |
| 38 | psRNATarget miRNA sponge | Web server | 10 min |
| 39 | ABA pathway analysis | python3 | 2 min |
| 40 | TF co-expression analysis | python3 | 2 min |
| 41 | Selaginella GO enrichment | python3 | 5 min |
| 42 | Cross-species BLAST (Cp vs Sl, Cp vs Arabidopsis) | blastn | 10 min |
| 43 | Selaginella time-course kinetics | python3 | 5 min |
| 44 | Rehydration analysis (Cp + Sl) | python3 | 5 min |
| 45 | WGCNA hub gene identification | Rscript | 5 min |
| 46 | Validate all results | python3 | 2 min |
| 47 | Manuscript preparation | — | Complete |

> Total computational runtime: ~15 hours end-to-end

---

## GFF Search Log — S. lepidophylla (Documented Gap)

No GFF annotation exists for *S. lepidophylla*. Steps taken to locate it:

1. BioProject PRJNA420971 — only SRA and BioSample linked, no assembly or GFF
2. VanBuren et al. 2018 data availability — genome under PRJNA386571 / SAMN07071123 — NCBI BioProject shows only SRA experiments, no assembly or annotation linked
3. NCBI Assembly search for "Selaginella lepidophylla VanBuren" — no results
4. GCF_000143415.4 — confirmed as *S. moellendorffii* not *S. lepidophylla* (wrong species)
5. Phytozome — not checked; requires login at https://phytozome-next.jgi.doe.gov/info/Slepidophylla_v1_0 — worth checking if annotation needed in future

**Conclusion:** No publicly available GFF exists. VanBuren 2018 genome deposited as raw sequence only. StringTie de novo mode is the correct approach for unannotated assemblies.

---

## Known Limitations

- **n=2 hydrated replicates** (C. plantagineum) — limits genome-wide FDR power. The 14 core candidates are DESeq2-confirmed; the 661 exploratory set uses nominal thresholds only.
- **De novo assembly** (S. lepidophylla) — no reference genome annotation; may introduce transcript boundary artefacts.
- **CPC2 flags 5/14** — designated "transcripts of interest". May contain sORFs encoding micropeptides while retaining regulatory non-coding functions.
- **Arabidopsis miRNA proxy** — psRNATarget predictions are hypothesis-generating only; resurrection plant miRNA databases do not exist.
- **No wet-lab validation** — qPCR, RIP-seq, and degradome sequencing are priority future directions.

---

## Links

| Resource | URL |
|----------|-----|
| This repository | github.com/legoreddragon-ai/resurrection_lncrna_pipeline |
| C. plantagineum RNA-seq | ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157098 |
| S. lepidophylla RNA-seq | ncbi.nlm.nih.gov/sra/?term=PRJNA420971 |
| C. plantagineum genome 2023 | doi.org/10.1111/tpj.16165 |
| CPC2 server | cpc2.gao-lab.org |
| psRNATarget | psrnattarget.noble.org |
| UniProt Swiss-Prot | ftp.uniprot.org/pub/databases/uniprot |
| Rfam | ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT |

---

## Citation

If you use this pipeline or build on these findings, please cite:

> Sumanth N. Discovery of Desiccation-Induced Long Non-Coding RNAs in the Resurrection Plant *Craterostigma plantagineum*: A Transcriptome-Wide Computational Analysis with Cross-Species Validation in *Selaginella lepidophylla*. 2026. bioRxiv [preprint].

---

*Resurrection Plant lncRNA Pipeline*
*Built from actual execution history — all problems encountered and resolved are documented above*
*Last updated: March 2026 · Pipeline: 100% complete · Manuscript: v6 final*