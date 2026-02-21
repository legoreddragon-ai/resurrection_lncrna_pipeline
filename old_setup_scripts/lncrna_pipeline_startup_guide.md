# lncRNA Discovery Pipeline — Complete Setup & Execution Guide
### Craterostigma plantagineum · RNA-Seq · Differential Expression · lncRNA Classification
**From Git Clone → Environment → Data Download → Alignment → DE Analysis → Results**

> **27 Steps · All Known Problems Documented · Scripts Included**
> Total runtime: 8–12 hours end-to-end

---

## Table of Contents

1. [System Requirements & Pre-Flight Check](#section-1)
2. [Conda Installation & Environment Setup](#section-2)
3. [Getting the Code](#section-3)
4. [Downloading RNA-Seq Data from NCBI](#section-4)
5. [Reference Genome & Index](#section-5)
6. [Running the Main Pipeline](#section-6)
7. [Expression Matrix, DE Analysis & Visualization](#section-7)
8. [Validating Your Results](#section-8)
9. [Troubleshooting Reference](#section-9)
10. [Quick Reference — All 27 Steps](#section-10)

---

<a name="section-1"></a>
## Section 1: System Requirements & Pre-Flight Check

Before installing anything, verify your environment meets the minimum requirements.

### Hardware Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| RAM | 8 GB | 16 GB | StringTie and HISAT2 are memory-intensive |
| Disk Space | 80 GB | 150 GB | Raw FASTQs (~4 GB each × 10 samples) + intermediates |
| CPU Cores | 4 | 8+ | All tools support multi-threading |
| OS | Ubuntu 20.04+ | Ubuntu 22.04+ | Also works on WSL2 on Windows |

---

### Step 1 — Verify OS and Architecture

```bash
uname -a
# Expected: Linux <hostname> ... x86_64 GNU/Linux

lsb_release -a
# Expected: Ubuntu 20.04 or 22.04

df -h /
# Expected: >80 GB available

free -h
# Expected: >8 GB total RAM
```

> ⚠️ **PROBLEM ENCOUNTERED:** Pipeline had hardcoded paths from the original development laptop (`/home/neil/resurrection_lncrna_pipeline` and `/mnt/p/legor/Documents/...`). These failed immediately on any other machine.
>
> ✅ **FIX:** Always use dynamic path resolution. Fixed in current scripts via `Path(__file__).parent.resolve()` in Python and `$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)` in Bash.

---

### Step 2 — Check Python

```bash
python3 --version
# Expected: Python 3.8 or higher

pip3 --version
# If missing:
sudo apt-get install python3-pip
```

---

### Step 3 — Check Git

```bash
git --version
# Expected: git version 2.x.x

# If missing:
sudo apt-get update && sudo apt-get install git -y
```

---

### Step 4 — Check SSL/Certificate Status

> This is the most common failure point. Fix proactively before attempting any NCBI downloads.

```bash
# Test NCBI connectivity
curl -v https://www.ncbi.nlm.nih.gov 2>&1 | grep -E '(SSL|TLS|certificate|Connected)'
# Expected: "SSL connection using TLSv1.3"
# If you see "certificate verify failed" — proceed to fix below

# Fix: update system CA certificates
sudo apt-get update
sudo apt-get install --reinstall ca-certificates -y
sudo update-ca-certificates

# Verify fix
curl -I https://www.ncbi.nlm.nih.gov
# Expected: HTTP/2 200
```

> ⚠️ **PROBLEM ENCOUNTERED:** SRA Toolkit version 2.9.6 uses `mbedtls` which rejects NCBI's current TLS certificate signing algorithm.
> ```
> mbedtls_ssl_handshake returned -9984
> X509 - Certificate verification failed
> ```
> System `curl` worked fine but `prefetch` and `fasterq-dump` both failed.
>
> ✅ **ROOT CAUSE:** Old SRA toolkit (2.9.6) bundled TLS library is incompatible with current NCBI certificates. Fix: upgrade the toolkit OR use direct `wget`/`curl` downloads. Detailed in Step 11.

---

<a name="section-2"></a>
## Section 2: Conda Installation & Environment Setup

### Step 5 — Install Miniconda

```bash
# Check if conda is already available
conda --version

# If not found, install Miniconda:
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3

# Initialize conda
~/miniconda3/bin/conda init bash
source ~/.bashrc

# Verify
conda --version
# Expected: conda 23.x.x or higher
```

---

### Step 6 — Install Bioinformatics Tools

```bash
# Install Python packages first
pip install pandas numpy scipy biopython matplotlib seaborn --break-system-packages

# Try conda for bioinformatics tools first
conda install -c bioconda trimmomatic hisat2 samtools stringtie -y

# If conda fails with LibMambaUnsatisfiableError:
sudo apt-get update
sudo apt-get install trimmomatic hisat2 samtools -y

# StringTie via conda only (not in apt):
conda install -c bioconda stringtie -y
```

> ⚠️ **PROBLEM ENCOUNTERED:** `conda install hisat2 samtools` failed with `LibMambaUnsatisfiableError` — dependency conflict involving `ossuuid` and other packages.
>
> ✅ **FIX:** Install `hisat2` and `samtools` via `apt-get`. Install `stringtie` separately via conda. Both coexist without conflict.

---

### Step 7 — Verify All Tools

```bash
trimmomatic -version   # Expected: 0.39
hisat2 --version       # Expected: HISAT2 version 2.2.x
samtools --version     # Expected: samtools 1.x
stringtie --version    # Expected: StringTie v2.x
python3 --version      # Expected: 3.8+
python3 -c 'import pandas, numpy, scipy; print("Python deps OK")'
```

All commands must return version numbers with no `command not found` errors before proceeding.

---

<a name="section-3"></a>
## Section 3: Getting the Code

### Step 8 — Clone the Repository

```bash
cd /mnt/neil   # Or your preferred location
git clone https://github.com/YOUR_USERNAME/resurrection_lncrna_pipeline.git
cd resurrection_lncrna_pipeline

# Verify structure
ls -la
# Expected directories: data/ references/ results/ scripts/ logs/
```

---

### Step 9 — Create All Required Directories

Certain directories are in `.gitignore` and will not be cloned. Create them manually:

```bash
mkdir -p data/fastq_raw
mkdir -p data/fastq_trimmed
mkdir -p results/hisat2_aligned
mkdir -p results/stringtie_output
mkdir -p results/figures
mkdir -p logs

# Verify
python3 run.py --check-structure
# Expected: All directories show as PRESENT
```

> ⚠️ **PROBLEM ENCOUNTERED:** `python3 run.py --check-structure` showed all directories as MISSING even though they existed. Root cause: `data_loader.py` still had hardcoded paths from the original machine.
>
> ✅ **FIX:** `config.py` was updated to use `Path(__file__).parent.resolve()` so all paths are computed relative to the script location. If you still see this error, verify `config.py` does NOT contain absolute paths.

---

### Step 10 — Test Python Configuration

```bash
python3 scripts/config.py

# Expected output:
# ================================================================================
# PROJECT CONFIGURATION
# ================================================================================
# Directories:
#   Project: /mnt/neil/resurrection_lncrna_pipeline
#   Data:    /mnt/neil/resurrection_lncrna_pipeline/data
#   Results: /mnt/neil/resurrection_lncrna_pipeline/results
# ...
# Directories created/verified.
```

---

<a name="section-4"></a>
## Section 4: Downloading RNA-Seq Data from NCBI

The data comes from GEO accession **GSE157098** / SRA project **PRJNA660052**. There are 10 paired-end samples totalling approximately 40 GB.

### Step 11 — Upgrade or Bypass the SRA Toolkit

```bash
# Try upgrading via conda first
conda install -c bioconda sra-tools=3.1.1 -y
fasterq-dump --version
# Expected: fasterq-dump.3.x.x

# If conda upgrade fails (missing ossuuid dependency),
# download the binary directly from NCBI:
cd ~
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=~/sratoolkit.3.2.1-ubuntu64/bin:$PATH
echo 'export PATH=~/sratoolkit.3.2.1-ubuntu64/bin:$PATH' >> ~/.bashrc

# Verify
fasterq-dump --version
# Expected: fasterq-dump.3.x.x

# Configure to accept all certificates
vdb-config --set /tls/allow-all-certs=true

# Test with a single accession
fasterq-dump SRR12542206 --outdir /tmp/ --split-files --progress
# Expected: join 100%, concat 100%, spots read: ~1,198,909
rm /tmp/SRR12542206_*.fastq
```

> ⚠️ **KNOWN FAILURE — SRA File Conversion:** Even with the toolkit upgraded, `prefetch` sometimes downloads `.sra` files that `fasterq-dump` cannot convert. Error: `"column undefined while opening cursor"`. If you see this, use the wget fallback in the download script below.

---

### Step 12 — Download All 10 Samples

Save this as `download_samples.sh` and run it:

```bash
#!/bin/bash
# download_samples.sh
# Downloads all 10 Craterostigma plantagineum RNA-seq samples
# Safe to re-run — skips already completed samples

BASE=/mnt/neil/resurrection_lncrna_pipeline
OUT=${BASE}/data/fastq_raw
mkdir -p ${OUT}

declare -A SAMPLES=(
    ["SRR12542161"]="cp_hyd_r1"
    ["SRR12542166"]="cp_hyd_r3"
    ["SRR12542171"]="cp_rehyd_r1"
    ["SRR12542176"]="cp_rehyd_r2"
    ["SRR12542181"]="cp_rehyd_r3"
    ["SRR12542186"]="cp_wc2_r2"
    ["SRR12542191"]="cp_wc2_r3"
    ["SRR12542196"]="cp_wc60_r1"
    ["SRR12542206"]="cp_wc60_r2"
    ["SRR12542211"]="cp_wc60_r3"
)

for ACC in "${!SAMPLES[@]}"; do
    NAME="${SAMPLES[$ACC]}"
    R1="${OUT}/${NAME}_1.fastq.gz"
    R2="${OUT}/${NAME}_2.fastq.gz"

    # Skip if already downloaded and non-empty
    if [ -s "$R1" ] && [ -s "$R2" ]; then
        echo "SKIP: $NAME already downloaded"
        continue
    fi

    echo "Downloading $NAME ($ACC)..."
    fasterq-dump "$ACC" \
        --outdir "${OUT}/" \
        --outfile "${NAME}.fastq" \
        --split-files --threads 8 --progress

    if [ $? -ne 0 ]; then
        echo "ERROR: fasterq-dump failed for $ACC — trying wget fallback"
        SUBDIR=$(echo $ACC | rev | cut -c1-3 | rev)
        BASE_URL="https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/SRR012/${ACC:3:5}/${ACC}/${ACC}.lite.1"
        wget --no-check-certificate -c "${BASE_URL}" -O "${OUT}/${ACC}.sra"
        fasterq-dump "${OUT}/${ACC}.sra" \
            --outdir "${OUT}/" \
            --outfile "${NAME}.fastq" \
            --split-files
        continue
    fi

    echo "Compressing $NAME..."
    gzip ${OUT}/${NAME}_1.fastq
    gzip ${OUT}/${NAME}_2.fastq
    echo "Done: $NAME"
done

echo "All downloads complete."
ls -lh ${OUT}/
```

```bash
bash download_samples.sh
```

---

### Step 13 — Verify Downloads

```bash
# List all FASTQ files
ls -lh data/fastq_raw/ | grep -E '\.(fastq|fastq\.gz)$'
# Expected: 20 files (2 per sample × 10 samples)
# Each file should be >10 MB

# Flag any 0-byte files
find data/fastq_raw/ -name '*.fastq.gz' -size 0 -print
# Expected: no output

# Count total files
ls data/fastq_raw/*.fastq.gz | wc -l
# Expected: 20
```

> ⚠️ **PROBLEM ENCOUNTERED:** Batch scripts created 0-byte files when downloads failed silently. `fasterq-dump` exited with code 0 but created empty output files.
>
> ✅ **FIX:** Always use `[ -s "$file" ]` (size > 0) not `[ -f "$file" ]` (file exists) to detect complete downloads. The script above handles this correctly. If you find 0-byte files: delete them and re-run — the skip logic will only re-download what is missing.

---

<a name="section-5"></a>
## Section 5: Reference Genome & Index

### Step 14 — Download Reference Files

```bash
cd references/

# Download transcriptome FASTA from GEO supplementary
wget --no-check-certificate \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_transcriptome_assembly_V2.fa.gz'

# Download annotation
wget --no-check-certificate \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz'

# Decompress
gunzip *.gz

# Verify
head -2 GSE157098_Cp_transcriptome_assembly_V2.fa
# Expected: >Cp_V2_contig_1

grep -c '>' GSE157098_Cp_transcriptome_assembly_V2.fa
# Expected: ~48045
```

---

### Step 15 — Build HISAT2 Index

This takes 5–15 minutes.

```bash
hisat2-build \
    references/GSE157098_Cp_transcriptome_assembly_V2.fa \
    references/Cp_transcriptome_index \
    -p 8

# Verify — should see 8 index files
ls references/Cp_transcriptome_index*.ht2 | wc -l
# Expected: 8
```

---

### Step 16 — Prepare GTF Annotation

```bash
# Check if annotation.gtf exists
ls references/annotation.gtf

# If not, convert from gff3
conda install -c bioconda gffread -y
gffread references/annotation.gff3 -T -o references/annotation.gtf

# Verify
head -5 references/annotation.gtf
wc -l references/annotation.gtf
# Expected: >50,000 lines
```

---

<a name="section-6"></a>
## Section 6: Running the Main Pipeline

The main pipeline runs Trim → Align → Assemble across all 10 samples. Run with `nohup` since this takes 8+ hours.

```bash
nohup bash run_remaining_pipeline.sh > logs/pipeline_run.log 2>&1 &
echo "PID: $!"

# Monitor progress
tail -f logs/pipeline_run.log
```

---

### Step 17 — Quality Trimming (Trimmomatic)

**Estimated time: 2–3 hours for 10 samples.**

What to expect per sample:

```
[2026-02-21 02:18:40] Trimming cp_hyd_r1...
[3.2s][warning][gc,alloc] pool-7-thread-1: Retried waiting for GCLocker...
# NOTE: GCLocker warnings are harmless Java GC messages — NOT errors
# Trimmomatic will still complete successfully
[2026-02-21 02:25:10] Done trimming: cp_hyd_r1
```

After trimming, verify:

```bash
ls -lh data/fastq_trimmed/ | grep '_paired'
# Expected: 20 paired files, each >10 MB
```

> ⚠️ **PROBLEM ENCOUNTERED:** `set -euo pipefail` in the pipeline script caused it to exit immediately when Trimmomatic returned a non-zero exit code due to GC warnings. The script appeared to run but stopped after the first sample.
>
> ✅ **FIX:** Remove the `-e` flag:
> ```bash
> sed -i 's/set -euo pipefail/set -uo pipefail/' run_remaining_pipeline.sh
> ```

---

### Step 18 — HISAT2 Alignment

**Estimated time: 3–4 hours for all 10 samples.**

The pipeline script handles this automatically. To run a single sample manually:

```bash
hisat2 \
    -x references/Cp_transcriptome_index \
    -1 data/fastq_trimmed/SAMPLE_1_paired.fastq.gz \
    -2 data/fastq_trimmed/SAMPLE_2_paired.fastq.gz \
    -p 8 --dta \
    2> logs/SAMPLE_hisat2.log \
| samtools sort -@ 8 -o results/hisat2_aligned/SAMPLE.sorted.bam

samtools index results/hisat2_aligned/SAMPLE.sorted.bam

# Check alignment rates
grep 'overall alignment rate' logs/*_hisat2.log
# Expected: 60-95% per sample. <40% indicates a problem.
```

> ⚠️ **PROBLEM ENCOUNTERED:** `cp_hyd_r1` BAM file did not exist because it was processed before the pipeline script existed. StringTie was run without `-G` (no annotation), producing 6 STRG.x transcripts instead of 48,045.
>
> ✅ **FIX:** Always verify abundance file line count after StringTie:
> ```bash
> wc -l results/stringtie_output/SAMPLE_abundance.txt
> # Must be ~48,046. If you see 7 or fewer: -G flag was missing, rerun.
> ```

---

### Step 19 — StringTie Assembly

**Estimated time: 2–3 hours for all 10 samples.**

```bash
# Per sample:
stringtie results/hisat2_aligned/SAMPLE.sorted.bam \
    -G references/annotation.gtf \
    -o results/stringtie_output/SAMPLE.gtf \
    -A results/stringtie_output/SAMPLE_abundance.txt \
    -p 8 -e \
    2> logs/SAMPLE_stringtie.log

# Critical verification
wc -l results/stringtie_output/SAMPLE_abundance.txt
# Expected: 48046 (48045 genes + 1 header)
# If you see 7 lines: the -G flag was missing. Delete and rerun.
```

---

<a name="section-7"></a>
## Section 7: Expression Matrix, DE Analysis & Visualization

### Step 20 — Verify All StringTie Outputs

```bash
ls results/stringtie_output/*_abundance.txt
# Expected: 10 files

# All must have ~48046 lines
wc -l results/stringtie_output/*_abundance.txt | sort -n | head -15
# Expected: all values near 48046
```

---

### Step 21 — Merge TPM Matrix

```bash
python3 scripts/task4_merge.py

# Expected output:
# Loaded: cp_hyd_r1 (48045 genes)
# Loaded: cp_hyd_r3 (48045 genes)
# ... (all 10 samples)
# Saved: results/stringtie_tpm_matrix.csv
# Shape: (48045, 12)
# 12 columns = Gene ID + Gene Name + 10 sample columns
```

---

### Step 22 — Differential Expression Analysis

```bash
python3 scripts/task5_de.py

# Expected output:
# Loading TPM matrix...
# Shape: (48045, 12)
# lncRNA candidates: 38978
# Expressed lncRNAs: 17574
# Saved full DE: 17574 lncRNAs
# Candidates (pval<0.05, |log2FC|>2, max TPM>=5): ~777
```

> ⚠️ **IMPORTANT — Statistical Limitation:**
> With only **n=2 hydrated replicates**, BH FDR correction across 17,574 tests yields **0 significant genes**. This is mathematically correct, not a pipeline error. The t-test with n=2 has only 1 degree of freedom — insufficient statistical power for genome-wide correction.
>
> The analysis therefore selects top 14 candidates using: **raw p < 0.05 + |log2FC| > 2 + expression score**. This limitation must be documented in any write-up or presentation.

---

### Step 23 — Generate Visualizations

```bash
python3 scripts/task6_visualize.py

# Expected output files:
# results/volcano_plot.png
# results/heatmap_significant_lncrnas.png
# results/boxplots_top_lncrnas.png

ls -lh results/*.png
# Each file should be >100 KB
```

---

<a name="section-8"></a>
## Section 8: Validating Your Results

### Step 24 — Sanity Check Results

```bash
# Check TPM matrix shape
python3 -c "
import pandas as pd
df = pd.read_csv('results/stringtie_tpm_matrix.csv')
print(f'Shape: {df.shape}')           # Expected: (48045, 12)
print(f'Columns: {list(df.columns)}') # Expected: Gene ID, Gene Name, cp_hyd_r1...
"

# Check DE results
python3 -c "
import pandas as pd
df = pd.read_csv('results/lncrna_differential_expression.csv')
print(f'DE results: {len(df)} lncRNAs')  # Expected: 17574
print(df.head(3))
"

# Check top 14
python3 -c "
import pandas as pd
df = pd.read_csv('results/significant_stress_responsive_lncrnas.csv')
print(f'Top candidates: {len(df)}')  # Expected: 14
print(df[['gene_id','mean_hydrated','mean_dehydrated','log2FC','pvalue']])
"
```

---

### Step 25 — Verify Top Candidate

```bash
python3 -c "
import pandas as pd
df = pd.read_csv('results/significant_stress_responsive_lncrnas.csv')
top = df.iloc[0]
print(f'Top candidate: {top.gene_id}')
print(f'Hydrated:    {top.mean_hydrated:.1f} TPM')
print(f'Dehydrated:  {top.mean_dehydrated:.1f} TPM')
print(f'Fold change: {2**top.log2FC:.1f}x')
print(f'p-value:     {top.pvalue:.4f}')
"
# Expected: top candidate shows >5-fold increase in dehydrated vs hydrated
```

---

### Expected Results Summary

| Output File | Expected Content | Key Numbers |
|-------------|-----------------|-------------|
| `stringtie_tpm_matrix.csv` | Full expression matrix | 48,045 genes × 10 samples |
| `coding_potential.txt` | lncRNA classification | 38,978 lncRNAs / 9,067 coding (81%/19%) |
| `lncrna_differential_expression.csv` | All lncRNA DE stats | 17,574 expressed lncRNAs tested |
| `significant_stress_responsive_lncrnas.csv` | Top candidates | 14 stress-induced lncRNAs |
| `volcano_plot.png` | DE visualization | Scatter of log2FC vs -log10(p) |
| `heatmap_significant_lncrnas.png` | Expression heatmap | 14 genes × 10 samples, z-scored |
| `boxplots_top_lncrnas.png` | Per-gene boxplots | Expression by condition for top 14 |

---

<a name="section-9"></a>
## Section 9: Troubleshooting Reference

Every problem encountered in this project, with root cause and fix:

| Error / Symptom | Root Cause | Fix |
|----------------|------------|-----|
| `mbedtls_ssl_handshake returned -9984` | SRA Toolkit 2.9.6 uses old TLS library incompatible with NCBI certs | Upgrade to SRA Toolkit 3.x or use `wget --no-check-certificate` |
| `invalid accession 'SRR12542206'` | `fasterq-dump` cannot reach NCBI due to TLS failure | Upgrade SRA toolkit first (Step 11) |
| All directories MISSING in `--check-structure` | `data_loader.py` has hardcoded paths from original machine | Use dynamic paths: `Path(__file__).parent.resolve()` |
| `LibMambaUnsatisfiableError: ossuuid` | `sra-tools 3.x` requires `ossuuid` package not in default channels | Download SRA toolkit binary directly from NCBI |
| `LibMambaUnsatisfiableError` (hisat2/samtools) | Conda dependency conflicts in base environment | Install via `apt-get` instead |
| 0-byte FASTQ files after download script | Script created file before downloading; failed download left 0-byte file | Use `[ -s "$file" ]` (size > 0) not `[ -f "$file" ]` in skip checks |
| `cp_hyd_r1_abundance.txt` has 7 lines | StringTie run without `-G` annotation flag | Rerun with `-G references/annotation.gtf -e` flags |
| `padj < 0.05: 0 genes` | Only 2 hydrated replicates — insufficient power for FDR across 17,574 tests | Expected behavior. Use raw p-value + fold-change + expression score |
| GCLocker warnings from Trimmomatic | Harmless Java GC messages — not errors | Ignore. Remove `-e` from `set -euo pipefail` |
| Pipeline exits after first sample | `set -euo pipefail` exits on GCLocker non-zero return code | Change to `set -uo pipefail` |
| EBI FTP paths 404 / 0 bytes | EBI mirror path structure varies by accession | Use NCBI direct download or `fasterq-dump` instead |

---

<a name="section-10"></a>
## Section 10: Quick Reference — All 27 Steps

| Step | Action | Est. Time | Command / Script |
|------|--------|-----------|-----------------|
| 1 | Check OS, disk, RAM | 2 min | `uname -a`, `df -h`, `free -h` |
| 2 | Check Python 3.7+ | 1 min | `python3 --version` |
| 3 | Check / install Git | 2 min | `git --version` |
| 4 | Fix SSL certificates | 5 min | `sudo apt-get install ca-certificates` |
| 5 | Install Miniconda | 10 min | `Miniconda3-latest-Linux-x86_64.sh` |
| 6 | Install pipeline tools | 15 min | `conda` + `apt-get` |
| 7 | Verify all tools | 2 min | `trimmomatic`, `hisat2`, `samtools`, `stringtie --version` |
| 8 | Clone repository | 2 min | `git clone <repo_url>` |
| 9 | Create gitignored directories | 1 min | `mkdir -p data/fastq_raw ...` |
| 10 | Test Python config | 1 min | `python3 scripts/config.py` |
| 11 | Upgrade SRA Toolkit | 10 min | `conda install sra-tools=3.1.1` or binary download |
| 12 | Download 10 samples | 2–4 hrs | `bash download_samples.sh` |
| 13 | Verify downloads (no 0-byte) | 2 min | `find data/fastq_raw/ -size 0` |
| 14 | Download reference FASTA | 20 min | `wget` from GEO FTP |
| 15 | Build HISAT2 index | 15 min | `hisat2-build ...` |
| 16 | Prepare GTF annotation | 5 min | `gffread annotation.gff3 -T -o annotation.gtf` |
| 17 | Trimmomatic — all 10 samples | 2–3 hrs | `bash run_remaining_pipeline.sh` |
| 18 | HISAT2 alignment — all 10 | 3–4 hrs | Auto via pipeline script |
| 19 | StringTie assembly — all 10 | 2–3 hrs | Auto via pipeline script |
| 20 | Verify StringTie outputs | 2 min | `wc -l results/stringtie_output/*_abundance.txt` |
| 21 | Merge TPM expression matrix | 5 min | `python3 scripts/task4_merge.py` |
| 22 | Differential expression analysis | 5 min | `python3 scripts/task5_de.py` |
| 23 | Generate visualizations | 5 min | `python3 scripts/task6_visualize.py` |
| 24 | Sanity check results | 5 min | `python3 -c 'import pandas; ...'` |
| 25 | Verify top candidate | 2 min | `python3 -c '...'` |
| 26 | Co-expression network *(next)* | 30 min | `scripts/task7_coexpression.py` |
| 27 | GO enrichment analysis *(next)* | 30 min | `scripts/task8_go_enrichment.py` |

---

> **Total Pipeline Runtime: 8–12 hours end-to-end**
>
> Steps 1–10: ~1 hour setup | Steps 11–16: ~3 hours data prep | Steps 17–25: ~8 hours compute

---

*Resurrection Plant lncRNA Pipeline — Setup Guide*
*Built from actual execution history including all problems encountered and resolved*
