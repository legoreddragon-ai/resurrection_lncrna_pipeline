# lncRNA Discovery Pipeline — Complete Setup & Execution Guide
### Craterostigma plantagineum · RNA-Seq · Differential Expression · lncRNA Classification
**From Git Clone → Environment → Data Download → Alignment → DE Analysis → Results**

> **27 Steps · All Known Problems Documented · Scripts Included**
> Total runtime: 8–12 hours end-to-end

---

## Setup Folder — Script Index

All scripts live in the `setup/` directory. The master script runs the full pipeline end-to-end.

```
setup/
├── run_pipeline.sh        ← MASTER: runs all 7 stages in sequence
├── 00_preflight.sh        ← Stage 0: environment & prerequisites check
├── 01_download.sh         ← Stage 1: download FASTQs + reference + build index
├── 02_trim.sh             ← Stage 2: Trimmomatic quality trimming
├── 03_align.sh            ← Stage 3: HISAT2 alignment + SAMtools sort/index
├── 04_stringtie.sh        ← Stage 4: StringTie transcript quantification
├── 05_merge_matrix.py     ← Stage 5: merge abundance files into TPM matrix
├── 06_de_analysis.py      ← Stage 6: differential expression + candidate selection
└── 07_visualize.py        ← Stage 7: volcano plot, heatmap, boxplots
```

### Running the Full Pipeline

```bash
# Full run — all 7 stages
nohup bash setup/run_pipeline.sh > logs/master.log 2>&1 &
echo "PID: $!"

# Monitor
tail -f logs/master.log

# Skip stages already completed
bash setup/run_pipeline.sh --skip-download --skip-trim --skip-align
```

### Running Individual Stages

```bash
BASE=/mnt/neil/resurrection_lncrna_pipeline

bash setup/00_preflight.sh $BASE      # Check environment
bash setup/01_download.sh $BASE       # Download data + references
bash setup/02_trim.sh $BASE           # Trim reads
bash setup/03_align.sh $BASE          # Align reads
bash setup/04_stringtie.sh $BASE      # Assemble transcripts
python3 setup/05_merge_matrix.py $BASE  # Merge TPM matrix
python3 setup/06_de_analysis.py $BASE   # DE analysis
python3 setup/07_visualize.py $BASE     # Visualizations
```

All scripts are **idempotent** — safe to re-run. Each stage skips samples already completed.

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

> **Script:** `bash setup/00_preflight.sh`
> Automatically checks all requirements below and reports PASS / WARN / FAIL for each.

### Hardware Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| RAM | 8 GB | 16 GB | StringTie and HISAT2 are memory-intensive |
| Disk Space | 80 GB | 150 GB | Raw FASTQs (~4 GB each × 10 samples) + intermediates |
| CPU Cores | 4 | 8+ | All tools support multi-threading (`THREADS=8` in all scripts) |
| OS | Ubuntu 20.04+ | Ubuntu 22.04+ | Also works on WSL2 on Windows |

---

### Step 1 — Verify OS and Architecture

```bash
bash setup/00_preflight.sh
# OR manually:
uname -a            # Expected: Linux <hostname> ... x86_64 GNU/Linux
lsb_release -a      # Expected: Ubuntu 20.04 or 22.04
df -h /             # Expected: >80 GB available
free -h             # Expected: >8 GB total RAM
```

**Expected preflight output:**
```
============================================================
  PREFLIGHT CHECK
============================================================
[ System ]
  [PASS] OS: Linux x86_64
  [PASS] RAM: 16 GB
  [PASS] Disk available: 120 GB
[ Required Tools ]
  [PASS] python3: Python 3.10.x
  [PASS] hisat2: HISAT2 version 2.2.1
  ...
============================================================
  PREFLIGHT SUMMARY: 18 passed | 0 warnings | 0 failed
============================================================
```

> ⚠️ **PROBLEM ENCOUNTERED:** Pipeline had hardcoded paths from the original development laptop (`/home/neil/resurrection_lncrna_pipeline` and `/mnt/p/legor/Documents/...`). These failed immediately on any other machine.
>
> ✅ **FIX:** All scripts in `setup/` use dynamic path resolution:
> - Python: `Path(__file__).parent.parent.resolve()`
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

> This is the most common failure point. `00_preflight.sh` checks this automatically.

```bash
# Test NCBI connectivity
curl -v https://www.ncbi.nlm.nih.gov 2>&1 | grep -E '(SSL|TLS|certificate|Connected)'
# Expected: "SSL connection using TLSv1.3"

# If you see "certificate verify failed":
sudo apt-get update
sudo apt-get install --reinstall ca-certificates -y
sudo update-ca-certificates

# Verify
curl -I https://www.ncbi.nlm.nih.gov   # Expected: HTTP/2 200
```

> ⚠️ **PROBLEM ENCOUNTERED:** SRA Toolkit 2.9.6 uses `mbedtls` which rejects NCBI's TLS certificate:
> ```
> mbedtls_ssl_handshake returned -9984
> X509 - Certificate verification failed
> ```
> System `curl` worked fine but `prefetch` and `fasterq-dump` both failed.
>
> ✅ **FIX:** Upgrade to SRA Toolkit 3.x (Step 11) or use `wget --no-check-certificate`. The `01_download.sh` script handles this automatically with a fallback.

---

<a name="section-2"></a>
## Section 2: Conda Installation & Environment Setup

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
pip install pandas numpy scipy biopython matplotlib seaborn --break-system-packages

# Bioinformatics tools — try conda first
conda install -c bioconda trimmomatic hisat2 samtools stringtie -y

# If conda fails with LibMambaUnsatisfiableError:
sudo apt-get update && sudo apt-get install trimmomatic hisat2 samtools -y
conda install -c bioconda stringtie -y   # StringTie via conda only
```

> ⚠️ **PROBLEM ENCOUNTERED:** `conda install hisat2 samtools` failed with `LibMambaUnsatisfiableError` — dependency conflict involving `ossuuid`.
>
> ✅ **FIX:** Install `hisat2` and `samtools` via `apt-get`. Install `stringtie` separately via conda. Both coexist fine.

---

### Step 7 — Verify All Tools

```bash
# 00_preflight.sh checks all of these automatically. Or run manually:
trimmomatic -version    # Expected: 0.39
hisat2 --version        # Expected: HISAT2 version 2.2.x
samtools --version      # Expected: samtools 1.x
stringtie --version     # Expected: StringTie v2.x
python3 -c 'import pandas, numpy, scipy; print("Python deps OK")'
```

---

<a name="section-3"></a>
## Section 3: Getting the Code

### Step 8 — Clone the Repository

```bash
cd /mnt/neil
git clone https://github.com/YOUR_USERNAME/resurrection_lncrna_pipeline.git
cd resurrection_lncrna_pipeline
ls -la
# Expected directories: data/ references/ results/ scripts/ setup/ logs/
```

---

### Step 9 — Create All Required Directories

Directories in `.gitignore` are not cloned. Create them with the preflight script, or manually:

```bash
# Option A — automatic (recommended):
bash setup/00_preflight.sh
# Preflight auto-creates any missing directories

# Option B — manual:
mkdir -p data/fastq_raw data/fastq_trimmed
mkdir -p results/hisat2_aligned results/stringtie_output results/figures
mkdir -p logs

python3 run.py --check-structure   # All should show PRESENT
```

> ⚠️ **PROBLEM ENCOUNTERED:** `--check-structure` showed all directories MISSING even though they existed. Root cause: `data_loader.py` had hardcoded absolute paths.
>
> ✅ **FIX:** `config.py` now uses `Path(__file__).parent.resolve()`. If you still see MISSING, verify `config.py` has no hardcoded paths.

---

### Step 10 — Test Python Configuration

```bash
python3 scripts/config.py
# Expected:
# PROJECT CONFIGURATION
# Directories:
#   Project: /mnt/neil/resurrection_lncrna_pipeline
#   Data:    .../data
#   Results: .../results
# Directories created/verified.
```

---

<a name="section-4"></a>
## Section 4: Downloading RNA-Seq Data from NCBI

> **Script:** `bash setup/01_download.sh`
>
> Downloads all 10 samples, reference FASTA, builds HISAT2 index, prepares GTF.
> Safe to re-run — skips completed downloads using `-s` (size > 0) file checks.

**GEO accession:** GSE157098 | **SRA project:** PRJNA660052 | **~40 GB total**

### Sample Manifest

| SRA Accession | Sample Name | Condition |
|---------------|-------------|-----------|
| SRR12542161 | cp_hyd_r1 | Hydrated replicate 1 |
| SRR12542166 | cp_hyd_r3 | Hydrated replicate 3 |
| SRR12542171 | cp_rehyd_r1 | Rehydrated replicate 1 |
| SRR12542176 | cp_rehyd_r2 | Rehydrated replicate 2 |
| SRR12542181 | cp_rehyd_r3 | Rehydrated replicate 3 |
| SRR12542186 | cp_wc2_r2 | Dehydrated WC2% replicate 2 |
| SRR12542191 | cp_wc2_r3 | Dehydrated WC2% replicate 3 |
| SRR12542196 | cp_wc60_r1 | Dehydrated WC60% replicate 1 |
| SRR12542206 | cp_wc60_r2 | Dehydrated WC60% replicate 2 |
| SRR12542211 | cp_wc60_r3 | Dehydrated WC60% replicate 3 |

---

### Step 11 — Upgrade or Bypass the SRA Toolkit

> `01_download.sh` handles this automatically with a `fasterq-dump` → `wget` fallback chain. If you need to fix it manually:

```bash
# Try conda upgrade first
conda install -c bioconda sra-tools=3.1.1 -y
fasterq-dump --version   # Expected: 3.x.x

# If conda fails (missing ossuuid), download binary directly:
cd ~
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
echo 'export PATH=~/sratoolkit.3.2.1-ubuntu64/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

vdb-config --set /tls/allow-all-certs=true

# Test
fasterq-dump SRR12542206 --outdir /tmp/ --split-files --progress
# Expected: join 100%, concat 100%, spots read: ~1,198,909
rm /tmp/SRR12542206_*.fastq
```

> ⚠️ **KNOWN FAILURE:** Even after upgrade, `prefetch` may download `.sra` files that `fasterq-dump` cannot convert (`"column undefined while opening cursor"`). `01_download.sh` automatically falls back to direct `wget` from NCBI if this happens.

---

### Step 12 — Run the Download Script

```bash
bash setup/01_download.sh /mnt/neil/resurrection_lncrna_pipeline

# Or as part of full pipeline:
bash setup/run_pipeline.sh
```

**What `01_download.sh` does:**
1. Downloads all 10 FASTQ pairs via `fasterq-dump` (falls back to `wget` if needed)
2. Compresses to `.fastq.gz`
3. Downloads transcriptome FASTA from GEO
4. Builds HISAT2 index (`references/Cp_transcriptome_index`)
5. Converts annotation GFF3 → GTF

---

### Step 13 — Verify Downloads

```bash
find data/fastq_raw/ -name '*.fastq.gz' -size 0    # Expected: no output
ls data/fastq_raw/*.fastq.gz | wc -l               # Expected: 20
ls -lh data/fastq_raw/ | grep fastq.gz             # All files should be >10 MB
```

> ⚠️ **PROBLEM ENCOUNTERED:** Batch scripts created 0-byte files when downloads failed silently. `fasterq-dump` returned exit code 0 but produced empty output.
>
> ✅ **FIX:** `01_download.sh` uses `[ -s "$file" ]` (size > 0) not `[ -f "$file" ]` (exists). If 0-byte files are found, delete them and re-run — completed samples are skipped automatically.

---

<a name="section-5"></a>
## Section 5: Reference Genome & Index

> Handled automatically by `bash setup/01_download.sh`. Manual steps below for reference.

### Step 14 — Download Reference FASTA

```bash
cd references/
wget --no-check-certificate \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_transcriptome_assembly_V2.fa.gz'
gunzip GSE157098_Cp_transcriptome_assembly_V2.fa.gz
grep -c '>' GSE157098_Cp_transcriptome_assembly_V2.fa   # Expected: ~48045
```

---

### Step 15 — Build HISAT2 Index

```bash
hisat2-build \
    references/GSE157098_Cp_transcriptome_assembly_V2.fa \
    references/Cp_transcriptome_index \
    -p 8
ls references/Cp_transcriptome_index*.ht2 | wc -l   # Expected: 8
```

---

### Step 16 — Prepare GTF Annotation

```bash
# If annotation.gtf doesn't exist, convert from gff3:
conda install -c bioconda gffread -y
gffread references/annotation.gff3 -T -o references/annotation.gtf
wc -l references/annotation.gtf   # Expected: >50,000 lines
```

---

<a name="section-6"></a>
## Section 6: Running the Main Pipeline

> **Scripts:** `02_trim.sh` → `03_align.sh` → `04_stringtie.sh`
> Run individually or via the master: `bash setup/run_pipeline.sh --skip-download`

### Step 17 — Quality Trimming (`02_trim.sh`)

**Estimated time: 2–3 hours for 10 samples.**

```bash
bash setup/02_trim.sh /mnt/neil/resurrection_lncrna_pipeline
# Monitor:
tail -f logs/02_trim.log
```

**What to expect:**
```
[2026-02-21 02:18:40] Trimming cp_hyd_r1...
[3.2s][warning][gc,alloc] pool-7-thread-1: Retried waiting for GCLocker...
# GCLocker warnings are harmless Java GC messages — NOT errors
[2026-02-21 02:25:10] Done: cp_hyd_r1 (40M + 40M)
```

**Verify:**
```bash
ls -lh data/fastq_trimmed/ | grep '_paired'   # Expected: 20 files, all >10 MB
```

> ⚠️ **PROBLEM ENCOUNTERED:** `set -euo pipefail` caused the pipeline to exit after the first sample when Trimmomatic's GC warning returned non-zero.
>
> ✅ **FIX:** `02_trim.sh` uses `set -uo pipefail` (no `-e`). GCLocker warnings are logged but do not abort the run. Real errors are caught by the `-s` output file check.

---

### Step 18 — HISAT2 Alignment (`03_align.sh`)

**Estimated time: 3–4 hours for 10 samples.**

```bash
bash setup/03_align.sh /mnt/neil/resurrection_lncrna_pipeline
# Monitor:
tail -f logs/03_align.log

# Check alignment rates after:
grep 'overall alignment rate' logs/03_align_*.log
# Expected: 60-95% per sample. <40% indicates a problem.
```

> ⚠️ **PROBLEM ENCOUNTERED:** `cp_hyd_r1` BAM was missing because it was processed early without the pipeline script. StringTie ran without `-G`, producing 6 lines in the abundance file instead of 48,045.
>
> ✅ **FIX:** `04_stringtie.sh` checks abundance file line count before skipping. Files with fewer than 1,000 lines are automatically re-run regardless of whether they exist.

---

### Step 19 — StringTie Assembly (`04_stringtie.sh`)

**Estimated time: 2–3 hours for 10 samples.**

```bash
bash setup/04_stringtie.sh /mnt/neil/resurrection_lncrna_pipeline
# Monitor:
tail -f logs/04_stringtie.log
```

**Critical verification:**
```bash
wc -l results/stringtie_output/*_abundance.txt | sort -n
# ALL files must have ~48046 lines
# If any file shows 7 lines: -G flag was missing. 04_stringtie.sh will auto-detect and rerun.
```

---

<a name="section-7"></a>
## Section 7: Expression Matrix, DE Analysis & Visualization

> **Scripts:** `05_merge_matrix.py` → `06_de_analysis.py` → `07_visualize.py`

### Step 20 — Verify All StringTie Outputs

```bash
ls results/stringtie_output/*_abundance.txt | wc -l   # Expected: 10
wc -l results/stringtie_output/*_abundance.txt | sort -n | head -12
# All values should be near 48046
```

---

### Step 21 — Merge TPM Matrix (`05_merge_matrix.py`)

```bash
python3 setup/05_merge_matrix.py /mnt/neil/resurrection_lncrna_pipeline

# Expected output:
# Loaded: cp_hyd_r1            (48,045 genes)
# Loaded: cp_hyd_r3            (48,045 genes)
# ... (all 10 samples)
# Saved: results/stringtie_tpm_matrix.csv
# Shape: 48,045 genes × 12 columns
```

The script exits with an error if any abundance file has fewer than 1,000 rows, preventing a silently broken matrix.

---

### Step 22 — Differential Expression (`06_de_analysis.py`)

```bash
python3 setup/06_de_analysis.py /mnt/neil/resurrection_lncrna_pipeline

# Expected output:
# lncRNA candidates: 38,978
# Expressed lncRNAs: 17,574
# P-value distribution:
#   pval < 0.05: 780
# After FDR correction:
#   padj < 0.05: 0        ← expected with n=2 hydrated replicates
# Candidates (pval<0.05, |log2FC|>2, max TPM>=5): 777
# Top 14 stress-responsive lncRNAs: [table]
```

> ⚠️ **IMPORTANT — Statistical Limitation:**
> With only **n=2 hydrated replicates**, BH FDR correction across 17,574 tests yields **0 significant genes**. This is mathematically expected — the t-test has 1 degree of freedom, insufficient for genome-wide correction.
>
> `06_de_analysis.py` selects top 14 candidates using: **raw p < 0.05 + |log2FC| > 2 + expression magnitude score**. This is defensible for a pilot study but must be disclosed in any write-up.

---

### Step 23 — Visualizations (`07_visualize.py`)

```bash
python3 setup/07_visualize.py /mnt/neil/resurrection_lncrna_pipeline

# Expected:
# Saved: results/volcano_plot.png
# Saved: results/heatmap_significant_lncrnas.png
# Saved: results/boxplots_top_lncrnas.png
ls -lh results/*.png   # Each file should be >100 KB
```

---

<a name="section-8"></a>
## Section 8: Validating Your Results

### Step 24 — Sanity Check Results

```bash
python3 -c "
import pandas as pd
tpm = pd.read_csv('results/stringtie_tpm_matrix.csv')
de  = pd.read_csv('results/lncrna_differential_expression.csv')
sig = pd.read_csv('results/significant_stress_responsive_lncrnas.csv')
print(f'TPM matrix:    {tpm.shape}')   # Expected: (48045, 12)
print(f'DE results:    {len(de)}')      # Expected: 17574
print(f'Top candidates:{len(sig)}')     # Expected: 14
print(f'Columns: {list(tpm.columns[:4])}...')
"
```

---

### Step 25 — Verify Top Candidate

```bash
python3 -c "
import pandas as pd
sig = pd.read_csv('results/significant_stress_responsive_lncrnas.csv')
top = sig.iloc[0]
print(f'Top candidate: {top.gene_id}')
print(f'Hydrated:      {top.mean_hydrated:.1f} TPM')
print(f'Dehydrated:    {top.mean_dehydrated:.1f} TPM')
print(f'Fold change:   {2**top.log2FC:.1f}x')
print(f'p-value:       {top.pvalue:.4f}')
"
# Expected: >5-fold increase in dehydrated vs hydrated
```

---

### Expected Results Summary

| Output File | Expected Shape | Key Numbers |
|-------------|---------------|-------------|
| `stringtie_tpm_matrix.csv` | 48,045 × 12 | 10 sample columns + Gene ID + Gene Name |
| `coding_potential.txt` | 48,045 rows | 38,978 lncRNAs / 9,067 coding (81%/19%) |
| `lncrna_differential_expression.csv` | 17,574 rows | All expressed lncRNAs with stats |
| `significant_stress_responsive_lncrnas.csv` | 14 rows | Top stress-induced candidates |
| `volcano_plot.png` | — | >100 KB |
| `heatmap_significant_lncrnas.png` | — | >100 KB |
| `boxplots_top_lncrnas.png` | — | >100 KB |

---

<a name="section-9"></a>
## Section 9: Troubleshooting Reference

| Error / Symptom | Root Cause | Script That Handles It | Manual Fix |
|----------------|------------|----------------------|------------|
| `mbedtls_ssl_handshake returned -9984` | SRA Toolkit 2.9.6 incompatible TLS | `01_download.sh` (wget fallback) | Upgrade SRA Toolkit to 3.x |
| `invalid accession 'SRR...'` | Same TLS failure, fasterq-dump can't reach NCBI | `01_download.sh` | Upgrade SRA Toolkit first |
| All directories MISSING in `--check-structure` | Hardcoded paths in data_loader.py | `00_preflight.sh` (auto-creates dirs) | Use `Path(__file__).parent.resolve()` |
| `LibMambaUnsatisfiableError: ossuuid` | sra-tools 3.x needs ossuuid not in conda | — | Download SRA toolkit binary from NCBI |
| `LibMambaUnsatisfiableError` (hisat2/samtools) | Conda dependency conflicts | — | `sudo apt-get install hisat2 samtools` |
| 0-byte FASTQ files | Download failed silently, left empty file | `01_download.sh` uses `-s` check | Delete 0-byte files, re-run `01_download.sh` |
| `_abundance.txt` has 7 lines | StringTie ran without `-G` flag | `04_stringtie.sh` auto-detects & reruns | Rerun with `-G annotation.gtf -e` |
| `padj < 0.05: 0 genes` | n=2 hydrated replicates — insufficient FDR power | `06_de_analysis.py` uses raw p-value fallback | Expected — document as limitation |
| GCLocker warnings (Trimmomatic) | Harmless Java GC messages | `02_trim.sh` uses `set -uo` not `set -euo` | Ignore |
| Pipeline exits after first sample | `set -euo pipefail` + GCLocker non-zero return | Fixed in `02_trim.sh` | `sed -i 's/set -euo/set -uo/' script.sh` |
| EBI FTP paths return 404 | EBI mirror path structure varies by accession | `01_download.sh` uses NCBI direct | Use fasterq-dump or NCBI wget |

---

<a name="section-10"></a>
## Section 10: Quick Reference — All 27 Steps

| Step | Action | Script | Est. Time |
|------|--------|--------|-----------|
| 1 | Check OS, disk, RAM | `00_preflight.sh` | 1 min |
| 2 | Check Python 3.7+ | `00_preflight.sh` | auto |
| 3 | Check / install Git | `00_preflight.sh` | auto |
| 4 | Fix SSL certificates | `00_preflight.sh` (check) / manual fix | 5 min |
| 5 | Install Miniconda | manual | 10 min |
| 6 | Install pipeline tools | manual (conda + apt-get) | 15 min |
| 7 | Verify all tools | `00_preflight.sh` | auto |
| 8 | Clone repository | `git clone` | 2 min |
| 9 | Create gitignored directories | `00_preflight.sh` (auto-creates) | auto |
| 10 | Test Python config | `python3 scripts/config.py` | 1 min |
| 11 | Upgrade SRA Toolkit | manual / `01_download.sh` handles fallback | 10 min |
| 12 | Download 10 samples | `01_download.sh` | 2–4 hrs |
| 13 | Verify downloads (no 0-byte) | `01_download.sh` (built-in check) | auto |
| 14 | Download reference FASTA | `01_download.sh` | 20 min |
| 15 | Build HISAT2 index | `01_download.sh` | 15 min |
| 16 | Prepare GTF annotation | `01_download.sh` | 5 min |
| 17 | Trimmomatic — all 10 samples | `02_trim.sh` | 2–3 hrs |
| 18 | HISAT2 alignment — all 10 | `03_align.sh` | 3–4 hrs |
| 19 | StringTie assembly — all 10 | `04_stringtie.sh` | 2–3 hrs |
| 20 | Verify StringTie outputs | `04_stringtie.sh` (built-in check) | auto |
| 21 | Merge TPM expression matrix | `05_merge_matrix.py` | 5 min |
| 22 | Differential expression | `06_de_analysis.py` | 5 min |
| 23 | Visualizations | `07_visualize.py` | 5 min |
| 24 | Sanity check results | manual `python3 -c '...'` | 5 min |
| 25 | Verify top candidate | manual `python3 -c '...'` | 2 min |
| 26 | Co-expression network *(next)* | `scripts/task7_coexpression.py` | 30 min |
| 27 | GO enrichment *(next)* | `scripts/task8_go_enrichment.py` | 30 min |

---

> **Total Pipeline Runtime: 8–12 hours end-to-end**
>
> Steps 1–10: ~1 hour setup | Steps 11–16: ~3 hours data prep | Steps 17–25: ~8 hours compute

---

*Resurrection Plant lncRNA Pipeline · Setup Guide*
*Built from actual execution history — all problems encountered and resolved are documented above*
