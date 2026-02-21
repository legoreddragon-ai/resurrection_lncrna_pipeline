#!/bin/bash

echo "Downloading 10 samples from EBI mirror..."
echo ""

# SRR to sample name mapping with correct EBI paths
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
    echo "Downloading $NAME ($ACC)..."
    
    # EBI path format: /vol1/fastq/SRR125/[dir]/[SRR]/[SRR]_[1/2].fastq.gz
    DIR=$(printf "%03d" $(((${ACC: -3:3} - 161) / 100 + 1)))
    
    wget --no-check-certificate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/${DIR}/${ACC}/${ACC}_1.fastq.gz -O data/fastq_raw/${NAME}_1.fastq.gz 2>&1 | grep -v "^--"
    wget --no-check-certificate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/${DIR}/${ACC}/${ACC}_2.fastq.gz -O data/fastq_raw/${NAME}_2.fastq.gz 2>&1 | grep -v "^--"
done

echo "Download complete!"
ls -lh data/fastq_raw/ | grep -v "^total"
