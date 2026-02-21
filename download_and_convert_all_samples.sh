#!/bin/bash

echo "Downloading and converting all 9 remaining samples..."
echo ""

# Array of SRR accessions and sample names
declare -A SAMPLES=(
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
    
    wget --no-check-certificate https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/SRR012/12542/${ACC}/${ACC}.lite.1 -O data/fastq_raw/${NAME}.sra
    
    echo "Converting $NAME to FASTQ..."
    fasterq-dump data/fastq_raw/${NAME}.sra -O data/fastq_raw/
    
    echo "Compressing FASTQ files for $NAME..."
    gzip data/fastq_raw/${NAME}_1.fastq
    gzip data/fastq_raw/${NAME}_2.fastq
    
    echo "Cleaning up SRA file..."
    rm data/fastq_raw/${NAME}.sra
    
    echo "Complete: $NAME"
    echo ""
done

echo "All downloads and conversions complete!"
ls -lh data/fastq_raw/ | grep cp_
