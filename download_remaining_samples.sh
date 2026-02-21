#!/bin/bash

# Array of remaining SRR accessions and sample names
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

echo "Downloading remaining 9 samples..."
echo ""

for ACC in "${!SAMPLES[@]}"; do
    NAME="${SAMPLES[$ACC]}"
    echo "Downloading $NAME ($ACC)..."
    
    # NCBI path
    PREFIX=$(echo $ACC | cut -c1-6)
    
    wget --no-check-certificate https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/${PREFIX:0:3}/${PREFIX:3:3}/${ACC}/${ACC}.lite.1 -O data/fastq_raw/${NAME}.sra 2>&1 | grep -E "saved|failed"
done

echo ""
echo "Downloads complete!"
ls -lh data/fastq_raw/ | grep sra
