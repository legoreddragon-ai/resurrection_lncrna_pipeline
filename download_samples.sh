#!/bin/bash

echo "Downloading 10 samples (1 technical replicate each) from EBI mirror..."
echo ""

# One representative SRR per GSM sample
SAMPLES=(
    "SRR12542161:cp_hyd_r1"
    "SRR12542166:cp_hyd_r3"
    "SRR12542171:cp_rehyd_r1"
    "SRR12542176:cp_rehyd_r2"
    "SRR12542181:cp_rehyd_r3"
    "SRR12542186:cp_wc2_r2"
    "SRR12542191:cp_wc2_r3"
    "SRR12542196:cp_wc60_r1"
    "SRR12542206:cp_wc60_r2"
    "SRR12542211:cp_wc60_r3"
)

for SAMPLE in "${SAMPLES[@]}"; do
    ACC="${SAMPLE%:*}"
    NAME="${SAMPLE#*:}"
    PREFIX=$(printf "%03d" $((${ACC: -3} / 100)))
    
    echo "Downloading $NAME ($ACC)..."
    wget --no-check-certificate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/${PREFIX}/${ACC}/${ACC}_1.fastq.gz -O data/fastq_raw/${NAME}_1.fastq.gz
    wget --no-check-certificate ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR125/${PREFIX}/${ACC}/${ACC}_2.fastq.gz -O data/fastq_raw/${NAME}_2.fastq.gz
    echo ""
done

echo "All downloads complete!"
ls -lh data/fastq_raw/
