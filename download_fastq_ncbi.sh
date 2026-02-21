#!/bin/bash

# Array of SRR accessions
ACCESSIONS=(
    "SRR12542161"
    "SRR12542166"
    "SRR12542171"
    "SRR12542176"
    "SRR12542181"
    "SRR12542186"
    "SRR12542191"
    "SRR12542196"
    "SRR12542206"
    "SRR12542211"
)

NAMES=(
    "cp_hyd_r1"
    "cp_hyd_r3"
    "cp_rehyd_r1"
    "cp_rehyd_r2"
    "cp_rehyd_r3"
    "cp_wc2_r2"
    "cp_wc2_r3"
    "cp_wc60_r1"
    "cp_wc60_r2"
    "cp_wc60_r3"
)

for i in ${!ACCESSIONS[@]}; do
    ACC=${ACCESSIONS[$i]}
    NAME=${NAMES[$i]}
    
    echo "Downloading $NAME ($ACC)..."
    
    # NCBI direct download path
    PREFIX=$(echo $ACC | cut -c1-6)
    wget --no-check-certificate https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/${PREFIX:0:3}/${PREFIX:3:3}/${ACC}/${ACC}_1.fastq.gz -O data/fastq_raw/${NAME}_1.fastq.gz
    wget --no-check-certificate https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-924/${PREFIX:0:3}/${PREFIX:3:3}/${ACC}/${ACC}_2.fastq.gz -O data/fastq_raw/${NAME}_2.fastq.gz
done

echo "Done!"
ls -lh data/fastq_raw/ | grep cp_
