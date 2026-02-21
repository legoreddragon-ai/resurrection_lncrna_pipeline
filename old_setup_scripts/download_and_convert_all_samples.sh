#!/bin/bash
echo "Downloading and converting all 9 remaining samples..."
echo ""

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

mkdir -p data/fastq_raw

for ACC in "${!SAMPLES[@]}"; do
    NAME="${SAMPLES[$ACC]}"
    echo "Downloading and converting $NAME ($ACC)..."

    fasterq-dump "$ACC" \
        --outdir data/fastq_raw/ \
        --outfile "${NAME}.fastq" \
        --split-files \
        --threads 4 \
        --progress

    if [ $? -ne 0 ]; then
        echo "ERROR: fasterq-dump failed for $ACC, skipping..."
        continue
    fi

    echo "Compressing FASTQ files for $NAME..."
    gzip data/fastq_raw/${NAME}_1.fastq
    gzip data/fastq_raw/${NAME}_2.fastq

    echo "Complete: $NAME"
    echo ""
done

echo "All downloads and conversions complete!"
ls -lh data/fastq_raw/ | grep cp_
