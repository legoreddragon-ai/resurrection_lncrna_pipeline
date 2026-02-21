#!/bin/bash

set -e

PROJECT_DIR="$HOME/resurrection_lncrna_pipeline"
cd "$PROJECT_DIR"

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}=========================================================================${NC}"
echo -e "${BLUE}RESURRECTION PLANT lncRNA DISCOVERY PIPELINE${NC}"
echo -e "${BLUE}=========================================================================${NC}"

echo -e "\n${YELLOW}[STEP 0] Creating directories...${NC}"
mkdir -p data/fastq_raw data/fastq_trimmed references results/stringtie_output results/hisat2_aligned scripts
echo -e "${GREEN}Directories created${NC}"

echo -e "\n${YELLOW}[STEP 1] Downloading annotation file...${NC}"

if [ ! -f "references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt" ]; then
    echo "Downloading annotation from GEO..."
    wget -q ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz \
        -O references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz
    gunzip -f references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz
    echo -e "${GREEN}Annotation downloaded${NC}"
else
    echo -e "${GREEN}Annotation already exists${NC}"
fi

echo -e "\n${YELLOW}[STEP 2] Converting annotation to GFF3/GTF...${NC}"

if [ ! -f "references/annotation.gff3" ]; then
    python3 << 'PYEOF'
import pandas as pd

df = pd.read_csv('references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt', sep='\t')

with open('references/annotation.gff3', 'w') as out:
    out.write("##gff-version 3\n")
    for idx, row in df.iterrows():
        transcript_id = row['Feature ID']
        annotation = str(row.get('Annotations (viridiplantae)', 'unknown')).replace('"', '')
        best_hit = str(row.get('Best hit (viridiplantae)', '')).replace('"', '')[:50]
        
        attributes = f"ID={transcript_id};Name={transcript_id};product={annotation};note={best_hit}"
        out.write(f"{transcript_id}\tGEO\tmRNA\t1\t1000\t.\t+\t.\t{attributes}\n")

with open('references/annotation.gtf', 'w') as out:
    for idx, row in df.iterrows():
        transcript_id = row['Feature ID']
        annotation = str(row.get('Annotations (viridiplantae)', 'unknown')).replace('"', '')
        
        attributes = f'gene_id "{transcript_id}"; transcript_id "{transcript_id}"; product "{annotation}";'
        out.write(f"{transcript_id}\tGEO\tmRNA\t1\t1000\t.\t+\t.\t{attributes}\n")

print(f"Converted {len(df)} transcripts to GFF3 and GTF")
PYEOF
    echo -e "${GREEN}Annotation converted${NC}"
else
    echo -e "${GREEN}GFF3/GTF already exist${NC}"
fi

echo -e "\n${BLUE}=========================================================================${NC}"
echo -e "${GREEN}Pipeline setup complete${NC}"
echo -e "${BLUE}=========================================================================${NC}"
echo ""
echo "Next steps:"
echo "1. Download FASTQ files from SRA (PRJNA660052)"
echo "2. Place trimmed FASTQ files in data/fastq_trimmed/"
echo "3. Run alignment and assembly steps"
echo ""

