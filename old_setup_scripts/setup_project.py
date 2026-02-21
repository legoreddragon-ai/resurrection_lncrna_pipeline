#!/usr/bin/env python3
"""
Setup script to download all necessary files and build project structure
for Resurrection Plant lncRNA Discovery Pipeline
"""

import os
import urllib.request
import subprocess
import sys
from pathlib import Path

def create_directories():
    """Create all necessary directories"""
    dirs = [
        'data/fastq_raw',
        'data/fastq_trimmed',
        'references',
        'results/stringtie_output',
        'results/hisat2_aligned',
        'scripts',
        'docs'
    ]
    
    for dir_path in dirs:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    print("Created directory structure")

def download_annotation():
    """Download annotation file from GEO"""
    print("\nDownloading annotation file...")
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz"
    output = "references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt.gz"
    
    try:
        urllib.request.urlretrieve(url, output)
        print(f"Downloaded: {output}")
        
        # Decompress
        subprocess.run(['gunzip', '-f', output], check=True)
        print("Decompressed annotation file")
    except Exception as e:
        print(f"Error downloading annotation: {e}")
        print("Download manually from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157098")

def download_transcriptome():
    """Download transcriptome FASTA from GEO"""
    print("\nDownloading transcriptome FASTA...")
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157098/suppl/GSE157098_Cp_transcriptome_assembly_V2.fa.gz"
    output = "references/GSE157098_Cp_transcriptome_assembly_V2.fa.gz"
    
    try:
        urllib.request.urlretrieve(url, output)
        print(f"Downloaded: {output}")
        
        # Decompress
        subprocess.run(['gunzip', '-f', output], check=True)
        print("Decompressed transcriptome file")
    except Exception as e:
        print(f"Error downloading transcriptome: {e}")
        print("Download manually from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157098")

def download_fastq_files():
    """
    Guide user to download FASTQ files from SRA
    """
    print("\nFASTQ FILES - Download Instructions")
    print("="*60)
    print("You need to download 60 FASTQ files from SRA")
    print("Project ID: PRJNA660052")
    print("Series ID: GSE157098")
    print("\nOptions:")
    print("1. Use SRA Toolkit:")
    print("   conda install -c bioconda sra-tools")
    print("   prefetch SRR12345678 (for each run)")
    print("   fasterq-dump SRR12345678")
    print("\n2. Download from NCBI directly:")
    print("   https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA660052")
    print("\n3. Use Aspera (fastest):")
    print("   ascp -i ~/aspera-key.openssh -Tr -Q -l 300m ...")
    print("\nPlace all FASTQ files in: data/fastq_raw/")
    print("="*60)

def create_annotation_files():
    """Convert TSV annotation to GFF3 and GTF"""
    print("\nConverting annotation files...")
    
    try:
        import pandas as pd
        
        df = pd.read_csv('references/GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt', sep='\t')
        
        # Create GFF3
        with open('references/annotation.gff3', 'w') as out:
            out.write("##gff-version 3\n")
            for idx, row in df.iterrows():
                transcript_id = row['Feature ID']
                annotation = str(row.get('Annotations (viridiplantae)', 'unknown')).replace('"', '')
                best_hit = str(row.get('Best hit (viridiplantae)', '')).replace('"', '')[:50]
                
                attributes = f"ID={transcript_id};Name={transcript_id};product={annotation};note={best_hit}"
                out.write(f"{transcript_id}\tGEO\tmRNA\t1\t1000\t.\t+\t.\t{attributes}\n")
        
        print("Created: references/annotation.gff3")
        
        # Create GTF
        with open('references/annotation.gtf', 'w') as out:
            for idx, row in df.iterrows():
                transcript_id = row['Feature ID']
                annotation = str(row.get('Annotations (viridiplantae)', 'unknown')).replace('"', '')
                
                attributes = f'gene_id "{transcript_id}"; transcript_id "{transcript_id}"; product "{annotation}";'
                out.write(f"{transcript_id}\tGEO\tmRNA\t1\t1000\t.\t+\t.\t{attributes}\n")
        
        print("Created: references/annotation.gtf")
        
    except Exception as e:
        print(f"Error creating annotation files: {e}")

def check_dependencies():
    """Check for required bioinformatics tools"""
    print("\nChecking dependencies...")
    
    tools = ['hisat2', 'samtools', 'stringtie', 'trimmomatic']
    missing = []
    
    for tool in tools:
        try:
            subprocess.run([tool, '--version'], capture_output=True, check=True)
            print(f"Found: {tool}")
        except:
            missing.append(tool)
            print(f"Missing: {tool}")
    
    if missing:
        print(f"\nInstall missing tools with:")
        print(f"conda install -c bioconda {' '.join(missing)}")
    
    return len(missing) == 0

def main():
    print("="*60)
    print("RESURRECTION PLANT lncRNA PIPELINE - PROJECT SETUP")
    print("="*60)
    
    create_directories()
    download_annotation()
    download_transcriptome()
    download_fastq_files()
    create_annotation_files()
    has_deps = check_dependencies()
    
    print("\n" + "="*60)
    print("SETUP COMPLETE")
    print("="*60)
    print("\nNext steps:")
    print("1. Download FASTQ files from SRA (see instructions above)")
    print("2. Run: python3 pipeline.py")
    print("="*60)

if __name__ == '__main__':
    main()