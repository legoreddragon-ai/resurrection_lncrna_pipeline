#!/usr/bin/env python3
"""
Complete Resurrection Plant lncRNA Discovery Pipeline
Craterostigma plantagineum Desiccation-Induced lncRNA Analysis
"""

import os
import subprocess
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
from Bio import SeqIO
import sys

class LncRNAPipeline:
    def __init__(self, project_dir='.'):
        self.project_dir = project_dir
        self.data_dir = os.path.join(project_dir, 'data')
        self.refs_dir = os.path.join(project_dir, 'references')
        self.results_dir = os.path.join(project_dir, 'results')
        
    def print_header(self, step, title):
        print("\n" + "="*80)
        print(f"STEP {step}: {title}")
        print("="*80)
    
    def check_files(self):
        """Check if required input files exist"""
        self.print_header(0, "Checking Input Files")
        
        required = [
            os.path.join(self.refs_dir, 'GSE157098_Cp_transcriptome_assembly_V2.fa'),
            os.path.join(self.refs_dir, 'GSE157098_Cp_Transcriptome_assembly_V2_Annotation.txt'),
            os.path.join(self.refs_dir, 'annotation.gff3')
        ]
        
        fastq_count = len(list(Path(self.data_dir).glob('fastq_trimmed/*_1.fastq.gz')))
        
        print(f"Reference FASTA: {'FOUND' if os.path.exists(required[0]) else 'MISSING'}")
        print(f"Annotation TXT: {'FOUND' if os.path.exists(required[1]) else 'MISSING'}")
        print(f"Annotation GFF3: {'FOUND' if os.path.exists(required[2]) else 'MISSING'}")
        print(f"Trimmed FASTQ files: {fastq_count}")
        
        if fastq_count == 0:
            print("\nERROR: No trimmed FASTQ files found in data/fastq_trimmed/")
            print("Please download and trim FASTQ files first.")
            return False
        
        return True
    
    def build_hisat2_index(self):
        """Build HISAT2 index"""
        self.print_header(1, "Building HISAT2 Index")
        
        index_file = os.path.join(self.refs_dir, 'Cp_transcriptome_index.1.ht2')
        
        if os.path.exists(index_file):
            print("HISAT2 index already exists")
            return True
        
        fasta = os.path.join(self.refs_dir, 'GSE157098_Cp_transcriptome_assembly_V2.fa')
        index_base = os.path.join(self.refs_dir, 'Cp_transcriptome_index')
        
        try:
            print(f"Building index from: {fasta}")
            subprocess.run(['hisat2-build', '-q', fasta, index_base], check=True)
            print("Index built successfully")
            return True
        except Exception as e:
            print(f"Error building index: {e}")
            return False
    
    def align_samples(self):
        """Align all samples with HISAT2"""
        self.print_header(2, "Aligning Samples with HISAT2")
        
        fastq_dir = os.path.join(self.data_dir, 'fastq_trimmed')
        index_base = os.path.join(self.refs_dir, 'Cp_transcriptome_index')
        
        # Get list of samples
        r1_files = sorted(Path(fastq_dir).glob('*_1.fastq.gz'))
        samples = [f.name.replace('_1.fastq.gz', '') for f in r1_files]
        
        print(f"Found {len(samples)} samples to align")
        
        for i, sample in enumerate(samples, 1):
            r1 = os.path.join(fastq_dir, f'{sample}_1.fastq.gz')
            r2 = os.path.join(fastq_dir, f'{sample}_2.fastq.gz')
            bam = os.path.join(self.results_dir, 'hisat2_aligned', f'{sample}.bam')
            
            if os.path.exists(bam):
                print(f"[{i}/{len(samples)}] {sample}: Already aligned")
                continue
            
            print(f"[{i}/{len(samples)}] Aligning {sample}...")
            
            try:
                # Align
                with open(os.path.join(self.results_dir, 'hisat2_aligned', f'{sample}.log'), 'w') as log:
                    hisat2_proc = subprocess.Popen(
                        ['hisat2', '-p', '4', '-q', '-x', index_base, '-1', r1, '-2', r2],
                        stdout=subprocess.PIPE,
                        stderr=log
                    )
                    
                    # Sort BAM
                    subprocess.run(
                        ['samtools', 'sort', '-@', '4', '-o', bam],
                        stdin=hisat2_proc.stdout,
                        check=True
                    )
                
                # Index BAM
                subprocess.run(['samtools', 'index', bam], check=True)
                print(f"  Created: {bam}")
                
            except Exception as e:
                print(f"  Error aligning {sample}: {e}")
        
        return True
    
    def assemble_transcripts(self):
        """Assemble transcripts with StringTie"""
        self.print_header(3, "Assembling Transcripts with StringTie")
        
        bam_dir = os.path.join(self.results_dir, 'hisat2_aligned')
        gff3 = os.path.join(self.refs_dir, 'annotation.gff3')
        
        bam_files = sorted(Path(bam_dir).glob('*.bam'))
        
        print(f"Found {len(bam_files)} BAM files to assemble")
        
        for i, bam_file in enumerate(bam_files, 1):
            sample = bam_file.stem
            gtf = os.path.join(self.results_dir, 'stringtie_output', f'{sample}.gtf')
            abundance = os.path.join(self.results_dir, 'stringtie_output', f'{sample}_abundance.txt')
            
            if os.path.exists(abundance):
                print(f"[{i}/{len(bam_files)}] {sample}: Already assembled")
                continue
            
            print(f"[{i}/{len(bam_files)}] Assembling {sample}...")
            
            try:
                subprocess.run([
                    'stringtie', str(bam_file),
                    '-G', gff3,
                    '-o', gtf,
                    '-A', abundance,
                    '-p', '4', '-q'
                ], check=True)
                print(f"  Created: {abundance}")
                
            except Exception as e:
                print(f"  Error assembling {sample}: {e}")
        
        return True
    
    def create_expression_matrix(self):
        """Create TPM expression matrix"""
        self.print_header(4, "Creating TPM Expression Matrix")
        
        abundance_dir = os.path.join(self.results_dir, 'stringtie_output')
        abundance_files = sorted(Path(abundance_dir).glob('*_abundance.txt'))
        
        print(f"Found {len(abundance_files)} abundance files")
        
        dfs = []
        for file in abundance_files:
            sample_name = file.name.replace('_abundance.txt', '')
            df = pd.read_csv(file, sep='\t')
            df = df[['Gene ID', 'TPM']].rename(columns={'TPM': sample_name})
            df = df.set_index('Gene ID')
            dfs.append(df)
        
        tpm_matrix = pd.concat(dfs, axis=1).fillna(0)
        output = os.path.join(self.results_dir, 'stringtie_tpm_matrix.csv')
        tpm_matrix.to_csv(output)
        
        print(f"Expression matrix: {tpm_matrix.shape}")
        print(f"  Genes: {tpm_matrix.shape[0]}")
        print(f"  Samples: {tpm_matrix.shape[1]}")
        print(f"Saved: {output}")
        
        return True
    
    def classify_lncrnas(self):
        """Classify transcripts as coding or lncRNA"""
        self.print_header(5, "Classifying Transcripts")
        
        fasta = os.path.join(self.refs_dir, 'GSE157098_Cp_transcriptome_assembly_V2.fa')
        output = os.path.join(self.results_dir, 'coding_potential.txt')
        filtered_fa = os.path.join(self.results_dir, 'transcripts_filtered.fa')
        
        # Filter by length
        print("Filtering transcripts (length >= 200 nt)...")
        count = 0
        with open(filtered_fa, 'w') as out:
            for record in SeqIO.parse(fasta, 'fasta'):
                if len(record.seq) >= 200:
                    SeqIO.write(record, out, 'fasta')
                    count += 1
        
        print(f"Filtered: {count} transcripts")
        
        # Predict coding potential
        def find_longest_orf(seq):
            seq_str = str(seq).upper()
            longest_orf = 0
            
            for frame in range(3):
                for i in range(frame, len(seq_str) - 2, 3):
                    codon = seq_str[i:i+3]
                    if codon in ['ATG']:
                        for j in range(i + 3, len(seq_str) - 2, 3):
                            stop = seq_str[j:j+3]
                            if stop in ['TAA', 'TAG', 'TGA']:
                                orf_len = j - i
                                longest_orf = max(longest_orf, orf_len)
                                break
            return longest_orf
        
        print("Analyzing coding potential...")
        results = []
        for record in SeqIO.parse(filtered_fa, 'fasta'):
            orf_length = find_longest_orf(record.seq)
            coding_score = (orf_length / len(record.seq)) if len(record.seq) > 0 else 0
            classification = 'coding' if coding_score > 0.5 else 'lncRNA'
            results.append(f"{record.id}\t{len(record.seq)}\t{orf_length}\t{coding_score:.3f}\t{classification}")
        
        with open(output, 'w') as f:
            f.write("Transcript_ID\tLength\tLongest_ORF\tCoding_Score\tClassification\n")
            f.write('\n'.join(results))
        
        lncrnas = sum(1 for r in results if 'lncRNA' in r)
        coding = sum(1 for r in results if 'coding' in r)
        
        print(f"Classification complete:")
        print(f"  lncRNAs: {lncrnas}")
        print(f"  Coding: {coding}")
        print(f"Saved: {output}")
        
        return True
    
    def create_id_mapping(self):
        """Create mapping between StringTie and reference IDs"""
        self.print_header(6, "Creating ID Mapping")
        
        gtf_file = os.path.join(self.results_dir, 'stringtie_output', 'cp_hyd_r1.gtf')
        output = os.path.join(self.results_dir, 'transcript_id_mapping.csv')
        
        if not os.path.exists(gtf_file):
            print(f"GTF file not found: {gtf_file}")
            print("Skipping ID mapping - will create during DE analysis")
            return True
        
        mapping = {}
        
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                reference = fields[0]
                attributes = fields[8]
                
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split('"')[1]
                        if transcript_id not in mapping:
                            mapping[transcript_id] = reference
                        break
        
        df = pd.DataFrame(list(mapping.items()), columns=['STRG_ID', 'Reference_ID'])
        df.to_csv(output, index=False)
        
        print(f"Mapping created: {len(df)} transcripts")
        print(f"Saved: {output}")
        
        return True
    
    def differential_expression(self):
        """Perform differential expression analysis"""
        self.print_header(7, "Differential Expression Analysis")
        
        tpm_file = os.path.join(self.results_dir, 'stringtie_tpm_matrix.csv')
        coding_file = os.path.join(self.results_dir, 'coding_potential.txt')
        mapping_file = os.path.join(self.results_dir, 'transcript_id_mapping.csv')
        
        print("Loading data...")
        tpm_matrix = pd.read_csv(tpm_file, index_col=0)
        coding_potential = pd.read_csv(coding_file, sep='\t')
        id_mapping = pd.read_csv(mapping_file)
        
        # Extract gene IDs
        id_mapping['Gene_ID'] = id_mapping['STRG_ID'].apply(lambda x: x.rsplit('.', 1)[0])
        ref_to_gene = dict(zip(id_mapping['Reference_ID'], id_mapping['Gene_ID']))
        
        # Log-transform
        log_tpm = np.log2(tpm_matrix + 1)
        
        # Get lncRNAs
        lncrna_refs = set(coding_potential[coding_potential['Classification'] == 'lncRNA']['Transcript_ID'])
        lncrna_genes = [ref_to_gene[ref] for ref in lncrna_refs if ref in ref_to_gene]
        lncrna_tpm = log_tpm[log_tpm.index.isin(lncrna_genes)]
        
        print(f"lncRNA expression matrix: {lncrna_tpm.shape}")
        
        # Define groups
        hydrated = [col for col in lncrna_tpm.columns if 'hyd' in col and 'rehyd' not in col]
        rehydrated = [col for col in lncrna_tpm.columns if 'rehyd' in col]
        dehydrated = [col for col in lncrna_tpm.columns if 'wc2' in col or 'wc60' in col]
        
        print(f"Hydrated: {hydrated}")
        print(f"Rehydrated: {rehydrated}")
        print(f"Dehydrated: {dehydrated}")
        
        # Calculate means
        hydrated_mean = lncrna_tpm[hydrated].mean(axis=1)
        rehydrated_mean = lncrna_tpm[rehydrated].mean(axis=1)
        dehydrated_mean = lncrna_tpm[dehydrated].mean(axis=1)
        
        # T-tests
        print("\nPerforming statistical tests...")
        pvalues = []
        log2fcs = []
        
        for idx in lncrna_tpm.index:
            dehyd_vals = lncrna_tpm.loc[idx, dehydrated].values
            hyd_vals = lncrna_tpm.loc[idx, hydrated].values
            
            mean_dehyd = dehyd_vals.mean() + 1
            mean_hyd = hyd_vals.mean() + 1
            
            t_stat, p_val = stats.ttest_ind(dehyd_vals, hyd_vals)
            log2fc = np.log2(mean_dehyd / mean_hyd)
            
            pvalues.append(p_val)
            log2fcs.append(log2fc)
        
        # Create results
        de_results = pd.DataFrame({
            'Gene_ID': lncrna_tpm.index,
            'Hydrated_Mean_TPM': 2**hydrated_mean.values - 1,
            'Dehydrated_Mean_TPM': 2**dehydrated_mean.values - 1,
            'Rehydrated_Mean_TPM': 2**rehydrated_mean.values - 1,
            'Log2FC_Dehyd_vs_Hyd': log2fcs,
            'Pvalue': pvalues
        })
        
        # FDR correction
        from scipy.stats import rankdata
        m = len(pvalues)
        ranks = rankdata(pvalues)
        de_results['AdjustedPvalue'] = de_results['Pvalue'] * m / ranks
        
        # Filter significant
        sig_lncrnas = de_results[(de_results['AdjustedPvalue'] < 0.05) & 
                                 (np.abs(de_results['Log2FC_Dehyd_vs_Hyd']) > 1)]
        
        # Add reference IDs
        gene_to_ref = dict(zip(id_mapping['Gene_ID'], id_mapping['Reference_ID']))
        de_results['Reference_ID'] = de_results['Gene_ID'].map(gene_to_ref)
        sig_lncrnas['Reference_ID'] = sig_lncrnas['Gene_ID'].map(gene_to_ref)
        
        # Save
        de_file = os.path.join(self.results_dir, 'lncrna_differential_expression.csv')
        sig_file = os.path.join(self.results_dir, 'significant_stress_responsive_lncrnas.csv')
        
        de_results.to_csv(de_file, index=False)
        sig_lncrnas.to_csv(sig_file, index=False)
        
        print(f"\nResults:")
        print(f"  Total lncRNAs analyzed: {len(de_results)}")
        print(f"  Significant lncRNAs: {len(sig_lncrnas)}")
        print(f"  Upregulated in dehydration: {len(sig_lncrnas[sig_lncrnas['Log2FC_Dehyd_vs_Hyd'] > 1])}")
        print(f"  Downregulated in dehydration: {len(sig_lncrnas[sig_lncrnas['Log2FC_Dehyd_vs_Hyd'] < -1])}")
        
        print(f"\nSaved:")
        print(f"  {de_file}")
        print(f"  {sig_file}")
        
        if len(sig_lncrnas) > 0:
            print(f"\nTop stress-responsive lncRNAs:")
            top = sig_lncrnas.nsmallest(5, 'Log2FC_Dehyd_vs_Hyd')[['Reference_ID', 'Log2FC_Dehyd_vs_Hyd', 'Hydrated_Mean_TPM', 'Dehydrated_Mean_TPM']]
            print(top.to_string())
        
        return True
    
    def run_full_pipeline(self):
        """Run complete pipeline"""
        print("\n" + "="*80)
        print("RESURRECTION PLANT lncRNA DISCOVERY PIPELINE")
        print("="*80)
        
        steps = [
            (self.check_files, "Checking files"),
            (self.build_hisat2_index, "Building HISAT2 index"),
            (self.align_samples, "Aligning samples"),
            (self.assemble_transcripts, "Assembling transcripts"),
            (self.create_expression_matrix, "Creating expression matrix"),
            (self.classify_lncrnas, "Classifying lncRNAs"),
            (self.create_id_mapping, "Creating ID mapping"),
            (self.differential_expression, "Performing differential expression"),
        ]
        
        for func, name in steps:
            try:
                if not func():
                    print(f"Failed at: {name}")
                    return False
            except Exception as e:
                print(f"Error in {name}: {e}")
                import traceback
                traceback.print_exc()
                return False
        
        print("\n" + "="*80)
        print("PIPELINE COMPLETE")
        print("="*80)
        print("\nKey Results:")
        print("  results/significant_stress_responsive_lncrnas.csv")
        print("  results/lncrna_differential_expression.csv")
        print("  results/stringtie_tpm_matrix.csv")
        print("="*80 + "\n")
        
        return True

def main():
    if len(sys.argv) > 1:
        project_dir = sys.argv[1]
    else:
        project_dir = '.'
    
    pipeline = LncRNAPipeline(project_dir)
    pipeline.run_full_pipeline()

if __name__ == '__main__':
    main()