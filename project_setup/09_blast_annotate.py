#!/usr/bin/env python3
"""
09_blast_annotate.py — BLAST co-expression partners against UniProt Swiss-Prot
Extracts the 467 unique partner contig sequences and runs blastx to get gene names/functions.

Usage:
    python3 project_setup/09_blast_annotate.py

Outputs:
    results/partner_sequences.fa        — FASTA of 467 partner contigs
    results/blast_results.txt           — raw blastx output
    results/blast_annotated_partners.csv — partners with gene names and functions
"""

import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTNERS    = os.path.join(BASE_DIR, "results", "coexpression_partners.csv")
TRANSCRIPTS = os.path.join(BASE_DIR, "references", "GSE157098_Cp_transcriptome_assembly_V2.fa")
BLAST_DB    = os.path.join(BASE_DIR, "references", "blast_db", "uniprot_sprot")
OUT_FASTA   = os.path.join(BASE_DIR, "results", "partner_sequences.fa")
OUT_BLAST   = os.path.join(BASE_DIR, "results", "blast_results.txt")
OUT_ANNOT   = os.path.join(BASE_DIR, "results", "blast_annotated_partners.csv")

THREADS     = 4      # increase if you have more cores
E_VALUE     = "1e-5" # standard threshold
MAX_TARGETS = 1      # top hit only

print("=" * 60)
print("Step 4B/5A — BLAST Annotation of Co-expression Partners")
print("=" * 60)

# ── Get unique partner IDs ────────────────────────────────────────────────────
print("\n[1/4] Loading partner gene IDs...")
pairs = pd.read_csv(PARTNERS)
partner_ids = set(pairs["partner_id"].unique())
print(f"      {len(partner_ids)} unique partner contigs to annotate")

# ── Extract sequences from transcriptome FASTA ────────────────────────────────
print("[2/4] Extracting sequences from transcriptome FASTA...")
found = 0
with open(OUT_FASTA, "w") as out_f:
    for record in SeqIO.parse(TRANSCRIPTS, "fasta"):
        # contig IDs in FASTA may be like "Cp_V2_contig_10455" or just "contig_10455"
        rec_id = record.id
        # try direct match first, then with prefix
        if rec_id in partner_ids:
            SeqIO.write(record, out_f, "fasta")
            found += 1
        elif "Cp_V2_" + rec_id in partner_ids:
            record.id = "Cp_V2_" + rec_id
            record.description = ""
            SeqIO.write(record, out_f, "fasta")
            found += 1

print(f"      {found} sequences extracted → {OUT_FASTA}")
if found == 0:
    print("  ⚠  No sequences found. Checking FASTA ID format...")
    # Print first few IDs from FASTA to debug
    for i, rec in enumerate(SeqIO.parse(TRANSCRIPTS, "fasta")):
        print(f"      FASTA ID example: '{rec.id}'")
        if i >= 4:
            break
    print("      Partner ID example:", list(partner_ids)[:3])
    sys.exit(1)

# ── Run blastx ───────────────────────────────────────────────────────────────
print(f"[3/4] Running blastx against UniProt Swiss-Prot...")
print(f"      {found} queries × 574,627 proteins — this may take 10-30 minutes...")

blast_cmd = [
    "blastx",
    "-query",    OUT_FASTA,
    "-db",       BLAST_DB,
    "-out",      OUT_BLAST,
    "-outfmt",   "6 qseqid sseqid pident length evalue bitscore stitle",
    "-evalue",   E_VALUE,
    "-max_target_seqs", str(MAX_TARGETS),
    "-num_threads", str(THREADS),
]

print(f"      Command: {' '.join(blast_cmd)}")
result = subprocess.run(blast_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  ✗ blastx failed:\n{result.stderr}")
    sys.exit(1)
print(f"      Done → {OUT_BLAST}")

# ── Parse BLAST results ───────────────────────────────────────────────────────
print("[4/4] Parsing BLAST results and annotating partners...")

blast_cols = ["query_id", "subject_id", "pident", "length", "evalue", "bitscore", "stitle"]
try:
    blast = pd.read_csv(OUT_BLAST, sep="\t", header=None, names=blast_cols)
    # Keep best hit per query (already max_target_seqs=1 but just in case)
    blast = blast.sort_values("bitscore", ascending=False).drop_duplicates("query_id")
    print(f"      {len(blast)} contigs got BLAST hits")
    print(f"      {len(partner_ids) - len(blast)} contigs had no hit (novel/plant-specific)")
except Exception:
    blast = pd.DataFrame(columns=blast_cols)
    print("      No BLAST hits found at e-value < 1e-5")

# Parse gene name and organism from stitle
# Swiss-Prot format: "sp|Q9FVA1|LEA_ARATH Late embryogenesis abundant protein ... OS=Arabidopsis thaliana"
def parse_gene_name(stitle):
    if pd.isna(stitle):
        return "unknown", "unknown", "unknown"
    parts = stitle.split(" ", 2)
    gene_desc = parts[2] if len(parts) > 2 else stitle
    # Extract OS (organism)
    os_part = ""
    if "OS=" in gene_desc:
        os_part = gene_desc.split("OS=")[1].split(" OX=")[0].strip()
        gene_desc = gene_desc.split("OS=")[0].strip()
    # Extract GN (gene name)
    gn_part = ""
    if "GN=" in stitle:
        gn_part = stitle.split("GN=")[1].split(" ")[0].strip()
    return gene_desc[:80], gn_part, os_part

blast["gene_description"], blast["gene_name"], blast["organism"] = zip(
    *blast["stitle"].apply(parse_gene_name)
)

# Merge with pairs
annotated = pairs.merge(
    blast[["query_id", "gene_description", "gene_name", "organism", "pident", "evalue", "bitscore"]],
    left_on="partner_id", right_on="query_id", how="left"
).drop(columns=["query_id"])

annotated.to_csv(OUT_ANNOT, index=False)
print(f"      Saved → {OUT_ANNOT}")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n── Top Co-expressed Partners with Annotations ─────────────")
top = (annotated.dropna(subset=["gene_description"])
       .sort_values("pearson_r", ascending=False)
       .drop_duplicates("partner_id")
       .head(20))

for _, row in top.iterrows():
    print(f"  r={row['pearson_r']:.4f}  {row['partner_id']:<25}  {row['gene_description'][:55]}")

# Stress-relevant keyword scan
print("\n── Stress-Pathway Hits ─────────────────────────────────────")
keywords = ["LEA", "dehydrin", "ABA", "abscisic", "trehalose", "antioxidant",
            "superoxide", "catalase", "peroxidase", "heat shock", "HSP",
            "late embryogenesis", "desiccation", "drought", "osmotic",
            "aquaporin", "chaperone", "ubiquitin", "proteasome"]

for kw in keywords:
    hits = annotated[annotated["gene_description"].str.contains(kw, case=False, na=False)]
    if len(hits) > 0:
        print(f"  {kw:<20} {len(hits):>4} pairs  —  e.g. {hits.iloc[0]['gene_description'][:45]}")

print("\n" + "=" * 60)
print("Step 4B/5A Complete")
print(f"  Total pairs annotated : {len(annotated[annotated['gene_description'].notna()]):,}")
print(f"  Novel (no BLAST hit)  : {len(annotated[annotated['gene_description'].isna()]):,}")
print(f"  Next step: python3 project_setup/10_go_enrichment.py")
print("=" * 60)
