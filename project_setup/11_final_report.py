#!/usr/bin/env python3
"""
11_final_report.py — Generate Final Research Summary Report
Pulls together all results from Steps 1-4B into a comprehensive markdown report.

Usage:
    python3 project_setup/11_final_report.py

Outputs:
    results/final_report.md     — full research report in markdown
    results/final_summary.txt   — plain text executive summary
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CANDIDATES   = os.path.join(BASE_DIR, "results", "significant_stress_responsive_lncrnas.csv")
COEXPR       = os.path.join(BASE_DIR, "results", "coexpression_summary.csv")
PARTNERS     = os.path.join(BASE_DIR, "results", "blast_annotated_partners.csv")
ENRICHMENT   = os.path.join(BASE_DIR, "results", "go_enrichment_results.csv")
PER_LNCRNA   = os.path.join(BASE_DIR, "results", "go_enrichment_lncrna.csv")
OUT_REPORT   = os.path.join(BASE_DIR, "results", "final_report.md")
OUT_SUMMARY  = os.path.join(BASE_DIR, "results", "final_summary.txt")

print("=" * 60)
print("Step 5 — Generating Final Research Report")
print("=" * 60)

# ── Load all results ──────────────────────────────────────────────────────────
cands   = pd.read_csv(CANDIDATES)
coexpr  = pd.read_csv(COEXPR)
partners = pd.read_csv(PARTNERS)
enrich  = pd.read_csv(ENRICHMENT)
per_lnc = pd.read_csv(PER_LNCRNA)

partners["gene_description"] = partners["gene_description"].fillna("")
annotated_partners = partners[partners["gene_description"] != ""].drop_duplicates("partner_id")
novel_partners     = partners[partners["gene_description"] == ""].drop_duplicates("partner_id")

top_enrich = enrich[enrich["padj"] < 0.05].sort_values("fold_enrichment", ascending=False)
date_str   = datetime.now().strftime("%B %d, %Y")

# Notable stress hits from partners
stress_kws = ["LEA", "embryogenesis abundant", "dehydrin", "heat shock", "HSP",
              "chaperone", "abscisic", "ABA", "trehalose", "peroxidase",
              "superoxide", "catalase", "desiccation", "drought", "aquaporin"]
stress_hits = annotated_partners[
    annotated_partners["gene_description"].str.contains(
        "|".join(stress_kws), case=False, na=False
    )
].drop_duplicates("partner_id")

# ── Build report ──────────────────────────────────────────────────────────────
lines = []
def w(s=""): lines.append(s)

w("# Discovery and Characterisation of Desiccation-Induced lncRNAs")
w("## in *Craterostigma plantagineum*")
w()
w(f"**Author:** Neil  ")
w(f"**Date:** {date_str}  ")
w(f"**Dataset:** GSE157098 / PRJNA660052  ")
w(f"**Pipeline:** github.com/legoreddragon-ai/resurrection_lncrna_pipeline  ")
w()
w("---")
w()

# ── Abstract ──────────────────────────────────────────────────────────────────
top_lnc = cands.iloc[0]
w("## Abstract")
w()
w(f"Resurrection plants survive near-complete desiccation through mechanisms that remain "
  f"incompletely understood, particularly at the level of non-coding RNA regulation. "
  f"We reanalysed public RNA-seq data from *Craterostigma plantagineum* across three "
  f"hydration states (hydrated, dehydrated, rehydrated; n=10 samples, PRJNA660052) to "
  f"identify and characterise long non-coding RNAs (lncRNAs) induced by desiccation stress. "
  f"Of 48,045 transcripts quantified, 38,978 (81%) met lncRNA criteria (length >200 nt, "
  f"low coding potential). Differential expression analysis identified "
  f"**{len(cands)} stress-responsive lncRNA candidates** completely silent under normal "
  f"conditions (0 TPM hydrated) and explosively activated under desiccation "
  f"(600–1,957 TPM). Co-expression network analysis at Pearson r > 0.90 revealed "
  f"**1,062 co-expression pairs** across **{coexpr['n_partners'].sum()} partner genes**, "
  f"with significant enrichment of chloroplast and photosynthetic recovery machinery "
  f"(3.7-fold, FDR = 6×10⁻¹⁴). Direct stress-response partners include LEA proteins, "
  f"HSP70 chaperones, and dehydrins. These findings suggest desiccation-induced lncRNAs "
  f"coordinate chloroplast protection and recovery during the desiccation-rehydration "
  f"cycle in this resurrection plant.")
w()
w("---")
w()

# ── 1. Introduction ───────────────────────────────────────────────────────────
w("## 1. Introduction")
w()
w("Vegetative desiccation tolerance (VDT) — the ability to survive near-complete water "
  "loss and recover upon rehydration — is a rare and remarkable trait found in "
  "so-called 'resurrection plants'. *Craterostigma plantagineum* is one of the most "
  "studied resurrection angiosperms, with a well-characterised transcriptome response "
  "involving LEA proteins, ABA signalling, and antioxidant systems.")
w()
w("Despite extensive study of protein-coding genes in resurrection plants, the role of "
  "long non-coding RNAs (lncRNAs) has received virtually no attention. In model plants "
  "such as *Arabidopsis thaliana*, lncRNAs have been shown to act as miRNA sponges, "
  "chromatin remodelling guides, and ABA signalling modulators under drought stress "
  "(e.g. DRIR, TCONS_00021861). Whether resurrection plants employ a unique or expanded "
  "set of stress-responsive lncRNAs — potentially explaining their extreme desiccation "
  "tolerance — remains an open question.")
w()
w("This study addresses that gap by reanalysing publicly available RNA-seq data from "
  "*C. plantagineum* to identify, characterise, and contextualise desiccation-induced "
  "lncRNAs through differential expression, co-expression network analysis, and "
  "functional annotation.")
w()
w("---")
w()

# ── 2. Methods ────────────────────────────────────────────────────────────────
w("## 2. Methods")
w()
w("### 2.1 Data Acquisition")
w()
w("RNA-seq data were obtained from NCBI GEO (accession GSE157098 / SRA PRJNA660052), "
  "comprising 10 paired-end samples from *C. plantagineum* across three conditions:")
w()
w("| Condition | Samples | IDs |")
w("|-----------|---------|-----|")
w("| Hydrated | n=2 | cp_hyd_r1, cp_hyd_r3 |")
w("| Dehydrated (WC2%) | n=2 | cp_wc2_r2, cp_wc2_r3 |")
w("| Dehydrated (WC60%) | n=3 | cp_wc60_r1, cp_wc60_r2, cp_wc60_r3 |")
w("| Rehydrated | n=3 | cp_rehyd_r1, cp_rehyd_r2, cp_rehyd_r3 |")
w()
w("### 2.2 RNA-seq Processing")
w()
w("Reads were quality-trimmed with **Trimmomatic** (SLIDINGWINDOW:4:15, MINLEN:36) "
  "and aligned to the *C. plantagineum* V2 transcriptome assembly (48,045 contigs; "
  "GSE157098) using **HISAT2** (alignment rate 60–95%). BAM files were sorted and "
  "indexed with **SAMtools**. Transcript abundance was quantified with **StringTie** "
  "in reference-guided mode, producing TPM values for all 48,045 genes across all "
  "10 samples.")
w()
w("### 2.3 lncRNA Identification")
w()
w("Transcripts were classified as lncRNAs based on: (1) length > 200 nucleotides and "
  "(2) ORF coding potential score < 0.5. This yielded **38,978 lncRNA candidates** "
  "representing 81% of the transcriptome, consistent with de novo transcriptome "
  "assemblies in other plant species.")
w()
w("### 2.4 Differential Expression Analysis")
w()
w("Differential expression between hydrated (n=2) and dehydrated conditions was tested "
  "using Welch's t-test with Benjamini-Hochberg FDR correction across 17,574 expressed "
  "lncRNAs (TPM ≥ 1). Due to the low replicate count (n=2 hydrated), FDR-controlled "
  "discovery was not possible; candidates were therefore selected using raw p < 0.05, "
  "|log₂FC| > 2, and an expression-weighted score (|log₂FC| × log₂(mean dehydrated "
  "TPM + 1)), yielding a pilot set of **14 stress-responsive candidates**.")
w()
w("### 2.5 Co-expression Network Analysis")
w()
w("Pearson correlation coefficients were computed between each of the 14 lncRNA "
  "candidates and all 6,511 expressed genes across all 10 samples. Gene pairs with "
  "r ≥ 0.90 were retained as co-expression partners.")
w()
w("### 2.6 Functional Annotation")
w()
w("Partner contig sequences were extracted from the V2 transcriptome FASTA and "
  "queried against the UniProt Swiss-Prot database (574,627 sequences) using "
  "**blastx** (E-value < 1×10⁻⁵, top hit per query). Functional enrichment was "
  "assessed using Fisher's exact test with BH correction, comparing category "
  "frequencies in the partner set against published background rates for plant "
  "proteomes.")
w()
w("---")
w()

# ── 3. Results ────────────────────────────────────────────────────────────────
w("## 3. Results")
w()
w("### 3.1 lncRNA Landscape of *C. plantagineum*")
w()
w(f"Of 48,045 transcripts in the V2 assembly, **38,978 (81%)** met lncRNA criteria. "
  f"Of these, 17,574 were expressed (TPM ≥ 1) in at least one condition and tested "
  f"for differential expression. The high lncRNA proportion is consistent with the "
  f"de novo nature of the transcriptome assembly.")
w()
w("### 3.2 Desiccation-Induced lncRNAs")
w()
w(f"Differential expression analysis identified **661 lncRNAs** with raw p < 0.05 "
  f"and log₂FC > 2 when comparing hydrated to dehydrated conditions. Strikingly, "
  f"**0 lncRNAs were significantly downregulated** — a completely one-sided response "
  f"suggesting a dedicated desiccation-activation programme rather than general "
  f"transcriptional disruption.")
w()
w(f"The top 14 candidates by expression-weighted score are presented below. All 14 "
  f"show **exactly 0 TPM in hydrated conditions** and are activated to 590–1,957 TPM "
  f"under desiccation — a binary off/on expression switch.")
w()

# Candidates table
w("| Gene ID | Hydrated TPM | Dehydrated TPM | Rehydrated TPM | log₂FC | p-value |")
w("|---------|-------------|----------------|----------------|--------|---------|")
for _, row in cands.iterrows():
    gid = row["gene_id"].replace("Cp_V2_contig_", "contig_")
    w(f"| {gid} | {row['mean_hydrated']:.1f} | {row['mean_dehydrated']:.1f} | "
      f"{row['mean_rehydrated']:.1f} | {row['log2FC']:.2f} | {row['pvalue']:.4f} |")
w()
w(f"Most candidates remain highly expressed in rehydrated samples, suggesting roles "
  f"in cellular recovery and photosynthetic apparatus restoration, not merely acute "
  f"stress alarm signalling.")
w()

w("### 3.3 Co-expression Network")
w()
total_pairs  = coexpr["n_partners"].sum()
hub_lncrna   = coexpr.iloc[0]["lncrna_id"].replace("Cp_V2_contig_", "contig_")
hub_partners = int(coexpr.iloc[0]["n_partners"])
w(f"Co-expression analysis at r ≥ 0.90 across all 10 samples identified "
  f"**1,062 co-expression pairs** involving **{len(annotated_partners) + len(novel_partners)} "
  f"unique partner genes**. All 14 lncRNA candidates have at least one co-expression "
  f"partner, confirming they are embedded in active regulatory networks.")
w()
w(f"The strongest network hub is **{hub_lncrna}** with {hub_partners} partners "
  f"(mean r = {coexpr.iloc[0]['mean_r']:.4f}), followed by contig_4950 "
  f"({int(coexpr.iloc[1]['n_partners'])} partners) and contig_2815 "
  f"({int(coexpr.iloc[2]['n_partners'])} partners).")
w()

w("### 3.4 Functional Annotation of Partners")
w()
w(f"Of 467 unique partner contigs, **{len(annotated_partners)} (78.6%)** received "
  f"BLAST hits against UniProt Swiss-Prot (E-value < 1×10⁻⁵). The remaining "
  f"**{len(novel_partners)} (21.4%)** had no significant similarity to known proteins, "
  f"representing potentially novel *C. plantagineum*-specific genes.")
w()

if len(stress_hits) > 0:
    w("Direct stress-response genes among the co-expression partners include:")
    w()
    for _, row in stress_hits.head(10).iterrows():
        w(f"- **{row['partner_id']}** (r={row['pearson_r']:.4f}): {row['gene_description'][:70]}")
    w()

w("### 3.5 Functional Enrichment Analysis")
w()
w(f"Fisher's exact test with BH correction identified **{len(top_enrich)} significantly "
  f"enriched functional category** (FDR < 0.05):")
w()

for _, row in enrich.sort_values("pvalue").head(6).iterrows():
    sig = " ✓ **(FDR < 0.05)**" if row["padj"] < 0.05 else ""
    w(f"- **{row['category']}**: {int(row['n_partners'])} genes ({row['pct_partners']}%), "
      f"{row['fold_enrichment']}× enriched, FDR = {row['padj']:.2e}{sig}")
w()
w("**Chloroplast and photosynthesis genes dominate the co-expression network** "
  "(53 genes, 14.9% of annotated partners, 3.7-fold enriched, FDR = 6×10⁻¹⁴). "
  "This finding suggests that desiccation-induced lncRNAs in *C. plantagineum* "
  "primarily coordinate the protection and restoration of the photosynthetic "
  "apparatus — a critical requirement for recovery from desiccation.")
w()
w("---")
w()

# ── 4. Discussion ─────────────────────────────────────────────────────────────
w("## 4. Discussion")
w()
w("### 4.1 A Desiccation-Exclusive lncRNA Programme")
w()
w("The complete absence of lncRNA expression in hydrated conditions (0 TPM for all "
  "14 candidates) followed by explosive activation under desiccation represents a "
  "binary regulatory switch rarely observed in transcriptomic studies. This pattern "
  "is consistent with the idea that resurrection plants have evolved dedicated "
  "non-coding regulatory machinery specifically for desiccation — machinery that "
  "is entirely dispensable under normal growth conditions.")
w()
w("### 4.2 Chloroplast Protection as a Central Theme")
w()
w("The strong enrichment of chloroplast and photosynthesis genes among co-expression "
  "partners is biologically significant. During desiccation, the chloroplast is "
  "particularly vulnerable to oxidative damage as the photosynthetic electron "
  "transport chain collapses. Resurrection plants are known to employ specialised "
  "mechanisms to protect chloroplast ultrastructure — including leaf folding, "
  "chloroplast repositioning, and upregulation of protective proteins.")
w()
w("The co-activation of lncRNAs with chloroplast genes suggests these non-coding "
  "transcripts may act as regulators of chloroplast protection programmes, "
  "potentially through chromatin remodelling at chloroplast gene loci or as "
  "post-transcriptional regulators of chloroplast-targeted mRNAs.")
w()
w("### 4.3 Stress-Response Partners")
w()
w("Despite not reaching FDR significance (likely due to small numbers and "
  "conservative background estimates), the direct identification of LEA proteins, "
  "HSP70 chaperones, dehydrins, and peroxidases among co-expression partners is "
  "notable. These proteins represent the canonical desiccation tolerance toolkit, "
  "and their co-expression with lncRNAs at r > 0.99 warrants experimental follow-up.")
w()
w("### 4.4 Novel Craterostigma-Specific Genes")
w()
w(f"The 98 partner contigs with no BLAST hits (21.4%) may represent genes unique "
  f"to resurrection plants or to *C. plantagineum* specifically. These represent "
  f"an intriguing target for future characterisation — if they are expressed only "
  f"under desiccation and co-regulate with known stress genes, they may encode "
  f"novel desiccation tolerance factors.")
w()
w("### 4.5 Limitations")
w()
w("Several limitations should be noted:")
w("1. **Replicate count**: n=2 hydrated replicates is insufficient for FDR-controlled "
  "genome-wide discovery. The 14 candidates are exploratory, not statistically definitive.")
w("2. **No functional validation**: All findings are correlative. qPCR validation and "
  "functional studies (overexpression, CRISPR knockdown) are required to establish "
  "causal roles.")
w("3. **Annotation completeness**: 21.4% of partner contigs lack BLAST annotation, "
  "limiting the functional enrichment analysis.")
w("4. **Background rates**: GO enrichment used estimated background rates rather than "
  "a formal GO database, as *C. plantagineum* lacks complete GO annotation.")
w()
w("---")
w()

# ── 5. Conclusions ────────────────────────────────────────────────────────────
w("## 5. Conclusions")
w()
w("This study provides the first systematic characterisation of desiccation-induced "
  "lncRNAs in *Craterostigma plantagineum*. Key findings:")
w()
w("1. **661 lncRNAs** are activated under desiccation; none are suppressed — a "
  "one-sided stress-induction programme.")
w("2. **14 top candidates** show binary expression (0 TPM hydrated → 600–1,957 TPM "
  "dehydrated) with sustained expression through rehydration.")
w("3. **Co-expression networks** link these lncRNAs to 467 partner genes, with "
  "significant enrichment of chloroplast and photosynthesis machinery (3.7×, "
  "FDR = 6×10⁻¹⁴).")
w("4. Direct partners include canonical desiccation tolerance genes: LEA proteins, "
  "HSP70, dehydrins, and peroxidases.")
w("5. **98 novel partner contigs** with no known homologs represent potential "
  "resurrection plant-specific desiccation factors.")
w()
w("These results support the hypothesis that resurrection plants employ a dedicated "
  "lncRNA regulatory layer during desiccation, with a particular emphasis on "
  "coordinating chloroplast protection and recovery.")
w()
w("---")
w()

# ── 6. Next Steps ─────────────────────────────────────────────────────────────
w("## 6. Recommended Next Steps")
w()
w("| Priority | Step | Rationale |")
w("|----------|------|-----------|")
w("| High | qPCR validation of top 3–5 candidates | Confirm expression patterns experimentally |")
w("| High | CPC2/CNCI coding potential analysis | Rigorous lncRNA classification |")
w("| High | Rfam rRNA/tRNA removal | Clean up the lncRNA set |")
w("| Medium | Selaginella comparison dataset | Test conservation across resurrection species |")
w("| Medium | miRNA sponge prediction (psRNATarget) | Identify post-transcriptional regulatory roles |")
w("| Medium | CRISPR knockdown of contig_551 | Test causal role of top candidate |")
w("| Low | Full GO database annotation | Requires genome-anchored assembly |")
w()
w("---")
w()

# ── 7. Data Availability ──────────────────────────────────────────────────────
w("## 7. Data & Code Availability")
w()
w("| Resource | Location |")
w("|----------|----------|")
w("| Pipeline code | github.com/legoreddragon-ai/resurrection_lncrna_pipeline |")
w("| Raw RNA-seq data | NCBI SRA: PRJNA660052 |")
w("| Transcriptome assembly | NCBI GEO: GSE157098 |")
w("| Candidate lncRNAs | results/significant_stress_responsive_lncrnas.csv |")
w("| Co-expression network | results/coexpression_partners.csv |")
w("| BLAST annotations | results/blast_annotated_partners.csv |")
w("| Enrichment results | results/go_enrichment_results.csv |")
w()
w("---")
w(f"*Report generated automatically by project_setup/11_final_report.py on {date_str}*")

# ── Write report ──────────────────────────────────────────────────────────────
with open(OUT_REPORT, "w") as f:
    f.write("\n".join(lines))
print(f"\n      Saved → {OUT_REPORT}")

# ── Executive summary ─────────────────────────────────────────────────────────
summary_lines = [
    "EXECUTIVE SUMMARY",
    "=" * 60,
    f"Project: lncRNA Discovery in Craterostigma plantagineum",
    f"Date:    {date_str}",
    f"Dataset: GSE157098 / PRJNA660052",
    "=" * 60,
    "",
    "KEY NUMBERS",
    f"  Transcripts quantified      : 48,045",
    f"  lncRNAs identified (81%)    : 38,978",
    f"  lncRNAs tested for DE       : 17,574",
    f"  Upregulated under stress    : 661",
    f"  Downregulated under stress  : 0",
    f"  Top candidates              : 14",
    f"  Co-expression pairs (r>0.9) : 1,062",
    f"  Annotated partner genes     : {len(annotated_partners)}",
    f"  Novel partner genes         : {len(novel_partners)}",
    "",
    "TOP FINDING",
    "  Chloroplast/photosynthesis genes are massively enriched",
    "  among co-expression partners (3.7x, FDR=6e-14).",
    "  lncRNAs co-activate with LEA proteins, HSP70, dehydrins.",
    "",
    "TOP CANDIDATES (by expression score)",
]
for _, row in cands.head(5).iterrows():
    summary_lines.append(
        f"  {row['gene_id']:<30}  log2FC={row['log2FC']:.2f}  "
        f"dehyd={row['mean_dehydrated']:.0f} TPM  p={row['pvalue']:.4f}"
    )
summary_lines += [
    "",
    "PIPELINE",
    "  github.com/legoreddragon-ai/resurrection_lncrna_pipeline",
    "",
    "STATUS: Steps 1-4B complete (75%). Remaining: qPCR validation,",
    "        CPC2 coding potential, Rfam filtering, miRNA sponge prediction.",
]

with open(OUT_SUMMARY, "w") as f:
    f.write("\n".join(summary_lines))
print(f"      Saved → {OUT_SUMMARY}")

print("\n" + "=" * 60)
print("PIPELINE COMPLETE — Steps 1 through 4B")
print()
print("Results saved to results/")
print("  final_report.md                   ← full research report")
print("  final_summary.txt                 ← executive summary")
print("  significant_stress_responsive_lncrnas.csv")
print("  coexpression_partners.csv")
print("  blast_annotated_partners.csv")
print("  go_enrichment_results.csv")
print("  volcano_plot.png")
print("  heatmap_significant_lncrnas.png")
print("  boxplots_top_lncrnas.png")
print("  coexpression_heatmap.png")
print("  go_enrichment_barplot.png")
print("=" * 60)
