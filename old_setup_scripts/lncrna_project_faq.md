# FAQ: Discovery of Stress-Responsive lncRNAs in Resurrection Plants
### A Self-Study Guide — From Simple to Complex
*For an 11th-grade researcher working on Craterostigma plantagineum*

---

## SECTION 1: The Big Picture (Start Here)

**Q1. What is a resurrection plant, and why is it scientifically interesting?**

**Q2. What is RNA, and how is it different from DNA?**

**Q3. What is a gene, and what does "gene expression" mean?**

**Q4. If most genes make proteins, what does a non-coding RNA actually do?**

**Q5. What is a long non-coding RNA (lncRNA), and how is it defined?**

**Q6. Why would a plant produce RNA that doesn't make a protein — isn't that wasteful?**

**Q7. What is desiccation, and what happens to a plant's cells when it dries out?**

**Q8. What is the central dogma of molecular biology, and where do lncRNAs fit in?**

---

## SECTION 2: Understanding the Experiment

**Q9. What is RNA sequencing (RNA-seq), and what does it measure?**

**Q10. Why do we use paired-end sequencing, and what does that mean practically?**

**Q11. What is a biological replicate, and why do we need more than one?**

**Q12. What are the three conditions in this experiment (hydrated, dehydrated, rehydrated), and what biological question does each help answer?**

**Q13. What is a FASTQ file, and what information does it contain?**

**Q14. Why do we need to trim reads before analysis? What are adapters?**

**Q15. What is a reference transcriptome, and why do we align reads to it?**

**Q16. What does HISAT2 do, and why is alignment necessary?**

**Q17. What is a BAM file, and how is it different from a FASTQ file?**

**Q18. What does StringTie do, and what is transcript assembly?**

---

## SECTION 3: Understanding the Numbers

**Q19. What is TPM (Transcripts Per Million), and why do we use it instead of raw read counts?**

**Q20. Why does a gene with 0 TPM in dehydrated samples but 1,320 TPM in hydrated samples matter biologically?**

**Q21. What is fold change, and what does log₂ fold change mean?**

**Q22. What is a p-value, and why isn't a small p-value alone enough to call something significant?**

**Q23. What is the Benjamini-Hochberg correction, and what problem does it solve?**

**Q24. What is an adjusted p-value (padj), and what does padj < 0.05 mean in practice?**

**Q25. Why did we use the thresholds padj < 0.05 AND |log2FC| > 1 together?**

**Q26. With only 2 hydrated replicates, how reliable is our t-test statistically?**

---

## SECTION 4: The lncRNA Identification Process

**Q27. How do we distinguish a lncRNA from a regular protein-coding gene computationally?**

**Q28. What is an ORF (Open Reading Frame), and why does ORF length matter for classification?**

**Q29. What is coding potential, and how is it scored?**

**Q30. In our pipeline, 81% of transcripts were classified as lncRNA — does that seem too high, and why might that be?**

**Q31. What is differential expression analysis, and what question does it answer?**

**Q32. Why did we focus on lncRNAs downregulated during dehydration specifically?**

**Q33. What does it mean biologically if an lncRNA is highly expressed only when the plant is hydrated?**

---

## SECTION 5: Interpreting the Results

**Q34. What is Cp_V2_contig_159 (STRG.167), and why is it the top candidate?**

**Q35. What does a "hydration-maintenance" regulatory role mean, and how would you test that hypothesis?**

**Q36. What is a volcano plot, and how do you read one?**

**Q37. What is a heatmap in gene expression, and what does the z-score normalization show?**

**Q38. Why do we look at rehydrated samples — what can recovery tell us that dehydration alone cannot?**

**Q39. Could any of the 14 lncRNAs be false positives? How would you find out?**

---

## SECTION 6: Mechanisms — How lncRNAs Work

**Q40. What are the main known mechanisms by which lncRNAs regulate gene expression?**

**Q41. What is chromatin remodeling, and how could an lncRNA influence it?**

**Q42. What is a miRNA sponge (competing endogenous RNA), and how might an lncRNA act as one?**

**Q43. What is cis-regulation vs. trans-regulation, and which is more likely for our lncRNAs?**

**Q44. What are LEA proteins, and why are they important in desiccation tolerance?**

**Q45. What is ABA (abscisic acid) signaling, and how does it connect to stress response?**

**Q46. Could our lncRNAs be regulating ABA signaling? How would you design an experiment to test this?**

---

## SECTION 7: Co-expression and Network Analysis (Next Steps)

**Q47. What is a co-expression network, and what does it tell us?**

**Q48. What is Pearson correlation, and how is it used to connect lncRNAs to coding genes?**

**Q49. What is WGCNA, and why is it commonly used in transcriptomics?**

**Q50. If lncRNA X is highly co-expressed with a drought-stress gene Y, what can we conclude — and what can we NOT conclude?**

**Q51. What is guilt-by-association, and why is it a useful but limited inference strategy?**

---

## SECTION 8: Functional Enrichment

**Q52. What is Gene Ontology (GO), and how is it organized?**

**Q53. What does GO enrichment analysis tell us about a set of genes?**

**Q54. What is KEGG, and how is a KEGG pathway different from a GO term?**

**Q55. If our co-expressed coding genes are enriched for "response to water deprivation," what does that suggest about the lncRNAs?**

---

## SECTION 9: Validation and Experimental Biology

**Q56. What is qPCR, and how would you use it to validate RNA-seq results?**

**Q57. What is CRISPR, and how could you use it to study lncRNA function in plants?**

**Q58. What is RNA interference (RNAi), and is it an alternative to CRISPR for studying lncRNAs?**

**Q59. What would a loss-of-function experiment for an lncRNA look like, and what phenotype would you hope to see?**

**Q60. Why is wet-lab validation important even when computational evidence is strong?**

---

## SECTION 10: Novelty and Cross-Species Comparisons

**Q61. What is BLAST, and how do you use it to assess whether an lncRNA is novel?**

**Q62. Why do lncRNAs tend to be less conserved across species than protein-coding genes?**

**Q63. What other resurrection plants exist, and why would comparing lncRNAs across species be valuable?**

**Q64. What is Selaginella lepidophylla, and how does its desiccation tolerance compare to Craterostigma?**

**Q65. Could any of our 14 lncRNAs have homologs in Arabidopsis thaliana — and why would that matter?**

---

## SECTION 11: Broader Applications

**Q66. How could understanding desiccation-tolerance lncRNAs help in crop science?**

**Q67. What is drought stress in agriculture, and why is it one of the biggest threats to food security?**

**Q68. Could lncRNAs from resurrection plants be transferred into crops to improve drought tolerance?**

**Q69. How is lncRNA research being applied in human medicine (cancer, neurological diseases)?**

**Q70. What is the connection between plant stress biology and climate change research?**

---

## SECTION 12: Thinking Like a Scientist

**Q71. What are the biggest limitations of this study, and how would you address them?**

**Q72. What is the difference between correlation and causation in biology?**

**Q73. How do you decide which of the 14 lncRNAs to prioritize for follow-up study?**

**Q74. What makes a good scientific hypothesis, and how would you frame one for this project?**

**Q75. If you were to publish this research, what journal would be appropriate and why?**

**Q76. What is peer review, and why is it important for scientific credibility?**

**Q77. How would you present this project at a science fair or research symposium?**

**Q78. What computational skills should you learn next to advance this kind of research?**

**Q79. What biology courses or university programs build on this kind of work?**

**Q80. What does a career in bioinformatics, plant biology, or computational genomics look like?**

---

*Work through these questions section by section. Sections 1–5 cover what you've already done. Sections 6–10 cover where the project is going. Sections 11–12 will help you think and communicate like a real researcher.*
