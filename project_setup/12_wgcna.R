library(WGCNA)
library(impute)
library(preprocessCore)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ── Load TPM matrix ───────────────────────────────────────────────────────────
tpm <- read.csv("results/stringtie_tpm_matrix.csv", row.names = 1)
tpm <- tpm[, -1]  # drop Gene Name column

# Sample names and condition labels
colnames(tpm)
conditions <- c("hydrated","hydrated",
                "rehydrated","rehydrated","rehydrated",
                "desiccated","desiccated","desiccated","desiccated","desiccated")

# ── Filter: mean TPM >= 1 across all samples ──────────────────────────────────
tpm_filt <- tpm[rowMeans(tpm) >= 1, ]
cat("Genes after filtering:", nrow(tpm_filt), "\n")

# ── Transpose for WGCNA (samples x genes) ────────────────────────────────────
datExpr <- as.data.frame(t(log2(tpm_filt + 1)))

# Check for bad genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
cat("Final dimensions:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# ── Pick soft threshold ───────────────────────────────────────────────────────
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf("results/figures/wgcna_soft_threshold.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology R^2",
     main="Scale independence")
abline(h=0.85, col="red")
dev.off()

# Pick lowest power where R^2 >= 0.85
chosen_power <- sft$fitIndices[which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] >= 0.85)[1], 1]
cat("Chosen soft threshold power:", chosen_power, "\n")
if (is.na(chosen_power)) {
  chosen_power <- 12
  cat("R^2 threshold not reached, defaulting to power =", chosen_power, "\n")
}

# ── Build network and detect modules ─────────────────────────────────────────
net <- blockwiseModules(datExpr,
                        power = chosen_power,
                        minModuleSize = 30,
                        mergeCutHeight = 0.25,
                        numericLabels = FALSE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("Modules found:", table(net$colors), "\n")

# ── Module-trait correlation ──────────────────────────────────────────────────
traitData <- data.frame(
  hydrated    = as.numeric(conditions == "hydrated"),
  desiccated  = as.numeric(conditions == "desiccated"),
  rehydrated  = as.numeric(conditions == "rehydrated"),
  row.names   = rownames(datExpr)
)

MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Save correlation table
write.csv(as.data.frame(moduleTraitCor), "results/wgcna_module_trait_correlation.csv")
cat("Module-trait correlations saved.\n")

# ── Plot module-trait heatmap ─────────────────────────────────────────────────
pdf("results/figures/wgcna_module_trait_heatmap.pdf")
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData),
               yLabels = rownames(moduleTraitCor),
               ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()

# ── Find which module contains our 14 candidates ─────────────────────────────
candidates <- c("Cp_V2_contig_551","Cp_V2_contig_12238","Cp_V2_contig_4950",
                "Cp_V2_contig_2815","Cp_V2_contig_8442","Cp_V2_contig_15059",
                "Cp_V2_contig_43263","Cp_V2_contig_7701","Cp_V2_contig_9183",
                "Cp_V2_contig_22145","Cp_V2_contig_31089","Cp_V2_contig_5534",
                "Cp_V2_contig_18802","Cp_V2_contig_6671")

candidate_modules <- net$colors[names(net$colors) %in% candidates]
cat("\nModule assignments for 14 candidates:\n")
print(candidate_modules)

write.csv(as.data.frame(candidate_modules),
          "results/wgcna_candidate_modules.csv")

# ── Save full module membership ───────────────────────────────────────────────
module_membership <- data.frame(gene = names(net$colors), module = net$colors)
write.csv(module_membership, "results/wgcna_module_membership.csv", row.names = FALSE)
cat("All done. Results in results/wgcna_*.csv and results/figures/\n")
