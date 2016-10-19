# Sarting ####
source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

# input single expr ####
# Or any other expression one want to check if it holds.
load("../..../modules_MEs.RData", verbose = TRUE)
# Rename variables to avoid conflicts
singLabels <- moduleColors
singColors <- moduleColors
singTree <- geneTree # The tree
singMEs <- orderMEs(MEs)

#Load the Consensus data
load(file = "Consensus-module_MEs.RData", verbose = TRUE)
consLabels <- moduleColors
consMEs <- MEs

# This code assumes that the identifiers are in both datasets common
# and that are of the same size between them

# comparing modules ####
# Isolate the module labels in the order they appear in
# ordered module eigengenes
singModules <- substring(colnames(singMEs), 3)
# Compares just against the first set of data not the consensus¿?
consModules <- substring(colnames(consMEs[[1]]$data), 3)

# Numbers of single and consensus modules
nSingMods <- length(singModules)
nConsMods <- length(consModules)

# Initialize tables of p-values and of the corresponding counts
pTable <- matrix(0, nrow = nSingMods, ncol = nConsMods)
CountTbl <- matrix(0, nrow = nSingMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (smod in 1:nSingMods) {
  for (cmod in 1:nConsMods) {
    singMembers <- singLabels == singModules[fsod]
    consMembers <- moduleColors == consModules[cmod]
    pTable[smod, cmod] <- -log10(fisher.test(singMembers, consMembers,
                                             alternative = "greater")$p.value)
    CountTbl[smod, cmod] <- sum(singColors == singModules[smod] &
                                  moduleColors == consModules[cmod])
  }
}

# heatmap ####

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] <- 1.3*max(pTable[is.finite(pTable)])
pTable[pTable > 50 ] <- 50
# Marginal counts (really module sizes)
singModTotals <- apply(CountTbl, 1, sum)
consModTotals <- apply(CountTbl, 2, sum)

# Actual plotting
pdf(file = "heatmap_ConsensusVsSingleModules.pdf", width = 10, height = 7)
par(mfrow = c(1,1), cex = 1.0, mar = c(8, 10.4, 2.7, 1) + 0.3)
# Use function labeledHeatmap to produce the color-coded table with all the
# trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", singModules),
               colorLabels = TRUE,
               xSymbols = paste0("Cons ", consModules, ": ", consModTotals),
               ySymbols = paste0("Sing ", singModules, ": ", singModTotals),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of set-specific and consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
dev.off()
