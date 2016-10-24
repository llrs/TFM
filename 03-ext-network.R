# Sarting ####
source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

# Options to get the files
singFolder <- "../../RD/bicor"
consFolder <- "../../comparison_HA/bicor"

# input single expr ####
# Or any other expression one want to check if it holds.
load(file.path(singFolder, "modules_ME.RData"), verbose = TRUE)
# Rename variables to avoid conflicts
singLabels <- moduleColors
singColors <- moduleColors
singMEs <- orderMEs(MEs)
load(file.path(singFolder, "Input.RData"), verbose = TRUE)
singNames <- colnames(data.wgcna)
#Load the Consensus data # Consensus-modules_MEs.RData
load(file = file.path(consFolder, "modules_ME.RData"), verbose = TRUE)
load(file = file.path(consFolder, "Input.RData"), verbose = TRUE)
consNames <- colnames(data.wgcna) #[[1]]$data
consMEs <- MEs
if (length(singNames) > length(consNames)) {

  singNames <- select(org.Hs.eg.db, keys = singNames,
                      keytype = "REFSEQ", columns = "SYMBOL")[, "SYMBOL"]
  names(singLabels) <- singNames
  names(singColors) <- singNames
  comNames <- intersect(singNames, consNames)
  keep <- singNames %in% comNames
  singLabels <- singLabels[keep]
  singLabels <- singLabels[!duplicated(names(singLabels))]
  singColors <- singColors[keep]
  singColors <- singColors[!duplicated(names(singColors))]
}
if (length(singLabels) < length(moduleColors)){
  moduleColors <- moduleColors[consNames %in% comNames]
}
# This code assumes that the identifiers are in both datasets common
# and that are of the same size between them

# comparing modules ####
# Isolate the module labels in the order they appear in
# ordered module eigengenes
singModules <- substring(colnames(singMEs), 3)
singModules <- singModules[singModules %in% unique(singLabels)]
# Compares just against the first set of data not the consensus¿?
consModules <- substring(colnames(consMEs), 3) # [[1]]$data

# Numbers of single and consensus modules
nSingMods <- length(singModules)
nConsMods <- length(consModules)

# Initialize tables of p-values and of the corresponding counts
pTable <- matrix(0, nrow = nSingMods, ncol = nConsMods)
CountTbl <- matrix(0, nrow = nSingMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (smod in 1:nSingMods) {
  for (cmod in 1:nConsMods) {
    singMembers <- singLabels == singModules[smod]
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
pdf(file = "heatmap_RD_Vs_HA_Modules.pdf", width = 10, height = 7)
par(mfrow = c(1,1), cex = 1.0, mar = c(8, 10.4, 2.7, 1) + 0.3)
# Use function labeledHeatmap to produce the color-coded table with all the
# trimmings
colo <- log10(CountTbl)
colo[] <- ifelse(is.infinite(colo), 0, colo)
labeledHeatmap.multiPage(Matrix = colo,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", singModules),
               colorLabels = TRUE,
               xSymbols = paste0("HA ", consModules, ": ", consModTotals),
               ySymbols = paste0("RD ", singModules, ": ", singModTotals),
               textMatrix = CountTbl,
               signed = FALSE,
               # zlim = c(0, max(CountTbl)),
               addPageNumberToMain = FALSE,
               main = "Correspondence of modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE,
               maxRowsPerPage = 20, maxColsPerPage = 20)
dev.off()
