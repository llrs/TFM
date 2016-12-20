# Sarting ####
source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

# Options to get the files
singFolder <- "../../RD/signed_signed"
consFolder <- "../../comparison_HA/signed_signed"

# input single expr ####
# Or any other expression one want to check if it holds.
load(file.path(singFolder, "modules_ME.RData"), verbose = TRUE)
# Rename variables to avoid conflicts
singLabels <- moduleColors
singColors <- moduleColors
sing.modules <- moduleColors
singMEs <- orderMEs(MEs)
load(file.path(singFolder, "Input.RData"), verbose = TRUE)
singNames <- colnames(data.wgcna)
sing.expr <- data.wgcna
#Load the Consensus data # Consensus-modules_MEs.RData
load(file = file.path(consFolder, "modules_ME.RData"), verbose = TRUE)
load(file = file.path(consFolder, "Input.RData"), verbose = TRUE)
consNames <- colnames(data.wgcna) #[[1]]$data
cons.expr <- data.wgcna
cons.modules <- moduleColors
consMEs <- MEs
keytype <- "REFSEQ"
if (length(singNames) > length(consNames)) {
  # Convert all into SYMBOLS And only keep those
  names.genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
                                              keys = singNames,
                                              keytype = keytype,
                                              columns = "SYMBOL"))
  name <- data.frame(keytype = singNames, mod = singLabels)
  ids <- merge(name, names.genes, by.y = keytype, by.x = "keytype",
               sort = FALSE)
  singLabels <- ids$SYMBOL
  singColors <- ids$mod
  sing <- as.character(singColors)
  names(sing) <- singLabels
  cons <- moduleColors
  names(cons) <- consNames

  comNames <- intersect(names(sing), names(cons))
  keep.sing <- names(sing) %in% comNames
  keep.sing <- keep.sing & !duplicated(names(sing))
  sing <- sing[keep.sing]
  keep.cons <- names(cons) %in% comNames & !duplicated(names(cons))
  cons <- cons[keep.cons]
}
if (length(singLabels) < length(moduleColors)){
  moduleColors <- moduleColors[consNames %in% comNames]
}

# comparing modules ####
comparison <- table(sing, cons)
# comparison <- matrix(as.data.frame.matrix(comparison))
# Isolate the module labels in the order they appear in
# ordered module eigengenes
singModules <- unique(sing)
# Compares just against the first set of data not the consensus¿?
consModules <- unique(cons) # [[1]]$data

# Numbers of single and consensus modules
nSingMods <- length(singModules)
nConsMods <- length(consModules)

# Initialize tables of p-values and of the corresponding counts
pTable <- matrix(0, nrow = nSingMods, ncol = nConsMods,
                 dimnames = list(rownames = singModules,
                                 colnames = consModules))
# CountTbl <- matrix(0, nrow = nSingMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (smod in 1:nSingMods) {
  for (cmod in 1:nConsMods) {
    singMembers <- sing == singModules[smod]
    consMembers <- cons == consModules[cmod]
    pTable[smod, cmod] <- -log10(fisher.test(singMembers, consMembers,
                                             alternative = "greater")$p.value)
    # CountTbl[smod, cmod] <- sum(sing == singModules[smod] &
    #                               cons == consModules[cmod])
  }
}

# heatmap ####

# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] <- 1.3*max(pTable[is.finite(pTable)])
pTable[pTable > 50 ] <- 50
# Marginal counts (really module sizes)
CountTbl <- comparison
singModTotals <- apply(CountTbl, 1, sum)
consModTotals <- apply(CountTbl, 2, sum)

pTable <- pTable[match(rownames(comparison), rownames(pTable)),
                 match(colnames(comparison), colnames(pTable))]
# Actual plotting
pdf(file = "heatmap_RD_Vs_HA_Modules.pdf", width = 10, height = 7)
par(mfrow = c(1,1), cex = 1.0, mar = c(8, 10.4, 2.7, 1) + 0.3)
# Use function labeledHeatmap to produce the color-coded table with all the
# trimmings
# colo <- log10(CountTbl)
# colo[] <- ifelse(is.infinite(colo), 0, colo)
labeledHeatmap.multiPage(Matrix = pTable, #-log10 p-value
               xLabels = paste(" ", colnames(CountTbl)),
               yLabels = paste(" ", rownames(CountTbl)),
               colorLabels = TRUE,
               xSymbols = paste0("HA ", colnames(CountTbl), ": ", consModTotals),
               ySymbols = paste0("RD ", rownames(CountTbl), ": ", singModTotals),
               textMatrix = CountTbl,
               signed = FALSE,
               # zlim = c(0, max(CountTbl)),
               addPageNumberToMain = FALSE,
               main = "Correspondence of modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE,
               maxRowsPerPage = 20, maxColsPerPage = 20)
dev.off()

# Modules preservation ####
setLabels <- c("RD", "HA")
sing.expr <- sing.expr[, keep.sing]
colnames(sing.expr) <- names(sing)
cons.expr <- cons.expr[, keep.cons]
multiExpr <- list(RD = list(data = sing.expr), HA = list(data = cons.expr))
multiColor <- list(RD = sing, HA = cons)

# Calculate the preservation
# mp <- modulePreservation(multiExpr, multiColor,
#                         referenceNetworks = 2,
#                         nPermutations = 200,
#                         networkType = "signed",
#                         corFnc = "bicor",
#                         randomSeed = 1,
#                         parallelCalculation = TRUE,
#                         verbose = 3)
# save(mp, file = "modulePreservation.RData")
load(file = "modulePreservation.RData", verbose = TRUE)

# Select the output
ref <- 1 # Adjusted to work
test <- 1
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                  mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                mp$preservation$Z[[ref]][[test]][, -1])
print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)))

# Module labels and module sizes are also contained in the results
modColors <- rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1]

# leave grey and gold modules out
plotMods <- !(modColors %in% c("grey", "gold"))
# Text labels for points
text <- modColors[plotMods]
# Auxiliary convenience variable
plotData <- cbind(mp$preservation$observed[[ref]][[test]][, 2],
                  mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains <- c("Preservation Median rank", "Preservation Zsummary");
pdf("Preservation_plots.pdf")
# Start the plot
for (p in 1:2) {
  min <- min(plotData[, p], na.rm = TRUE)
  max <- max(plotData[, p], na.rm = TRUE)
  # Adjust ploting ranges appropriately
  if (p == 2) {
    if (min > -max/10) {
      min <- -max/10
    }
    ylim <- c(min - 0.1 * (max - min), max + 0.1 * (max - min))
  } else{
    ylim <- c(max + 0.1 * (max - min), min - 0.1 * (max - min))
  }
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1,
       bg = modColors[plotMods], pch = 21, main = mains[p],
       cex = 2.4, ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim, xlim = c(10, 2000), cex.lab = 1.2,
       cex.axis = 1.2, cex.main = 1.4)
  # labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1,
  #             offs = 0.08)
  # For Zsummary, add threshold lines
  if (p == 2) {
    abline(h = 0)
    abline(h = 2, col = "blue", lty = 2)
    abline(h = 10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off()
