# ==============================================================================
#
#  Code chunk 1: Starting from the preiously saved data
#
# ==============================================================================

pdfn <- function(...){
  # Close any device and open a pdfn with the same options
  if (length(dev.list()) > 1) {
    dev.off()
  }
  pdf(...)
}

library("WGCNA")
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
enableWGCNAThreads(6)
# Load the data saved in the first part
load(file = "InputWGCNA.RData", verbose = TRUE)
#The variable lnames contains the names of loaded variables.
library("biomaRt")
library("hgu133plus2.db")
library("GOstats")
library("graphite")
library("KEGGgraph")
library("KEGG.db")
library("RBGL")
source("bio_cor.R")

# ==============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================

ncol(data.wgcna)
bio_mat <- bio.cor(colnames(data.wgcna))
save(bio_mat, file = "bio_correlation.RData")
# Choose a set of soft-thresholding powers
powers = c(1:30)
# Call the network topology analysis function
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, verbose = 5,
                        networkType = "signed", corFnc = cor.all,
                        corOptions = list(use = "p", bio_mat = bio_mat))
# Plot the results:
pdfn(file = "Power_calculations_bio.cor.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,
     col = "red")

stop("Testing bio.all function")
# ==============================================================================
#
#  Code chunk 3: Automatic blocks creation using the power calculated
#
# ==============================================================================

# print(paste("Recomended power", sft$powerEstimate))
# net <- blockwiseModules(data.wgcna, power = 12, # the max connectivity of 0.77
#                 # Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
#                 # 12     0.7730  -1.490          0.738     602     400.0   2220
#                        TOMType = "signed", minModuleSize = 30,
#                        maxBlockSize = 8000, networkType = "signed",
#                        pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "AH_signed",
#                        verbose = 3)

# save(net, file = "net.RData")
load("net.RData")
# ==============================================================================
#
#  Code chunk 4: Plot how the modules correlate in a first dendrogram
#
# ==============================================================================


# open a graphics window
pdfn(file = "dendro.pdf", width = 12, height = 9)
# Convert labels to colors for plotting
mergedColors <- net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors of first dendrogram",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# ==============================================================================
#
#  Code chunk 5: See how are the modules found
#
# ==============================================================================

# Calculate eigengenes
MEList <- moduleEigengenes(data.wgcna, colors = net$colors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
corME <- cor(MEs)
MEDiss <- 1 - corME
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
pdfn("Modules_relationship.pdf")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

labeledHeatmap(corME, xLabels = colnames(corME), yLabels = colnames(corME),
               sSymbols = colnames(corME), ySymbols = colnames(corME))

gm <- table(net$colors)
gm
hist(gm)
count.p <- function(data, per){
  # Calculate the percentatge above per of data
  sum(data >= per)/length(data)
}
perc <- unlist(lapply(gm, count.p, data = gm))
plot(cbind(gm[order(perc)],perc[order(perc)]), type = "o",
     xlab = "Size of the modules",
     ylab = "Proportion of modules above the size",
     main = "Distribution of the size of the modules",
     col = names(gm[order(perc)]))





# Call an automatic merging function
merge <- mergeCloseModules(data.wgcna,
                          net$colors, cutHeight = MEDissThres, MEs = MEs,
                          verbose = 3)
# The merged module colors
mergedColors <- merge$colors;
# Eigengenes of the new merged modules:
mergedMEs <- merge$newME

plot(merge$dendro, main = "After merging the modules")

dev.off()
# ==============================================================================
#
#  Code chunk 6: Save the data
#
# ==============================================================================

moduleColors <- mergedColors
MEs <- mergedMEs
save(MEs, moduleColors, file = "TNF_AH-network-auto.RData")
