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

library(WGCNA)
options(stringsAsFactors = FALSE) 
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
enableWGCNAThreads(6)
# Load the data saved in the first part
load(file = "Input.RData", verbose = TRUE) 
#The variable lnames contains the names of loaded variables.


data.wgcna <- t(exprs[, pheno$group  == "AH"])

# ==============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================
message("Calculating the power to use")

# Choose a set of soft-thresholding powers
powers <- c(1:30)
# Call the network topology analysis function
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, verbose = 5,
                        networkType = "signed")
# Plot the results:
pdfn(file = "Power_calculations.pdf", width = 9, height = 5)
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

# ==============================================================================
#
#  Code chunk 3: Calculating adjacency
#
# ============================================================================== 
message(paste("Calculating adjacency with power", sft$powerEstimate))

if (1/sqrt(ncol(data.wgcna)) ^ sft$powerEstimate * ncol(data.wgcna) >= 0.1) {
  warning("Are you sure of this power?")
}
softPower <- sft$powerEstimate
pdf("scalefree.pdf")
k <- softConnectivity(data.wgcna, type = "signed", power = softPower)
hist(k, main = paste("Connectivity for the given power of", softPower))
scaleFreePlot(k, main = "Check scale free topology\n")
dev.off()
adjacency  <- adjacency(data.wgcna, power = softPower, type = "signed")


# ==============================================================================
#
#  Code chunk 4: Calculating TOM and TOMdiss
#
# ==============================================================================
message("Calculating TOM and TOMdiss")

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# ==============================================================================
#
#  Code chunk 5: Building the hierarchical clustering
#
# ==============================================================================
message("Building the hierarchical clustering")

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
png("TOMdiss_clustring.png")
plot(geneTree, xlab = "", sub = "", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# ==============================================================================
#
#  Code chunk 6: Merging modules
#
# ==============================================================================
message("Merging modules")

# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 5
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)


# ==============================================================================
#
#  Code chunk 7: Dendrogram of genes
#
# ==============================================================================
message("Building a dendrogram of genes with merged modules")

# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
png("dendrogram.png")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# ==============================================================================
#
#  Code chunk 8: Calculate the eigengenes
#
# ==============================================================================
message("Calculate the eigengenes")

# Calculate eigengenes
MEList <- moduleEigengenes(data.wgcna, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
png("EV_modules.png")
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

# ==============================================================================
#
#  Code chunk 9: Merge with fixed height
#
# ==============================================================================
message("Fixed merging")

MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
# Call an automatic merging function
merge <- mergeCloseModules(data.wgcna, dynamicColors, cutHeight = MEDissThres, 
                           verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <- merge$newMEs


# ==============================================================================
#
#  Code chunk 10: Rebuild the dendrogram
#
# ==============================================================================
message("Rebuilding the dendrogram")

png("gene_dendro.png")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ==============================================================================
#
#  Code chunk 11: Save the data 
#
# ==============================================================================

moduleColors <- mergedColors
MEs <- mergedMEs
save(MEs, moduleColors, file = "miRNA-network.RData")
