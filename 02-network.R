# ==============================================================================
#
#  Code chunk 1: Starting from the previously saved data
#
# ==============================================================================

source("00-general.R")

# Load the data saved in the first part
# load(file = "InputWGCNA.RData", verbose = TRUE)
load(file = "shared_genes.RData", verbose = TRUE)
#The variable lnames contains the names of loaded variables.

# ==============================================================================
#
#  Code chunk 1b: Calculate the biological information of the genes
#
# ==============================================================================


# bio_mat <- bio.cor(colnames(data.wgcna))
# save(bio_mat, file = "bio_correlation.RData")

# ==============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================

# Choose a set of soft-thresholding powers
powers = c(1:30)
# Call the network topology analysis function
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, verbose = 5,
                        networkType = "unsigned",
                        # corFnc = cor.all,
                        # corOptions = list(use = "p", bio_mat = bio_mat)
                        )
# Plot the results:
pdfn(file = "Power_calculations.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, unsigned R^2", type = "n",
     main = "Scale independence", ylim = c(0, 1))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red", ylim = c(0, 1))
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,
     col = "red")
dev.off()

# ==============================================================================
#
#  Code chunk 3: Automatic blocks creation using the power calculated
#
# ==============================================================================

print(paste("Recomended power", sft$powerEstimate))
if (is.na(sft$powerEstimate)){
  stop("Review the power manually")
}
net <- blockwiseModules(data.wgcna, power = 9, # the max connectivity of 0.73
                # Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
                #    9   0.7250 -1.180          0.660     389    148.00   1970
                TOMType = "unsigned",
                minModuleSize = 30,
                maxBlockSize = 8000,
                networkType = "unsigned",
                pamRespectsDendro = FALSE,
                saveTOMs = TRUE,
                saveTOMFileBase = "AH_unsig.unsig",
                verbose = 3,
                nThreads = nThreads)

save(net, file = "net_unsigned.RData")
load("net_unsigned.RData")
# ==============================================================================
#
#  Code chunk 4: Plot how the modules correlate in a first dendrogram
#
# ==============================================================================


# open a graphics window
pdfn(file = "dendro_unsig.pdf", width = 12, height = 9)
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
pdfn("Modules_relationship_unsig.pdf")
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

perc <- unlist(lapply(gm, count.p, data = gm))
plot(cbind(gm[order(perc)],perc[order(perc)]), type = "o",
     xlab = "Size of the modules",
     ylab = "Proportion of modules above the size",
     main = "Distribution of the size of the modules",
     col = names(gm[order(perc)]))

# Call an automatic merging function
# Without merging from now on!!
# merge <- mergeCloseModules(data.wgcna,
#                           net$colors, cutHeight = MEDissThres, MEs = MEs,
#                           verbose = 3)
# # The merged module colors
# mergedColors <- merge$colors;
# # Eigengenes of the new merged modules:
# MEs <- merge$newME
#
# plot(merge$dendro, main = "After merging the modules")

dev.off()
# ==============================================================================
#
#  Code chunk 6: Save the data
#
# ==============================================================================

moduleColors <- mergedColors

save(MEs, moduleColors, file = "TNF_AH-network-unsig.RData")
