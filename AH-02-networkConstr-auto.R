# ==============================================================================
#
#  Code chunk 1: Starting from the preiously saved data
#
# ==============================================================================


library(WGCNA)
options(stringsAsFactors = FALSE) 
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
enableWGCNAThreads(6)
# Load the data saved in the first part
lnames <- load(file = "Input.RData") 
#The variable lnames contains the names of loaded variables.
lnames


# ==============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
# Call the network topology analysis function
sft <- pickSoftThreshold(data.wgcna, powerVector = powers, verbose = 5, 
                        networkType = "signed")
# Plot the results:
pdf(file = "Power_calculations", width = 9, height = 5)
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
dev.off()

# ==============================================================================
#
#  Code chunk 3: Automatic blocks creation using the power calculated
#
# ==============================================================================

print(paste("Recomended power", sft$powerEstimate))
net <- blockwiseModules(data.wgcna, power = sft$powerEstimate, 
                       TOMType = "signed", minModuleSize = 30, 
                       maxBlockSize = 8000, networkType = "signed", 
                       pamRespectsDendro = FALSE, 
                       saveTOMs = TRUE, 
                       saveTOMFileBase = "AH_signed", 
                       verbose = 3)

save(net, file = "net.RData")
# load("net.RData")
# ==============================================================================
#
#  Code chunk 4: Plot How the modules correlate in a dendrogram
#
# ==============================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module colors of first dendrogram", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)


# ==============================================================================
#
#  Code chunk 5: Save data for the next step
#
# ==============================================================================



moduleColors = net$colors
MEs = net$MEs 
geneTree = net$dendrograms[[1]] 
dendro <- net$dendrograms
save(MEs, moduleColors, file = "TNF_AH-network-auto.RData")
