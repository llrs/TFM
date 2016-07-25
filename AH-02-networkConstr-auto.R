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
load(file = "InputWGCNA.RData", verbose = TRUE) 
#The variable lnames contains the names of loaded variables.


# ==============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================


# Choose a set of soft-thresholding powers
powers = c(1:30)
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
#  Code chunk 4: Plot how the modules correlate in a dendrogram
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
#  Code chunk 5: Save data for the next step
#
# ==============================================================================

gm <- table(net$colors)
hist(gm)
count.p <- function(data, per){
  # Calculate the percentatge above per of data
  sum(data >= per)/length(data)
}
perc <- unlist(lapply(gm, count.p, data = gm))
pdfn(file = "modules_distr.pdf")
plot(cbind(gm[order(perc)],perc[order(perc)]), type = "o",  
     xlab = "Size of the modules", 
     ylab = "Proportion of modules above the size",
     main = "Distribution of the size of the modules",
     col = names(gm[order(perc)]))


dev.off()
moduleColors = net$colors
MEs = net$MEs 
geneTree = net$dendrograms[[1]] 
dendro <- net$dendrograms
save(MEs, moduleColors, file = "TNF_AH-network-auto.RData")
