# 1 ============================================================================
#
#  Code chunk 1: Starting from the previously saved data
#
# ==============================================================================

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)
# Load the data saved in the first part
load(file = "Input.RData", verbose = TRUE)

nGenes <- ncol(data.wgcna)
nSamples <- nrow(data.wgcna)

# 1b ===========================================================================
#
#  Code chunk 1b: Calculate the biological information of the genes
#
# ==============================================================================

if (bio.corFnc) {
  tryCatch({load("bio_correlation.RData")},
           error = function(x){
             bio_mat <- bio.cor2(colnames(data.wgcna), ids = "Symbol",
                                 react = TRUE)
             save(bio_mat, file = "bio_correlation.RData")
           })
}


# 2 ============================================================================
#
#  Code chunk 2: Deciding the power to use
#
# ==============================================================================

# Choose a set of soft-thresholding powers

# Call the network topology analysis function
if (bio.corFnc) {
  sft <- pickSoftThreshold(data.wgcna,
                           powerVector = powers,
                           verbose = 5,
                           networkType = adj.opt,
                           corFnc = cor.all, corOptions(bio_mat = bio_mat,
                           w = c(0.5, 0.5)))
} else {
  sft <- pickSoftThreshold(data.wgcna,
                           powerVector = powers, verbose = 5,
                           networkType = adj.opt)
}


# load("sft.RData", verbose = TRUE)
cex1 <- 0.9
# Plot the results:
pdfn(file = "Network_building.pdf")
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, R^2", type = "n",
     main = "Scale independence", ylim = c(0, 1))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red", ylim = c(0, 1))
# this line corresponds to using an R^2 cut-off of h
abline(h = c(0.90, 0.85), col = c("red", "green"))
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1,
     col = "red")
abline(h = c(100, 1000), col = c("green", "red"))


print(paste("Recomended power", sft$powerEstimate))
if (is.na(sft$powerEstimate)) {
  stop("Estimated power, is NA\nReview the power manually!")
} else if (1/sqrt(nGenes) ^ sft$powerEstimate * nGenes >= 0.1) {
  warning("Are you sure of this power?")
}
sft$powerEstimate <- 3
print(paste("Using power", sft$powerEstimate))
save(sft, file = "sft.RData")

# Calculate connectivity and plot it
k <- softConnectivity(data.wgcna, type = adj.opt, power = sft$powerEstimate)
plot(density(k))
scaleFreePlot(k, main = paste0("Check scale free topology, power",
                               sft$powerEstimate))
dev.off()

# 3 ============================================================================
#
#  Code chunk 3: Automatic blocks creation using the power calculated
#
# ==============================================================================
net <- blockwiseModules(data.wgcna,
                        power = sft$powerEstimate,
                TOMType = TOM.opt,
                networkType = adj.opt,
                minModuleSize = 30,
                maxBlockSize = 8000,
                pamRespectsDendro = FALSE,
                saveTOMs = TRUE,
                saveTOMFileBase = "TOM",
                verbose = 3)

save(net, file = "net.RData")
load("net.RData", verbose = TRUE)

# 4 ============================================================================
#
#  Code chunk 4: Plot how the modules correlate in a first dendrogram
#
# ==============================================================================


pdf(file = "dendro.pdf", width = 12, height = 9)
# Convert labels to colors for plotting
mergedColors <- net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors of first dendrogram",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 5 ============================================================================
#
#  Code chunk 5: Explore the connectivity
#
# ==============================================================================
#
connect <- intramodularConnectivity.fromExpr(data.wgcna, colors = net$colors,
                                             networkType = adj.opt,
                                             power = sft$powerEstimate,
                                             scaleByMax = TRUE)
save(connect, file = "kIM.RData")
load(file = "kIM.RData", verbose = TRUE)


# Calculate eigengenes, it is already calculated
MEs <- net$MEs
MEs <- orderMEs(MEs)

# It is the same as MM == kME
kME <- signedKME(data.wgcna, MEs)
save(kME, file = "kME.RData")
# load("kME.RData", verbose = TRUE)

# 5b ===========================================================================
#
#  Code chunk 5b: See how are the modules found
#
# ==============================================================================


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

labeledHeatmap(MEDiss, xLabels = colnames(corME), yLabels = colnames(corME),
               xSymbols = colnames(corME), ySymbols = colnames(corME))

gm <- table(net$colors)
gm
gm <- gm[names(gm) != "grey"]
hist(gm, xlab = "Size of modules")

perc <- unlist(lapply(gm, count.p, data = gm))
plot(cbind(gm[order(perc)], perc[order(perc)]), type = "o",
     xlab = "Size of the modules",
     ylab = "Proportion of modules above the size",
     main = "Distribution of the size of the modules",
     col = names(gm[order(perc)]))

dev.off()

a <- sapply(unique(moduleColors), function(x){
  p <- module.expr(data.wgcna, moduleColors, x)
  ggsave(filename = name.file("module", x, ".png"),
         plot = p)
})

# 6 ============================================================================
#
#  Code chunk 6: Save the data for the next process
#
# ==============================================================================

moduleColors <- net$colors

save(MEs, moduleColors, file = "modules_ME.RData")
