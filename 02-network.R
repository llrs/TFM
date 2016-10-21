# Sarting ####
source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

# Input can be from those individual projects
load("../../isa_HA/unsigned_signed/Input.RData", verbose = TRUE)
isa.exprs <- data.wgcna
isa.disease <- vclin
load("../../silvia_HA/unsigned_signed/Input.RData", verbose = TRUE)
silvia.exprs <- data.wgcna
silvia.disease <- vclin
# Creating the multiData/multiExprs/... object
data.wgcna <- multiSet(isa = isa.exprs, silvia = silvia.exprs)
# Should list2multiData be used?
                         # No, it doesn't check or do anything else than messing

# Check if the genes are comparable according to WGCNA
gsg <- goodSamplesGenesMS(data.wgcna)
if (!gsg$allOK) {
  stop("check your data")
}

# MergeMaid ####
# Check if the genes follow the same correlations / remove unwanted noise
# data.merge <- sapply(data.wgcna, function(x){t(x$data)})
# isa <- data.merge$isa
# silvia <- data.merge$silvia
# mergm <- mergeExprs(isa, silvia)
# corcor <- intCor(mergm)
# pdf("mergmaid.pdf")
# plot(mergm, xlab = names(mergm)[1], ylab = names(mergm)[2],
#      main = "Integrative correlation of the top gene",
#      col = 3, pch = 4)
# hist(corcor, main = "Integrative correlation coeficient")
#
# intcor <- intcorDens(mergm)
# plot(intcor)
# dev.off()
# save(intcor, corcor, file = "mergemaid.RData")
load("../filtered/mergemaid.RData", verbose = TRUE)
coef <- as.vector(corcor@pairwise.cors)
names(coef) <- rownames(corcor@pairwise.cors)
comp.genes <- names(coef)[coef > 0] # Threshold of comparison
discutibles.genes <- names(coef)[coef <= 0]

data.wgcna2 <- lapply(data.wgcna, function(x, keep) {
  x$data[, colnames(x$data) %in% keep]
}, keep = comp.genes)
names(data.wgcna2) <- names(data.wgcna)
for (i in 1:length(data.wgcna)) {
  data.wgcna[[i]]$data <- data.wgcna2[[i]]
}

chSet <- checkSets(data.wgcna)
nGenes <- chSet$nGenes
nSamples <- chSet$nSamples
nSets <- chSet$nSets
save(data.wgcna, file = "Input.RData")
# bio.cor? ####
if (bio.corFnc) {
  bio_mat <- tryCatch({load("bio_correlation.RData")},
           warning = function(x){
             bio_mat <- bio.cor2(colnames(data.wgcna), ids = "Symbol",
                                 react = TRUE)
             save(bio_mat, file = "bio_correlation.RData")
             return(bio_mat)
           })
}


# power ####

# Choose a set of soft-thresholding powers
# powerTables <- vector(mode = "list", length = nSets)
# # Call the network topology analysis function for each set in turn
# for (set in 1:nSets) {
#   # Calculate the appropiate sft ####
#   if (bio.corFnc) {
#     sft <- pickSoftThreshold(data.wgcna[[set]]$data,
#                              powerVector = powers,
#                              verbose = 5,
#                              networkType = adj.opt,
#                              corFnc = cor.all,
#                              corOptions = list(bio_mat = bio_mat,
#                                                w = c(0.5, 0.5)))
#   } else {
#     sft <- pickSoftThreshold(data.wgcna[[set]]$data,
#                              powerVector = powers, verbose = 5,
#                              # corFnc = bicor,
#                              corOptions = list(nThreads = 6), #, maxPOutliers = 0.05),
#                              networkType = adj.opt)
#   }
#   # Set it as originally
#   powerTables[[set]] <- list(data = sft[[2]])
# }
# collectGarbage()
# save(powerTables, file = "powers_multiSet.RData")
load(file = "powers_multiSet.RData", verbose = TRUE)
power <- mean(multiple.softThreshold(powerTables))
pdf("Network_building.pdf")
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity",
             "Median connectivity", "Max connectivity")
# Get the minima and maxima of the plotted points
ylim <- matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets) {
  for (col in 1:length(plotCols)) {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]],
                       na.rm = TRUE)
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]],
                       na.rm = TRUE)
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

pars <- par(mfcol = c(2,2), mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 <- 0.7
for (col in 1:length(plotCols)) {
  for (set in 1:nSets) {
    if (set == 1) {
      plot(powerTables[[set]]$data[,1],
           -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           xlab = "Soft Threshold (power)",ylab = colNames[col],type = "n",
           ylim = ylim[, col],
           main = colNames[col])
      addGrid()
    }
    if (col == 1) {
      text(powerTables[[set]]$data[,1],
           -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
           labels = powers, cex=cex1,col=colors[set]);
    } else
      text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
           labels = powers,cex=cex1,col=colors[set]);
    if (col == 1) {
      legend("bottomright", legend = names(data.wgcna), col = colors, pch = 20)
    } else {
      legend("topright", legend = names(data.wgcna), col = colors, pch = 20)
    }
  }
}
dev.off()

message(paste("Recomended power", power))
if (is.na(power)) {
  stop("Estimated power, is NA\nReview the power manually!")
} else if (1/sqrt(nGenes) ^ power * nGenes >= 0.1) {
  warning("Are you sure of this power?")
}
message(paste("Using power", power))

# Not able to calculate on the multSet object
# Calculate connectivity and plot it
# k <- softConnectivity(data.wgcna, type = adj.opt, power = power,
#                       # corFnc = "bicor",
#                       )
# plot(density(k))
# scaleFreePlot(k, main = paste0("Check scale free topology, power",
#                                power))


# Network construction ####
net <- blockwiseConsensusModules(data.wgcna,
                        power = power,
                TOMType = TOM.opt,
                networkType = adj.opt,
                # corType = "bicor",
                # maxPOutliers = 0.05
                minModuleSize = 30,
                maxBlockSize = 8000,
                pamRespectsDendro = FALSE,
                saveTOMs = TRUE,
                saveTOMFileBase = "TOM",
                verbose = 3)

save(net, file = "net.RData")
load("net.RData", verbose = TRUE)

# Dendro ####

# pdf(file = "dendro.pdf", width = 12, height = 9)
# # Convert labels to colors for plotting
# mergedColors <- net$colors
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors of first dendrogram",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()

# # Connectivity ####
# connect <- intramodularConnectivity.fromExpr(data.wgcna, colors = net$colors,
#                                              networkType = adj.opt,
#                                              power = power,
#                                              scaleByMax = TRUE)
# save(connect, file = "kIM.RData")
# load(file = "kIM.RData", verbose = TRUE)


# Calculate eigengenes, it is already calculated
MEs <- consensusOrderMEs(net$multiMEs)

# Module exploring ####

# Calculate dissimilarity of module eigengenes
MEDiss <- consensusMEDissimilarity(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(8,10)
pdf("Modules_relationship.pdf")

plotEigengeneNetworks(MEs, names(data.wgcna))
par(pars)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")

labeledHeatmap(MEDiss,
               xLabels = substring(colnames(MEs), 3),
               yLabels = substring(colnames(MEs), 3),
               xSymbols = substring(colnames(MEs), 3),
               ySymbols = substring(colnames(MEs), 3))

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

moduleColors <- net$colors
# a <- sapply(unique(moduleColors), function(x){
#   p <- module.expr(data.wgcna, moduleColors, x)
#   ggsave(filename = name.file("module", x, ".png"),
#          plot = p)
# })

# save data ####
save(MEs, moduleColors, file = "Consensus-module_MEs.RData")
