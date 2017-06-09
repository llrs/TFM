# starting ####

source("~/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)
getwd()
consensus <- FALSE # If it from a consensus file
power <- TRUE # Calculate the power or reuse the existing in the folder
network <- TRUE # Build the network or reuse the existing in the folder
dendro <- TRUE # plot a dendro or not
connectivity <- FALSE # If consensus the connectivity shouldn't be calculated
# Input can be from those individual projects
if (consensus) {
  load("../Early_Network.RData", verbose = TRUE)
  isa.exprs <- data.wgcna
  isa.disease <- vclin
  load("../Late_Network.RData", verbose = TRUE)
  silvia.exprs <- data.wgcna
  silvia.disease <- vclin
  # Creating the multiData/multiExprs/... object
  data.wgcna <- multiSet(early = isa.exprs, late = silvia.exprs)
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
  # load("../filtered/mergemaid.RData", verbose = TRUE)
  # coef <- as.vector(corcor@pairwise.cors)
  # names(coef) <- rownames(corcor@pairwise.cors)
  # comp.genes <- names(coef)[coef > 0] # Threshold of comparison
  # discutibles.genes <- names(coef)[coef <= 0]
  # data.wgcna2 <- lapply(data.wgcna, function(x, keep) {
  #   x$data[, colnames(x$data) %in% keep]
  # }, keep = comp.genes)
  # names(data.wgcna2) <- names(data.wgcna)
  # for (i in 1:length(data.wgcna)) {
  #   data.wgcna[[i]]$data <- data.wgcna2[[i]]
  # }

  chSet <- checkSets(data.wgcna)
  nGenes <- chSet$nGenes
  nSamples <- chSet$nSamples
  nSets <- chSet$nSets
  save(data.wgcna, file = "Input.RData")
} else {
  # Load the data saved in the first part
  load(file = "../../ASH_AH.RData", verbose = TRUE)
  nGenes <- ncol(data.wgcna)
  nSamples <- nrow(data.wgcna)
}

# bio.cor ####
if (bio.corFnc) {
  if (file.exists("~/Documents/geneSim.RData")) {
    load("~/Documents/geneSim.RData", verbose = TRUE)
    entrez <- mapIds(org.Hs.eg.db, keys = colnames(data.wgcna),
                     keytype = "SYMBOL", column = "ENTREZID")
    m <- matrix(ncol = ncol(data.wgcna), nrow = ncol(data.wgcna),
                dimnames = list(entrez, entrez))
    bio_mat <- list(react = AintoB(o, m))
  } else {
    bio_mat <- bioCor(colnames(data.wgcna), ids = "SYMBOL",
                        all = TRUE, BPPARAM = SnowParam())
    save(bio_mat, file = "bio_correlation.RData")
  }
}

# power ####
if (power) {
  if (consensus) {
    # Choose a set of soft-thresholding powers
    powerTables <- vector(mode = "list", length = nSets)
    # Call the network topology analysis function for each set in turn
    for (set in 1:nSets) {
      # Calculate the appropiate sft ####
      sft <- pickSoftThreshold(data.wgcna[[set]]$data,
                               powerVector = powers,
                               verbose = 5,
                               corFnc = bicor,
                               corOptions = list(maxPOutliers = 0.10),
                               networkType = adj.opt,
                               RsquaredCut = 0.80)
      # Set it as originally
      powerTables[[set]] <- list(data = sft[[2]])
    }
    collectGarbage()
    save(powerTables, file = "powers_multiSet.RData")
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
               xlab = "Soft Threshold (power)", ylab = colNames[col], type = "n",
               ylim = ylim[, col],
               main = colNames[col])
          addGrid()
        }
        if (col == 1) {
          text(powerTables[[set]]$data[,1],
               -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
               labels = powers, cex = cex1, col = colors[set])
        } else
          text(powerTables[[set]]$data[,1],
               powerTables[[set]]$data[,plotCols[col]],
               labels = powers, cex = cex1, col = colors[set])
        if (col == 1) {
          legend("bottomright", legend = names(data.wgcna),
                 col = colors, pch = 20)
        } else {
          legend("topright", legend = names(data.wgcna),
                 col = colors, pch = 20)
        }
      }
    }
    dev.off()
  } else {
    # Call the network topology analysis function
    sft <- pickSoftThreshold(data.wgcna,
                             powerVector = powers, verbose = 5,
                             corFnc = bicor,
                             corOptions = list(nThreads = 6,
                                               maxPOutliers = 0.10),
                             networkType = adj.opt,
                             RsquaredCut = 0.80)
    save(sft, file = "sft.RData")
    power <- sft$powerEstimate
    # load("sft.RData", verbose = TRUE)
    cex1 <- 0.9
    # Plot the results:
    pdf(file = "Network_building.pdf")
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

  }

  message("The recommended power is ", power)
  if (is.na(power)) {
    stop("Estimated power, is NA\nReview the power manually!")
  } else if (1/sqrt(nGenes) ^ power * nGenes >= 0.1) {
    warning("Are you sure of this power?")
  } else  {
    opt1 <- (adj.opt == "signed hybrid" | adj.opt == "unsigned") & power >= 15
    opt2 <- adj.opt == "signed" & power >= 30
    if (opt1 | opt2) {
      warning("Check the FAQ question 6")
    }
  }


  # softConnectivity ####
  if (!consensus) {
    k <- softConnectivity(data.wgcna,
                          type = adj.opt,
                          power = power,
                          corFnc = "bicor",
                          corOptions = "maxPOutliers = 0.10")
    plot(density(k))
    scaleFreePlot(k, main = paste0("Check scale free topology, power",
                                   sft$powerEstimate))
    dev.off()
  }
} else {
  if (consensus) {
    load(file = "powers_multiSet.RData", verbose = TRUE)
    power <- mean(multiple.softThreshold(powerTables, min = 0.81))
    # power <- 10
  } else {
    # load("sft.RData", verbose = TRUE)
    # power <- sft$powerEstimate
    power <- 4 #manually so Rsft
  }
}
message(paste("Using power", power))


# Network construction ####
if (network) {
  if (bio.corFnc) {
    adj <- adjacency(data.wgcna, type = adj.opt, power = power)
    adj.bio <- addSimilarities(adj, bio_mat, weights = c(0.8, 0.2))
    TOM <- TOMsimilarity(adj.bio, TOMType = TOM.opt)
    dissTOM <- 1 - TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
    moduleColors <- labels2colors(dynamicMods$labels)
    MEs <- moduleEigengenes(data.wgcna, moduleColors)
    MEs <- MEs$eigengenes
  } else if (consensus) {
    net <- blockwiseConsensusModules(data.wgcna,
                                     power = power,
                                     TOMType = TOM.opt,
                                     networkType = adj.opt,
                                     corType = "bicor",
                                     maxPOutliers = 0.10,
                                     minModuleSize = 30,
                                     maxBlockSize = 8000,
                                     pamRespectsDendro = FALSE,
                                     saveTOMs = TRUE,
                                     saveTOMFileBase = "TOM",
                                     verbose = 3)
    save(net, file = "net.RData")
  } else {
    net <- blockwiseModules(data.wgcna,
                            power = power,
                            TOMType = TOM.opt,
                            networkType = adj.opt,
                            corType = "bicor",
                            maxPOutliers = 0.10,
                            minModuleSize = 30,
                            maxBlockSize = 8000,
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "TOM",
                            verbose = 3)
    save(net, file = "net.RData")
  }
} else {
  load("net.RData", verbose = TRUE)
}

if (connectivity) {
  # IntramodularConnectivity ####
  if (!consensus) {
    connect <- intramodularConnectivity.fromExpr(
      data.wgcna, colors = net$colors,
      networkType = adj.opt,
      power = power,
      scaleByMax = TRUE,
      corFnc = "bicor",
      corOptions = "maxPOutliers = 0.10")
    save(connect, file = "kIM.RData")
    load(file = "kIM.RData", verbose = TRUE)
  }
}
# Extract eigengenes ####
if (consensus) {
  MEs <- consensusOrderMEs(net$multiMEs)
} else if (!bio.corFnc) {
  MEs <- orderMEs(net$MEs)
  # Module exploring ####
  # It is the same as MM == kME
  kME <- signedKME(data.wgcna, MEs)
  save(kME, file = "kME.RData")
  # load("kME.RData", verbose = TRUE)
}

# Modules ####
# Calculate dissimilarity of module eigengenes
if (consensus) {
  MEDiss <- consensusMEDissimilarity(MEs)
} else {
  corME <- cor(MEs)
  MEDiss <- 1 - corME
}
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
if (!bio.corFnc) {
  moduleColors <- net$colors
}
pars <- par()
pdf("Modules_relationship.pdf")
if (consensus) {
  plotEigengeneNetworks(MEs, names(data.wgcna))
  par(pars)
} else {
  # Plot the dendrogram and the module colors underneath
  if (dendro & !bio.corFnc) {
    plotDendroAndColors(net$dendrograms[[1]],
                        moduleColors[net$blockGenes[[1]]],
                        "Module colors of first dendrogram",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
}

par(pars)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres <- 0.25
# Plot the cut line into the dendrogram
abline(h = MEDissThres, col = "red")
par(pars)
labeledHeatmap(MEDiss,
               xLabels = colnames(MEs),
               yLabels = colnames(MEs),
               xSymbols = colnames(MEs),
               ySymbols = colnames(MEs))

gm <- table(moduleColors)
gm
gm <- gm[names(gm) != "grey"]
hist(gm, xlab = "Size of modules")

perc <- unlist(lapply(gm, count.p, data = gm))
par(pars)
plot(cbind(gm[order(perc)], perc[order(perc)]), type = "o",
     xlab = "Size of the modules",
     ylab = "Proportion of modules above the size",
     main = "Distribution of the size of the modules",
     col = names(gm[order(perc)]))
dev.off()

# Plot the expression pattern among the samples
# a <- sapply(unique(moduleColors), function(x){
#   p <- module.expr(data.wgcna, moduleColors, x)
# #   p <- p + ggtitle("6 RD and 6 non RD")
#   ggsave(filename = name.file("module", x, ".png"),
#          plot = p)
# })

# save data ####
if (consensus) {
  save(MEs, moduleColors, file = "consensus-modules_MEs.RData")
} else {
  save(MEs, moduleColors, file = "modules_ME.RData")
}
