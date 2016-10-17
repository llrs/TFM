# Input files should be modified to address:
#  - more studies combined (Normalizations)
#  - more clinical information (This should be unified before)

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)

# Load input data ####
# Phenodata contains the information of the experiment space separated!
# not tab separated

pheno.isa <- read.affy(pheno1)
# pheno.silvia <- read.affy(pheno2)

# Download set ####
# geosupp <- getGEOSuppFiles(gse_number)
# geosupp
# Unpack the CEL files
# untar(path.raw, exdir = study.dir)
# cels <- list.files(study.dir, pattern = "[gz]")
# sapply(cels, gunzip)
# cels  <-  list.files(study.dir, pattern = "CEL")

# disease.silvia <- read.csv("clean_variables.csv")
disease.isa <- read.csv("samples_AH.csv")

setwd(data.files.out)
# save(pheno.isa, pheno.silvia, file = "pheno.RData")
# load("pheno.RData", verbose = TRUE)

# Adjusting intensity and making them comparable ####
c.isa <- rma(pheno.isa)
# c.silvia <- rma(pheno.silvia)
# save(c.isa, c.silvia, file = "rma.pheno.RData")
# load("rma.pheno.RData")
co.isa <- sum.e(c.isa)
# co.silvia <- sum.e(c.silvia)
# save(co.isa, co.silvia, file = "exprs.RData")
# load("exprs.RData", verbose = TRUE)

# Merging datasets ####
# Merge the data of each batch into a single matrix
# co.silvia.df <- as.data.frame(t(co.silvia), row.names = colnames(co.silvia))
# co.isa.df <- as.data.frame(t(co.isa), row.names = colnames(co.isa))
# merged <- rbind.fill(co.silvia.df, co.isa.df)
# rownames(merged) <- c(colnames(co.silvia), colnames(co.isa))
# save(merged, file = "collapsed.micro.RData")
# load("collapsed.micro.RData", verbose = TRUE)

# with.na <- apply(merged, 2, function(x){any(is.na(x))})
# merged.shared <- merged[, !with.na] # Keep the shared genes
# merged.pca <- t(merged)
# merged.shared.pca <- t(merged.shared)
# data.wgcna <- merged.shared[1:15, ]
# save(data.wgcna, file = "shared_genes.RData")
# load(file = "shared_genes.RData", verbose = TRUE)

# Merging with MergeMaid ####

# mergm <- mergeExprs(co.silvia, co.isa)
# corcor <- intCor(mergm)
# pdf("mergmaid.pdf")
# plot(mergm, xlab = names(mergm)[1], ylab = names(mergm)[2],
#      main = "Integrative correlation of the top gene",
#      col = 3, pch = 4)
# hist(corcor, main = "Integrative correlation coeficient")
#
# intcor <- intcorDens(mergm)
# plot(intcor)
# # cox.coeff <- modelOutcome(mergm, outcome = c(3, 3), # Obscure parameter
# #                           method = "linear")
# # plot(coeff(cox.coeff), main = "Coeficients")
# dev.off()
# save(intcor, corcor, file = "mergemaid.RData")
# load("mergemaid.RData", verbose = TRUE)

# coef <- as.vector(corcor@pairwise.cors)
# names(coef) <- rownames(corcor@pairwise.cors)
# comp.genes <- names(coef)[coef > 0] # Threshold of comparison
# discutibles.genes <- names(coef)[coef <= 0]

# Store the procedence of the data
# orig.data <- as.factor(c(rep("Silvia", ncol(co.silvia)),
#                          rep("Isa", ncol(co.isa))))
# new.subset <- merged.shared.pca[rownames(merged.shared.pca) %in% comp.genes, ]

# pca.graph(data = new.subset, file = "merged.shared_filtered.pca.pdf",
#           col = as.numeric(orig.data),
#           outcome = orig.data)
#
# pca.graph(data = merged.shared.pca, file = "merged.shared.pca.pdf",
#           col = as.numeric(orig.data),
#           outcome = orig.data)

# Code used for the DE analysis ######
# set colour palette
# cols <- brewer.pal(8, "Set1")
#
# # plot a boxplot of unnormalised and normalized intensity values
#
# png("normalization_boxplot.png", width = 1000, height = 750)
# pars <- par(mfrow = c(1,3), mar = c(4, 5, 3, 2))
# boxplot(celfiles, col = cols, title = "Unnormalised")
# boxplot(c.gcrma, col = cols, title = "Gcrma normalization")
# boxplot(c.rma, col = cols, titlw = "RMA normalization")
# dev.off()
#
# # Plot a density vs log intensity histogram for the unnormalised and normalised data
# png("normalization_histogram.png", width = 1500, height = 750)
# par(mfrow = c(1,3), mar = c(4, 5, 2, 2))
# hist(celfiles, col = cols, main = "without")
# hist(c.gcrma, col = cols, main = "gcrma") # This seem to be the best normalization algorithm
# hist(c.rma, col = cols, main = "rma")
# dev.off()
#
# # MA plots to see how well they are
# # add option plot.method = "smoothScatter" to get a fancier plot
# png("normalization_maplot.png", width = 1500, height = 750)
# par(mfrow = c(1,3), mar = c(4, 5, 2, 2))
# MAplot(celfiles, cex = 0.75, ref.title = "Raw")
# MAplot(c.gcrma, cex = 0.75, ref.title = "Gcrma")
# MAplot(c.rma, cex = 0.75, ref.title = "rma")
# dev.off()
# par(pars)
#
# # look at the relationships between the samples using heirarchical clustering
# eset <- exprs(c.gcrma)
# methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
# pars <- par(mfrow = c(2,3))
# pdf("hierarchical clustering.pdf")
# for (meth in methods) {
#   distance <- dist(t(eset), method = meth)
#   cluster <- hclust(distance)
#   plot(cluster, labels = phenoData(celfiles)$Target, main = paste(meth, "distance"))
# }
# par(pars)
#
# # Depending on the type of distance the clustering goes different
#
# # Check for batch effect (based on day performed)
# scanDate <- protocolData(c.gcrma)$ScanDate
# scanDate <- gsub(" .*", "", scanDate)
# scanDate <- as.Date(scanDate, "%m/%d/%Y")
# minscan <- min(scanDate)
# days <- scanDate - minscan
#
# days[days == 0] <- 1
# days[days == 2] <- 2
# days[days == 8] <- 3
# days[days == 13] <- 4
# batch <- as.numeric(days)
#
# # Most of the AH mirocarrays where done on 2 batches,
# # one of them without control, bias?
# table(data.frame(Outcome = c.gcrma$Type, Batch = batch))
#
# # Performs the clustering taking into account the date it was performed.
# d.celfiles <- as.dist(1 - cor(exprs(celfiles), method = "spearman"))
# d.gcrma <- as.dist(1 - cor(exprs(c.gcrma), method = "spearman"))
# sampleClustering <- hclust(d.gcrma)
# sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
# names(batch) <- phenoData(celfiles)$Type
# outcome <- as.character(phenoData(celfiles)$Type)
# names(outcome) <- sampleNames(celfiles)
# sampleDendrogram_c <- dendrapply(sampleDendrogram, function(x, batch, labels) {
#   ## for every node in the dendrogram if it is a leaf node
#   if (is.leaf(x)) {
#     attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))
#     ## color by batch
#     attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
#   }
#   x
# }, batch, outcome) ## these are the second and third arguments in the function
#
# png("hierarchical.png")
# # Plot dendogram with the information of the dates.
# plot(sampleDendrogram_c, main = "Hierarchical clustering of samples")
# legend("bottom", cex = 0.75,
#        paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))
# dev.off()
#
# # Multidimension scaling (PCA)
# pdf("new_PCA.pdf", onefile = TRUE)
# cmd <- cmdscale(d.celfiles)
# plot(cmd, type = "n", main = "PCA without normalization")
# text(cmd, outcome, col = batch, cex = 0.9)
# legend("top", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)
# # Doesn't show any problem between batch
# cmd <- cmdscale(d.gcrma)
# plot(cmd, type = "n", main = "PCA with gcrma normalization")
# text(cmd, outcome, col = batch, cex = 0.9)
# legend("top", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)
#
#

# Batch effect correction
# Quantifying the counfounding factor
# pdf("PCA.counfounding.merged.pdf")
# s <- fast.svd(t(scale(merged.shared.pca, center = TRUE, scale = TRUE)))
# PCA <- s$d ^ 2/sum(s$d ^ 2)
#
# plot(PCA, type = "b", lwd = 2, las = 1,
     # xlab = "Principal Component", ylab = "Proportion of variance",
     # main = "Principal components contributions")
# dev.off()
# save(c.gcrma, file = "corrected_exprs.RData")
# dev.off()
#

# ComBat ####

# combat.exp <- ComBat(merged.shared.pca, orig.data)

# pca.graph(data = combat.exp, file = "merged.shared.pca.combat.pdf",
#           col = as.numeric(orig.data),
#           outcome = as.character(orig.data))

# combat.subset <- ComBat(new.subset, orig.data)

# pca.graph(data = combat.subset, file = "merged.shared_filtered.pca.combat.pdf",
#           col = as.numeric(orig.data),
#           outcome = as.character(orig.data))

# QR decomposition ####

# qrexp <- removeBatchEffect(merged.shared.pca, orig.data)
# pca.graph(data = qrexp, file = "merged.shared.pca.rBE.pdf",
#           col = as.numeric(orig.data),
#           outcome = as.character(orig.data))

# Mergemaid after using Combat to make them comparable ####
# isa <- combat.exp[ ,1:15]
# silvia <- combat.exp[, 16:ncol(combat.exp)]
#
# mergm <- mergeExprs(isa, silvia)
# corcor <- intCor(mergm)
#
# pdf("mergmaid_combat.pdf")
# plot(mergm, xlab = names(mergm)[1], ylab = names(mergm)[2],
#      main = "Integrative correlation",
#      col = 3, pch = 4)
# hist(corcor, main = "Integrative correlation coeficient")
#
# intcor <- intcorDens(mergm)
# # cox.coeff <- modelOutcome(mergm, outcome = c(3, 3), # Obscure parameter
# #                           method = "linear")
# # plot(coeff(cox.coeff), main = "Coeficients")
# dev.off()
# save(intcor, corcor, file = "mergemaid_combat.RData")
# load("mergemaid_combat.RData", verbose = TRUE)
#
# coef <- as.vector(corcor@pairwise.cors)
# names(coef) <- rownames(corcor@pairwise.cors)
# comp.genes <- names(coef)[coef > 0]
# discutibles.genes <- names(coef)[coef <= 0]

# Prepare the variables to the right format ####

colnames(disease.isa) <- tolower(colnames(disease.isa))
clin.isa <- cbind("files" = rownames(pData(pheno.isa)),
                  pData(pheno.isa))
# clin.silvia <- cbind("files" = rownames(pData(pheno.silvia)),
#                      pData(pheno.silvia))
disease.isa <- rename.col(disease.isa, "codi_pacient", "id")
disease.isa <- rename.col(disease.isa, 'creat', 'creatinine')
disease.isa <- rename.col(disease.isa, 'leuc', 'leucos')
disease.isa <- rename.col(disease.isa, 'alb', 'albumin')
disease.isa <- rename.col(disease.isa, 'triglicerids', 'trigyicerides')
disease.isa <- rename.col(disease.isa, 'glucosa', 'glucose')
disease.isa <- rename.col(disease.isa, 'viu_3m', 'status_90')
disease.isa$status_90 <- fact2num(disease.isa$status_90, "si", "alive")
disease.isa$status_90 <- fact2num(disease.isa$status_90, "no", "exitus")

disease.isa$plaq <- disease.isa$plaq*1000
disease.isa <- rename.col(disease.isa, 'plaq', 'platelets')
disease.isa$hb <- disease.isa$hb/10
# disease.isa <- rename.col(disease.isa, 'hb', 'hb_g.dl')

# clin <- rbind.fill(clin.isa, clin.silvia)
# disease <- rbind.fill(disease.silvia, disease.isa)
disease <- disease.isa
clin <- clin.isa
disease <- as.data.frame(apply(disease, 2, trim))

disease$aki <- fact2num(disease$aki, "(yes)|(si)", 1)
disease$aki <- fact2num(disease$aki, "no", 0)
disease$aki <- level.na(disease$aki)

disease$status_90 <- fact2num(disease$status_90, "alive", 1)
disease$status_90 <- fact2num(disease$status_90, "exitus", 0)
disease$status_90 <- level.na(disease$status_90)

disease$infection_hospitalization <- fact2num(disease$status_90, "yes", 0)
disease$infection_hospitalization <- fact2num(disease$status_90, "no", 1)
disease$status_90 <- level.na(disease$status_90)

vclin <- merge(clin, disease, by.x = "Sample", by.y = "id", all.x = TRUE)
int.Var <- c("Sample", "files", "meld", "maddrey", "lille_corte", "lille",
            "status_90", "glucose",
             "trigyicerides", "ast", "alt", "bili_total", "creatinine",
             "albumin", "inr", "ggt", "ap", "leucos", "hb_g.dl", "hematocrit",
             "platelets", "hvpg_corte20", "hvpg", "aki",
             "infection_hospitalization")
vclin <- vclin[, colnames(vclin) %in% int.Var]
# 1 above, 0 below
if ("hvpg" %in% colnames(vclin)) {
  vclin$hvpg_corte20 <- as.numeric(as.numeric(vclin$hvpg) > 20)
}
# 1 below, 0 above
if ("lille" %in% colnames(vclin)) {
  vclin$lille_corte <- as.numeric(as.numeric(vclin$lille) < 0.45)
}

data.wgcna <- t(co.isa)
# Filtering those with low correlation between studies
# data.wgcna <- data.wgcna[ , colnames(data.wgcna) %in% comp.genes]
# Changing the name of the files by the samples name
rownames(data.wgcna) <- vclin$Sample[match(rownames(data.wgcna), vclin$files)]
gsg <- goodSamplesGenes(data.wgcna, verbose = 3)

if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing", sum(!gsg$goodGenes),"genes"))
    # printFlush(paste("Removing genes:",
    #                  paste(names(data.wgcna)[!gsg$goodGenes],
    #                        collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing", sum(!gsg$goodSamples),"samples."))
    # printFlush(paste("Removing samples:",
    #                  paste(rownames(data.wgcna)[!gsg$goodSamples],
    #                        collapse = ", ")))
  }
  # Remove the offending genes and samples from the data:
  data.wgcna <- data.wgcna[gsg$goodSamples, gsg$goodGenes]
}

# samples.pre <- t(merged.shared.pca)
# rownames(samples.pre) <- vclin$Sample[match(rownames(samples.pre), vclin$files)]
pdf("Samples.pdf")
# sampleTree <- hclust(dist(samples.pre), method = "average")
# plot(sampleTree, main = "Sample clustering to detect outliers",
#      sub = "Before ComBat correction", xlab = "", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
sampleTree <- hclust(dist(data.wgcna), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "After corrections", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


# Re-cluster samples
sampleTree2 <- hclust(dist(data.wgcna),
                      method = "average")
# Convert traits to a color representation: white means low, red means high,
# grey means missing entry
v <- apply(vclin[, 3:ncol(vclin)], 2, as.numeric)
no.keep <- apply(v, 2, function(x){all(is.na(x))})
v <- v[, !no.keep]
rownames(v) <- vclin$Sample
traitColors <- numbers2colors(v, signed = FALSE)
dimnames(traitColors) <- dimnames(v)
# Plot the sample dendrogram and the colors underneath.
pdf("Dendro_traits.pdf")
plotDendroAndColors(sampleTree2, traitColors,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Keep those who are different to 1
rm.i <- apply(v, 2, function(x){length(unique(x[!is.na(x)]))}) == 1
v <- v[, !rm.i]

vclin <- v
vclin <- orderby(vclin, rownames(data.wgcna), names.x = TRUE)
save(data.wgcna, vclin, file = "Input.RData")
