# Analyse data with WGCNA from doi:10.1136/gutjnl-2011-301146.

library("GEOquery")
library("affy")
library("affyPLM")
library("simpleaffy")
library("Heatplus")
# library("affyQCReport")
library("corpcor")
library("sva")
library("annotate")
library("hgu133plus2.db")
library("data.table")
library("ggbiplot")
library("corpcor")
library("WGCNA")
enableWGCNAThreads(6)
library("edgeR")
library("limma")
library("gcrma")
library("RColorBrewer")
library("hgu219.db")
library("plyr")
# Prepare some variables
gse_number <- "GSE28619"
path_files <- "../Documents/data/GSE28619_RAW/"
test_file <- "GSM709348.CEL.gz"
path_file <- paste0(path_files, test_file)
origDir <- setwd("~/Documents/data")
raw.tar <- paste0(gse_number, "_RAW.tar")
path_raw <- paste(gse_number, raw.tar, sep = "/")
#By exploring the data on GEO I found that the platform is :
# GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# So we need to load the bioconductor library with the information about it
library("Affyhgu133Plus2Expr")
library("hgu133plus2probe")
library("hgu133plus2cdf")

# Tutorial on https://www.biostars.org/p/53870/
# Based on this tutorial:
# http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
# Download set
# geosupp <- getGEOSuppFiles(gse_number)
# geosupp
#Unpack the CEL files
# untar(path_raw, exdir = "data")
# cels <- list.files("data/", pattern = "[gz]")
# sapply(paste("data", cels, sep = "/"), gunzip)
# cels  <-  list.files("data/", pattern = "CEL")
setwd("../data/")

# Phenodata contains the information of the experiment space separated! 
# not tab separated
warning("phenodata is manually created! And not under git!")
pheno.isa <- read.affy("pheno.isa.txt")
pheno.silvia <- read.affy("pheno.silvia.txt")
# stop("readfiles!")
setwd(origDir)

pca.graph <- function(celfiles=NULL, data=NULL, file, outcome = NULL, 
                      col = NULL, ...){
  # Normalize data and plots PCA of samples
  # Data is the normalized, celfiles are the raw files
  if (is.null(data)) {
    if (!is.null(celfiles)) {
      data <- rma(celfiles)
    } else {
      stop("celfiles or data must be provided")
    }
  }
  
  if ("ExpressionSet" %in% is(data)) {
    if (is.null(outcome)) {
      outcome <- as.character(phenoData(data)$Type)
    }
    dists <- as.dist(1 - cor(exprs(data), method = "spearman"))
  } else {
    if (is.null(outcome)) {
      outcome <- rep("AH", ncol(data))
    }
    dists <- as.dist(1 - WGCNA::cor(data, method = "spearman", 
                                    use = "pairwise.complete.obs"))
  }
  
  cmd <- cmdscale(dists)
  pdf(file)
  pca.plo <- ggbiplot(prcomp(dists, scale. = TRUE), 
                      obs.scale = 1, var.scale = 1, ellipse = TRUE,
                      group = outcome, var.axes = FALSE,
                      # circle = TRUE
                      )
  plot(pca.plo)
  plot(cmd, type = "n", main = "PCA samples", ...)
  text(cmd, outcome, col = col, cex = 0.9)
  dev.off()
  invisible(data)
}

sum.e <- function(eset){
  # Given a expression set transforms it to the gene symbols
  ann <- annotation(eset)
  pkg <- paste0(ann, ".db")
  pkg.ann <- eval(parse(text = pkg))
  symbols <- select(pkg.ann, keys(pkg.ann), "SYMBOL")
  expr <- exprs(eset)
  probes <- rownames(expr)
  rowGroup <- symbols[match(probes, symbols$PROBE),"SYMBOL"]
  # the length of rowGroup and probes should be the same
  # even if there is a warning we cannot omit it 
  # length(unique(rowGroup)) == length(probes.isa)
  corr.e <- collapseRows(expr, rowGroup = rowGroup, rowID = probes)
  corr.sy <- corr.e$datETcollapsed
  return(corr.sy)
}

c.isa <- rma(pheno.isa)
# co.isa <- sum.e(c.isa)
c.silvia <- rma(pheno.silvia)
# co.silvia <- sum.e(c.silvia)
# save(co.silvia, co.isa, file = "collapsed.micro.RData")
load("collapsed.micro.RData", verbose = TRUE)

# Merge the data of each batch into a single matrix
co.silvia.df <- as.data.frame(t(co.silvia), row.names = colnames(co.silvia))
co.isa.df <- as.data.frame(t(co.isa), row.names = colnames(co.isa))
merged <- rbind.fill(co.silvia.df, co.isa.df)
rownames(merged) <- c(colnames(co.silvia), colnames(co.isa))
merged.pca <- t(merged)

orig.data <- as.factor(c(rep(1, ncol(co.silvia)), 
                         rep(2, ncol(co.isa))))
pca.graph(data = merged.pca, file = "merged.pca.pdf", 
          col = as.numeric(orig.data),
          outcome = orig.data)
stop("end interesting part")

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


# # Quantifying the counfounding factor
# s <- fast.svd(t(scale(t(exprs(celfiles)), center = TRUE, scale = TRUE)))
# PCA <- s$d ^ 2/sum(s$d ^ 2)
# plot(PCA, type = "b", lwd = 2, las = 1,
#      xlab = "Principal Component", ylab = "Proportion of variance",
#      main = "Principal components contributions")
# dev.off()
# save(c.gcrma, file = "corrected_exprs.RData")
# dev.off()
setwd(origDir)
load("corrected_exprs.RData", verbose = TRUE)


# Prepare the variables to the right format
disease <- read.csv("clean_variables.csv")
data.complete <- cbind("files" = rownames(pData(c.gcrma)), pData(c.gcrma))
vclin <- merge(data.complete, disease, by.x = "Sample", by.y = "id")
int.Var <- c("Sample", "files", "meld", "maddrey", "lille_corte", "lille", "status_90", "glucose",
             "trigycierides", "ast", "alt", "bili_total", "creatinine", 
             "albumin", "inr", "ggt", "ap", "leucos", "hb_g.dl", "hematocrit",
             "platelets", "tp_seg", "hvpg_corte20", "hvpg", "aki",
             "infection_hospitalization")
vclin <- vclin[, colnames(vclin) %in% int.Var]
disease.r <- apply(vclin, 2, as.numeric)
nam <- c("status_90", "infection_hospitalization", "aki", "hvpg_corte20",
         "hvpg_corte20", "lille_corte")
for (n in nam) {
  disease.r[,n] <- as.factor(vclin[,n])
}
disease <- disease.r[, -c(1, 2)]

exp <- exprs(c.gcrma)
# Subset just the AH samples
data.wgcna <- t(exp[, pData(c.gcrma)$Type == "AH"])
gsg <- goodSamplesGenes(data.wgcna, verbose = 3)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", 
                     paste(names(data.wgcna)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(data.wgcna)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  data.wgcna <- data.wgcna[gsg$goodSamples, gsg$goodGenes]
}

pdf("samples.pdf")
sampleTree <- hclust(dist(data.wgcna), method = "average")
pars <- par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

pdf("dendro_traits.pdf")
# Re-cluster samples
sampleTree2 <- hclust(dist(data.wgcna[rownames(data.wgcna) %in% vclin$files, ]),
                      method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(disease, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(disease),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(data.wgcna, vclin, file = "InputWGCNA.RData")


