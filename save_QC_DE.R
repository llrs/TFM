# Analyse data from article

# Load required libraries
library("GEOquery", quietly = T)
library("affy", quietly = T)
library("simpleaffy", quietly = T)
library("edgeR", quietly = T)
library("limma", quietly = T)
library("gcrma", quietly = T)
library("RColorBrewer", quietly = T)
library("affyPLM", quietly = T)
library("Heatplus")
library("affyQCReport")
library("corpcor")
library("sva")
library("annotate")
library("hgu133plus2.db")
library("data.table")
library("ggbiplot")
library("corpcor")
library("WGCNA")

# Prepare some variables
# gse_number <- "GSE28619"
# path_files <- "../Documents/GSE28619_RAW/"
# test_file <- "GSM709348.CEL.gz"
# path_file <- paste0(path_files, test_file)
setwd("~/Documents")
# raw.tar <- paste0(gse_number, "_RAW.tar")
# path_raw <- paste(gse_number, raw.tar, sep="/")
#By exploring the data on GEO I found that the platform is :
# GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# So we need to load the bioconductor library with the information about it
# library("Affyhgu133Plus2Expr")

# Tutorial on https://www.biostars.org/p/53870/
# Based on this tutorial:
# http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
# Download set
# getGEOSuppFiles(gse_number)

#Unpack the CEL files
# untar(path_raw, exdir="data")
# cels <- list.files("data/", pattern = "[gz]")
# sapply(paste("data", cels, sep="/"), gunzip)
# cels  <-  list.files("data/", pattern = "CEL")
setwd("data/")

# Phenodata contains the information of the experiment space separated! not tab separated
celfiles <- read.affy("phenodata.txt")
QCReport(celfiles, file="ExampleQC_data.pdf") 
# affyQAReport(celfiles, repName = "report1")

# Raw images of the microarrays to visualize, artifacts
image(celfiles)

# Normalize data
c.gcrma <- gcrma(celfiles)
pca.c.gcrma <- prcomp(c.gcrma, scale. = TRUE)
plot(pca.c.gcrma)
g <- ggbiplot(pca.c.gcrma, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
              circle = TRUE)
c.rma <- rma(celfiles)
plot(prcomp(c.rma, scale. = TRUE))

#Write RMA-normalized, mapped data to file
# write.table(grma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# set colour palette
cols <- brewer.pal(8, "Set1")

# plot a boxplot of unnormalised and normalized intensity values

png("normalization_boxplot.png", width=1000, height=750)
pars <- par(mfrow=c(1,3), mar = c(4, 5, 3, 2))
boxplot(celfiles, col=cols, title="Unnormalised")
boxplot(c.gcrma, col=cols, title="Gcrma normalization")
boxplot(c.rma, col=cols, titlw = "RMA normalization")
dev.off()

# Plot a density vs log intensity histogram for the unnormalised and normalised data
png("normalization_histogram.png", width=1500, height=750)
par(mfrow=c(1,3), mar = c(4, 5, 2, 2))
hist(celfiles, col=cols, main="without")
hist(c.gcrma, col=cols, main="gcrma") # This seem to be the best normalization algorithm
hist(c.rma, col=cols, main="rma")
dev.off()

# MA plots to see how weel they are
# add option plot.method = "smoothScatter" to get a fancier plot
png("normalization_maplot.png", width=1500, height=750)
par(mfrow=c(1,3), mar = c(4, 5, 2, 2))
MAplot(celfiles, cex = 0.75, ref.title = "Raw")
MAplot(c.gcrma, cex = 0.75, ref.title = "Gcrma")
MAplot(c.rma, cex = 0.75, ref.title = "rma")
dev.off()
par(pars)

# look at the relationships between the samples using heirarchical clustering
eset <- exprs(c.gcrma)
methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
pars <- par(mfrow=c(2,3))
for (meth in methods){
  distance <- dist(t(eset), method = meth)
  cluster <- hclust(distance)
  plot(cluster, labels = phenoData(celfiles)$Target, main = paste(meth, "distance"))
}
par(pars)

# Depending on the type of distance the clustering goes different

# Check for batch effect (based on day performed)
scanDate <- protocolData(c.gcrma)$ScanDate
scanDate <- gsub(" .*", "", scanDate)
scanDate <- as.Date(scanDate, "%m/%d/%Y")
minscan <- min(scanDate)
days <- scanDate - minscan

days[days == 0] <- 1
days[days == 2] <- 2
days[days == 8] <- 3
days[days == 13] <- 4
batch <- as.numeric(days)

# Most of the AH mirocarrays where done on 2 batches, one of them without control, bias?
table(data.frame(Outcome = c.gcrma$Target, Batch = batch))

# Performs the clustering taking into account the date it was performed.                                                                                                                                                                                                                                                                                                                                                                                        
d.celfiles <- as.dist(1 - cor(exprs(celfiles), method = "spearman"))
d.gcrma <- as.dist(1 - cor(exprs(c.gcrma), method = "spearman"))
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- phenoData(celfiles)$Target
outcome <- as.character(phenoData(celfiles)$Target)
names(outcome) <- sampleNames(celfiles)
sampleDendrogram_c <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  ## for every node in the dendrogram if it is a leaf node
  if (is.leaf(x)) {
    attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")]))## color by batch
    attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
  }
  x
}, batch, outcome) ## these are the second and third arguments in the function

# Plot dendogram with the information of the dates.
plot(sampleDendrogram_c, main = "Hierarchical clustering of samples")
legend("bottom", cex=0.75, paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

# Multidimension scaling (PCA) 
cmd <- cmdscale(d.celfiles)
plot(cmd, type = "n", main="PCA without normalization")
text(cmd, outcome, col = batch, cex = 0.9)
legend("top", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)
# Doesn't show any problem between batch
cmd <- cmdscale(d.gcrma)
plot(cmd, type = "n", main="PCA with gcrma normalization")
text(cmd, outcome, col = batch, cex = 0.9)
legend("top", paste("Batch", unique(batch)), fill = unique(batch), inset = 0.01)

# Quantifying the counfounding factor
s <- fast.svd(t(scale(t(exprs(celfiles)), center = TRUE, scale = TRUE)))
PCA <- s$d^2/sum(s$d^2)
plot(PCA, type="b", lwd=2, las=1,
     xlab="Principal Component", ylab="Proportion of variance",
     main="Principal components contributions")

#Surrogate variable analysis
# To slow in the laptop
# sv <- sva(exprs(celfiles), mod, mod0)
# par(mfrow = c(2, 5))
# for (i in 1:sv$n.sv) boxplot(sv$sv[, i] ~ batch, main = sprintf("SV %d", i), xlab = "Batch")

# Models of the disease
mod <- model.matrix(~Target, data = pData(celfiles))
mod0 <- model.matrix(~1, data = pData(celfiles))

pValues <- f.pvalue(exprs(celfiles), mod, mod0)
sum(p.adjust(pValues, method = "BH") < 0.05)
dim(celfiles)

# Using combat
combatexp <- ComBat(exprs(celfiles), batch, mod)
d <- as.dist(1 - cor(combatexp, method = "spearman"))
sampleClustering <- hclust(d)
sampleDendrogram <- as.dendrogram(sampleClustering, hang = 0.1)
names(batch) <- sampleNames(celfiles)
outcome <- as.character(celfiles$Target)
names(outcome) <- sampleNames(celfiles)
sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch, labels) {
  ## for every node in the dendrogram if it is a leaf node
  if (is.leaf(x)) {
    attr(x, "nodePar") <- list(lab.col = as.vector(batch[attr(x, "label")])) ## color by batch
    attr(x, "label") <- as.vector(labels[attr(x, "label")]) ## label by outcome
  }
  x
}, batch, outcome) ## these are the second and third arguments in the function

plot(sampleDendrogram, main = "Hierarchical clustering of samples with ComBat")
legend("topright", paste("Batch", sort(unique(batch))), fill = sort(unique(batch)))

### DE
# Variability between samples
IQRs <- esApply(c.gcrma, 1, IQR)
plot.ecdf(IQRs, xlab = "IQR", main = "Empirical CDF of IQR values")
abline(v = quantile(IQRs, prob = 0.3), col = "red", lwd = 2)

# Fitting model with the corrected data. 
fit <- lmFit(c.gcrma, mod)
fit <- eBayes(fit)
res <- decideTests(fit)
summary(res)

# Finding the DE genes
tt <- topTable(fit, coef = 2, n = Inf)
pars <- par(mfrow = c(2, 2), mar = c(4, 5, 2, 2))
hist(tt$P.Value, xlab = "Raw P-values", main = "")
hist(tt$P.Value, xlab = "Raw P-values", breaks = 1000, main = "")
hist(tt$adj.P.Val, xlab = "Adjusted P-values", main = "")
hist(tt$adj.P.Val, xlab = "Adjusted P-values", breaks = 1000, main = "")

par(mfrow=c(1,1))
# Volcano plot
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     cex.axis=1.2, las=1, cex.lab=1.5, xlab=expression(paste(log[2], " Fold change")),
     ylab=expression(paste(-log[10], " P-value")), main="Volcano plot FC > 5")
threshold <- log2(5)
points(tt[tt$P.Value < 0.001 & tt$logFC >= threshold, "logFC"], -log10(tt[tt$P.Value < 0.001 & tt$logFC >= threshold, "P.Value"]), pch=".", cex=4,
       col="red")
points(tt[tt$P.Value < 0.001 & tt$logFC <= -threshold, "logFC"], -log10(tt[tt$P.Value < 0.001 & tt$logFC <= -threshold, "P.Value"]), pch=".", cex=4,
       col="green")
abline(h=-log10(max(tt[tt$P.Value < 0.001, "P.Value"])), col=grey(0.5), lty=2)
abline(v=threshold,  col=grey(0.5), lty=2)
abline(v=-threshold,  col=grey(0.5), lty=2)

# Change the names of the probes
# keytypes(hgu133plus2.db)
# k <- head(keys(hgu133plus2.db, keytype="ALIAS"))
# select(hgu133plus2.db, keys=k, columns=c("SYMBOL","REFSEQ"), keytype="SYMBOL")
annots <- select(hgu133plus2.db, keys=rownames(tt),
                 columns=c("SYMBOL","GENENAME", "ACCNUM"), keytype="PROBEID") # Outputs a warning
resultTable0 <- merge(tt, annots, by.x = 0, by.y="PROBEID", all.y = FALSE, all.x = T)
resultTable <- resultTable0[!is.na(resultTable0$SYMBOL), ]
resultTable2 <- resultTable[resultTable$P.Value < 0.001 & abs(resultTable$logFC)  >= threshold,]

tt.accnum <- merge(tt, accnums, by.x = 0, by.y="PROBEID", all.y=FALSE, all.x=TRUE)

# Plot heatmap
par(mfrow=c(1,1))
heatmap.plot <- regHeatmap(exprs(c.gcrma[resultTable2$Row.names,]))
plot(heatmap.plot)

# Names of the table to compare values
extrac <- c("M83248", "AW190565", "NM_000089", "NM_005764", "NM_003247", "K01228", "N30339", "BC001388", "NM_001845", "NM_003254", "AU146808", "X05610", "AK026829", "NM_000393")
inflam <- c("NM_004591", "NM_002993", "NM_004221", "NM_000584", "NM_016639", "NM_001511", "NM_001565", "AI817041")
commun <- c("NM_015515", "BG327863", "BC002700", "NM_002276", "NM_000224")
other <- c("NM_020299", "NM_005564", "NM_006398", "NM_016548", "J04152", "NM_002354", "NM_000903", "AL51445", "NM_001442", "NM_005980")
table_l <- c(extrac, inflam, commun, other)

de <- resultTable2[resultTable2$ACCNUM %in% table_l, ]
de$FC <- 2^abs(de$logFC)

# Write the DE genes to a table and as a list for DAVID
write.table(resultTable2$Row.names, "DE_genes2.txt", col.names=F, row.names=F, quote=F)
write.table(unique(resultTable2[,-length(names(resultTable2))]), file="DE_genes.tsv", row.names = F, quote=F, sep="\t")

# Summarizing with collapseRows
annots.sg <- unique(annots[,-c(4)])
resultTable01 <- merge(tt, annots.sg, by.x = 0, by.y="PROBEID", all.y = FALSE, all.x = T)
collapsed <- collapseRows(resultTable2, rowGroup = resultTable2$ACCNUM, rowID = rownames(resultTable2))
dat1Collapsed <- data.frame(collapsed$group2row, collapsed$datETcollapsed)

# Using table.data package to summarize DE genes
simplyDE <- unique(resultTable2[,-length(names(resultTable2))])
DT <- data.table(simplyDE)

# Takes the max of duplicated
for (x in DT$SYMBOL[duplicated(DT$SYMBOL)]){
  DT[SYMBOL==x] <- DT[SYMBOL==x][which.max(abs(logFC))]
}
DE <- unique(DT)

# Store the data
write.table(DE, file="DE_genes_simplified.tsv", row.names = F, quote=F, sep="\t")
write.table(DE$Row.names, "DE_genes_simplified.txt", col.names=F, row.names=F, quote=F)
