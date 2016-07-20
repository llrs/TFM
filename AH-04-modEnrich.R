#  Analyse the modules with globaltest o GOstats, SPIA, hopach
# First steps similar to the DE analysis

lnames <- load(file = "TNF_AH-network-auto.RData") 
lnames
lnames <- load(file = "Input.RData") 
lnames
lnames <- load(file = "ME.RData")
lnames
lnames <- load(file = "Module_info.RData")
lnames

library("globaltest")
library("limma")
library("SPIA")
library("GOstats")
library("topGO")
library("ggbio")
library("affy")
library("RColorBrewer")
library("gcrma")
library("sva")
library("svd")

moduleName <- "lightcoral"
exprs <- t(data.wgcna[, net$colors  ==  moduleName])
colnames(exprs) <- samples

# Models of the disease
mod <- model.matrix(~., data = disease)
mod0 <- model.matrix(~type, data = celfiles)

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