# ==============================================================================
#
#  Code chunk 1: Read the saved data from previous steps
#
# ==============================================================================

pdfn <- function(...){
  # Close any device and open a pdfn with the same options
  if (length(dev.list()) > 1) {
    dev.off()
  }
  pdf(...)
}
library("ggplot2")

# Load the WGCNA package
library(WGCNA)
enableWGCNAThreads(6)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE) 
# Load the expression and trait data saved in the first part
lnames = load(file = "Input.RData") 
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "TNF_AH-network-auto.RData") 
lnames


# ==============================================================================
#
#  Code chunk 2: Correlate eigenvalues of modules with Variables
#
# ==============================================================================


# Define numbers of genes and samples
nGene <- ncol(data.wgcna) 
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(data.wgcna[samples %in% ids, ], moduleColors)
save(MEs0, file = "ME.RData")
# load("ME.RData")
MEs <- orderMEs(MEs0$eigengenes)
disease <- disease[disease$samples %in% ids, -1]
disease[, "type"] <- as.numeric(disease[, "type"]) - 1
# disease <- disease[, 2]
# disease <- as.matrix(as.numeric(disease) - 1)

nSamples <- ncol(disease)
moduleTraitCor <- cor(MEs, disease, use = "p") 

# Calculating the adjusted p-value
# moduleTraitPvalue <- p.adjust(corPvalueStudent(moduleTraitCor, nSamples), "fdr")
# dim(moduleTraitPvalue) <- dim(moduleTraitCor)
# dimnames(moduleTraitPvalue) <- dimnames(moduleTraitCor)

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# TODO: Correct the p-values for multiple testing using the cor.test
# TODO: Calculate the statistical power
# 
# result <- apply(MEs, 2, function(x){
#   apply(disease, 2, function(y){
#     cor.test(x, y, method = "pearson")
#   })
# })
# 
# extract <- function(core, value = c("statistic", "parameter", "p.value", 
#                                     "estimate", "null.value", 
#                                      "alternative", "method", "data.name")) {
#   # Extract the value of a correlation in a list of lists.
#   a <- sapply(core, function(x){
#     sapply(x, function(y){
#       y$value
#     })
#   })
#   t(a)
# }
# 
# p.val.cor <- extract(result, "p.value")
# adj.p.val.cor <- p.adjust(unlist(p.val.cor), "fdr")
# dim(adj.p.val.cor) <- dim(moduleTraitCor)
# dimnames(adj.p.val.cor) <- dimnames(moduleTraitCor)
# 
# result4 <- extract(result, "estimate")

# ==============================================================================
#
#  Code chunk 3: Display the correlations of modules and variables in a heatmap
#
# ==============================================================================


pdfn(file = "variables_heatmap.pdf", width = 10, height = 6, onefile = TRUE)
# Will display correlations and their p-values as text
textMatrix =  paste0(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 2), ")")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Coloring taking into account both the correlation value and the p-value
coloring <- sapply(colnames(moduleTraitCor), function(x){
  moduleTraitCor[, x]/(1 + moduleTraitPvalue[, x])})
coloring <- sapply(colnames(coloring), function(x){
  y <- coloring[,x]
2*(y - min(y))/(max(y) - min(y)) - 1
})

# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = coloring,
               xLabels = colnames(disease),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               addPageNumberToMain = FALSE,
               main = "Module-trait relationships")

save(moduleTraitCor, moduleTraitPvalue, file = "Module_info.RData")

# ==============================================================================
#
#  Code chunk 4: Calculates the correlation between gene significance of a 
#                variable and the module membership
#
# ==============================================================================


# Define variable weight containing the weight column of datTrait
#weight = as.data.frame(datTraits$weight_g) 
#names(weight) = "weight"
# names (colors) of the modules

modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(data.wgcna[samples %in% ids, ], 
                                         MEs, use = "p")) 

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples)) 

names(geneModuleMembership) <- paste0("MM", modNames) 
names(MMPvalue) <- paste0("p.MM", modNames) 

geneTraitSignificance <- as.data.frame(
  cor(data.wgcna[samples %in% ids, ], 
      disease, 
      use = "p"))
GSPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) 

names(geneTraitSignificance) <- paste("GS.", colnames(disease)) 
names(GSPvalue) <- paste("p.GS.", colnames(disease)) 


# ==============================================================================
#
#  Code chunk 5: Plots the relationship between GS and MM of a module
#
# ==============================================================================

module <- "lightcoral"
column <- match(module, modNames) 
moduleGenes <- moduleColors == module 
var <- "leucos"
varc <- match(var, colnames(disease))
# pdfn(file = "MM_GS_blue.pdf", width = 7, height = 7) 
# par(mfrow = c(1, 1)) 
# 
data <- cbind("MM" = geneModuleMembership[moduleGenes, column],
              "GS" = geneTraitSignificance[moduleGenes, varc],
              "GSP" = GSPvalue[moduleGenes, varc],
              "MMP" = MMPvalue[moduleGenes, column])
# head(data)

# TODO: Plot the values taking into account the p-value
# TODO: Make a correlation of those values taking into account their p-values
# g <- ggplot(as.data.frame(data), aes(MM, GS))
# g + geom_point(aes(colour = GSP, size = MMP)) + theme_bw()
# 
# cor(data[,"MM"]/(data[,"MMP"]/sum(data[,"MMP"])), 
#     data[,"GS"]/(data[,"GSP"]/sum(data[,"GSP"])), 
#     use = 'pairwise.complete.obs')

# Weighted correlation
library("boot")

w <- (1 - data[,"GSP"]) * (1 - data[,"MMP"])
# Taking into account if there are empty values
svar <- apply(data[,c("MM", "GS")], 1, function(x){sum(is.na(x))})
w.cor <- corr(data[!as.logical(svar), c("MM", "GS")], w[!as.logical(svar)])
pdfn(file = paste("MM_GS", module, var, ".pdf", sep = "_"), 
     width = 7, height = 7)
verboseScatterplot(geneModuleMembership[moduleGenes, column], 
 geneTraitSignificance[moduleGenes, varc], 
 xlab = paste("Module Membership in", module, "module"), 
 ylab = paste("Gene significance for", var), 
 main = paste0("Module membership vs. gene significance\nWeighted cor=", 
              signif(w.cor, digits = 2), ", unweighted"),
 abline = 1, 
 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# ==============================================================================
#
#  Code chunk 8: Annotate the probes with a gene name
#
# ==============================================================================

#
# annot = read.csv(file = "GeneAnnotation.csv") 
# dim(annot)
# names(annot)
# probes = names(data.wgcna)
# probes2annot = match(probes, annot$substanceBXH)
# # The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# # Should return 0.


# ==============================================================================
#
#  Code chunk 9
#
# ==============================================================================


# Create the starting data frame
# geneInfo0 <- data.frame(substanceBXH = probes, 
#                       geneSymbol = annot$gene_symbol[probes2annot], 
#                       LocusLinkID = annot$LocusLinkID[probes2annot], 
#                       moduleColor = moduleColors, 
#                       geneTraitSignificance, 
#                       GSPvalue)
# # Order modules by their significance for weight
# modOrder <- order(-abs(cor(MEs, weight, use = "p"))) 
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership)) {
#   oldNames = names(geneInfo0)
#   geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
#                          MMPvalue[, modOrder[mod]]) 
#   names(geneInfo0)<- c(oldNames, paste0("MM.", modNames[modOrder[mod]]), 
#                        paste0("p.MM.", modNames[modOrder[mod]]))
# }
# # Order the genes in the geneInfo variable first by module color, 
# # then by geneTraitSignificance
# geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight)) 
# geneInfo <- geneInfo0[geneOrder, ]
#
#
# ==============================================================================
# 
#  Code chunk 10: Save information about each gene in a csv file
# 
# ==============================================================================
#
#
# write.csv(geneInfo, file = "geneInfo.csv")

dev.off()
