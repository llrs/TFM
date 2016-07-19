# ==============================================================================
#
#  Code chunk 1: Read the saved data from previous steps
#
# ==============================================================================

pdfn <- function(...){
  # Close any device and open a pdfn with the same options
  dev.off()
  pdfnn(...)
}

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
nSamples <- length(ids) 
# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(data.wgcna[samples %in% ids, ], moduleColors)
MEs <- orderMEs(MEs0$eigengenes)
moduleTraitCor <- cor(MEs, disease, use = "p") 
# Calculating the adjusted p-value
moduleTraitPvalue <- p.adjust(corPvalueStudent(moduleTraitCor, nSamples), "fdr")
dim(moduleTraitPvalue) <- dim(moduleTraitCor)
dimnames(moduleTraitPvalue) <- dimnames(moduleTraitCor)
origPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# TODO: Correct the p-values for multiple testing using the cor.test
# TODO: Calculate the statistical power

result <- apply(MEs, 2, function(x){
  apply(disease, 2, function(y){
    cor.test(x, y, method = "pearson")
  })
})
result3 <- sapply(result, function(x){
  sapply(x, function(y){
    y$p.value
  })
})
result3 <- p.adjust(result3, "fdr")
result3 <- t(result3)
dim(result3) <- dim(moduleTraitCor)
dimnames(result3) <- dimnames(moduleTraitCor)

result4 <- sapply(result, function(x){
  sapply(x, function(y){
    y$estimate
  })
})

# ==============================================================================
#
#  Code chunk 3: Display the correlations of modules and variables in a heatmap
#
# ==============================================================================


pdfn(file = "variables_heatmap.pdfn", width = 10, height = 6, onefile = TRUE)
# Will display correlations and their p-values as text
textMatrix =  paste0(signif(moduleTraitCor, 2), "\n(", 
                           signif(moduleTraitPvalue, 1), ")") 
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3)) 

# Coloring taking into account both the correlation value and the p-value
coloring <- sapply(colnames(moduleTraitCor), function(x){
  moduleTraitCor[, x]*max(moduleTraitPvalue[, x])})

# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = coloring, 
               xLabels = colnames(disease), 
               yLabels = names(MEs), 
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = greenWhiteRed(50), 
               # textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.5, 
               zlim = c(-1, 1), 
               main = paste("Module-trait relationships"))
dev.off()

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
      disease[, "infection_hospitalization"], 
      use = "p"))
GSPvalue <- as.data.frame(
  corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)) 

names(geneTraitSignificance) <- paste("GS.", names(disease)) 
names(GSPvalue) <- paste("p.GS.", names(disease)) 


# ==============================================================================
#
#  Code chunk 5: Plots the relationship between GS and MM of a module
#
# ==============================================================================


module <- "blue"
column <- match(module, modNames) 
moduleGenes <- moduleColors == module 

pdfn(file = "MM_GS_blue.pdfn", width = 7, height = 7) 
par(mfrow = c(1, 1)) 
verboseScatterplot(geneModuleMembership[moduleGenes, column], 
                   geneTraitSignificance[moduleGenes, 1], 
                   xlab = paste("Module Membership in", module, "module"), 
                   ylab = "Gene significance for infection_hospitalization", 
                   main = paste("Module membership vs. gene significance\n"), 
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

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
