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
pngn <- function(...){
  # Close any device and open a pdfn with the same options
  if (length(dev.list()) > 1) {
    dev.off()
  }
  png(...)
}
library("ggplot2")

# Load the WGCNA package
library(WGCNA)
enableWGCNAThreads(6)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE) 
# Load the expression and trait data saved in the first part
load(file = "InputWGCNA.RData", verbose = TRUE) 
#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
load(file = "TNF_AH-network-auto.RData", verbose = TRUE) 
library("boot")

# ==============================================================================
#
#  Code chunk 2: Correlate eigenvalues of modules with Variables
#
# ==============================================================================


# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
# We don't need to recalculate as we store it previously
# # Recalculate MEs with color labels
# MEs0 <- moduleEigengenes(data.wgcna, moduleColors)
# save(MEs0, file = "ME.RData")
# # load("ME.RData")
# MEs <- orderMEs(MEs0$eigengenes)

disease.r <- apply(vclin, 2, as.numeric)
nam <- c("status_90", "infection_hospitalization", "aki", "hvpg_corte20",
         "hvpg_corte20", "lille_corte")
for (n in nam) {
  a <- as.factor(vclin[,n])
  levels(a)[levels(a) == ""] <- NA
  disease.r[, n] <- a
}
disease <- disease.r[, -c(1, 2)]

nSamples <- ncol(disease)
keepSamples <- rownames(data.wgcna) %in% vclin$files
moduleTraitCor <- cor(MEs[keepSamples, ], disease, 
                      use = "p") 

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
# Calculate the number of samples used for the correlation
n <- apply(disease, 2, function(x){sum(!is.na(x))})

# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = coloring,
               xLabels = paste0(colnames(disease), " (", n, ")"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               12,
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

geneModuleMembership <- as.data.frame(cor(data.wgcna[keepSamples, ], 
                                         MEs[keepSamples, ], use = "p")) 

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples)) 

names(geneModuleMembership) <- paste0("MM", modNames) 
names(MMPvalue) <- paste0("p.MM", modNames) 

geneTraitSignificance <- as.data.frame(
  cor(data.wgcna[keepSamples, ], 
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
warning("Manually made!")
IM <- list( meld = c("skyblue", "darkolivegreen", "midnightblue"),
            maddrey = c("skyblue", "darkolivegreen", "steelblue"),
            # lille_corte = c("darkturquoise", "royalblue", 
            #                 "darkgreen", "darkgrey"),
            lille = c("midnightblue", "salmon", "bisque4", "darkmagenta", 
                      "darkgreen", "lightcyan"),
            status_90 = c("salmon", "darkgreen"),
            glucose = c("brown", "lightsteelblue1"),
            ast = c("green", "white"),
            alt = c("green", "white", "floralwhite"),
            bili_total = c("magenta", "paleturquoise"),
            creatinine = c("plum1", "darkorange"),
            albumin = c("darkgreen", "steelblue"),
            inr = c("darkolivegreen", "plum1", "skyblue"),
            ggt = c("greenyellow", "skyblue3"),
            ap = c("bisque4", "skyblue3", "black"),
            leucos = c("darkmagenta", "red"),
            hb_g.dl = c("darkolivegreen", "darkorange"),
            hematocrit = c("darkorange", "steelblue", "darkolivegreen"),
            platelets = c("green", "white"),
            tp_seg = c("darkolivegreen", "grey"),
            hvpg_corte20 = c("orangered4", "skyblue3", "bisque4"),
            hvpg = c("sienna3", "orange"),
            aki = c("magenta", "paleturquoise"),
            infection_hospitalization = c("navajowhite2", "bisque4"))

GGMMfun <- function(x, var, MM, GS, GSP, MMP, moduleColors, modNames, 
                    disease){
  module <- x
  column <- match(module, modNames) 
  moduleGenes <- moduleColors == module 
  varc <- match(var, colnames(disease))
  
  data <- cbind("MM" = MM[moduleGenes, column],
                "GS" = GS[moduleGenes, varc],
                "GSP" = GSP[moduleGenes, varc],
                "MMP" = MMP[moduleGenes, column])
  
  # Weights of the correlation
  w <- (1 - data[,"GSP"]) * (1 - data[,"MMP"])
  # Taking into account if there are empty values
  svar <- apply(data[,c("MM", "GS")], 1, function(x){sum(is.na(x))})
  w.cor <- corr(data[!as.logical(svar), c("MM", "GS")], w[!as.logical(svar)])
  pngn(file = paste("MM_GS", var, module, ".png", sep = "_"), 
       width = 700, height = 700)
  verboseScatterplot(data[!as.logical(svar), "MM"], 
                     data[!as.logical(svar), "GS"], 
                     xlab = paste("Module Membership in", module, "module"), 
                     ylab = paste("Gene significance for", var), 
       main = paste0("Module membership vs. gene significance\nWeighted cor=", 
                                   signif(w.cor, digits = 2), ", unweighted"),
                     abline = 1, 
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                     col = ifelse(module %in% c("white", "floralwhite"),
                                  "black", module))
}

# GGMMfun("salmon", "status_90", geneModuleMembership, geneTraitSignificance,
#         GSPvalue, MMPvalue, moduleColors, modNames, disease)

a <- sapply(names(IM), function(y, d){
  sapply(d[[y]], 
         # function(x, var){print(paste("module", x, "var", var))},
         # var = y)
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease)
}, d = IM)

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
