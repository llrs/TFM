# ==============================================================================
#
#  Code chunk 1: Read the saved data from previous steps
#
# ==============================================================================


library("ggplot2")
# Load the WGCNA package
library("WGCNA")
library("dplyr")
enableWGCNAThreads(6)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
load(file = "InputWGCNA.RData", verbose = TRUE)
load(file = "shared_genes.RData", verbose = TRUE)
#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
load(file = "TNF_AH-network-unsig.RData", verbose = TRUE)
library("boot")
# library("hgu133plus2.db")


# ==============================================================================
#
#  Code chunk 2: Correlate eigenvalues of modules with Variables
#
# ==============================================================================


# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
nSamples <- nrow(vclin)
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


pdf(file = "variables_heatmap_sh.pdf", width = 10, height = 6, onefile = TRUE)
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
y <- table(moduleColors)
colors <- substring(names(MEs), 3)
ylabels <- paste0("ME", names(y[match(names(y), colors)]),
                  " (", y, ")")

# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = coloring,
               xLabels = paste0(colnames(disease), " (", n, ")"),
               yLabels = ylabels,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               12,
               addPageNumberToMain = FALSE,
               main = "Module-trait relationships")
dev.off()
save(moduleTraitCor, moduleTraitPvalue, file = "Module_info_sh.RData")

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

names(geneTraitSignificance) <- paste0("GS.", colnames(disease))
names(GSPvalue) <- paste0("p.GS.", colnames(disease))


# ==============================================================================
#
#  Code chunk 5: Plots the relationship between GS and MM of a module
#
# ==============================================================================
select.modules <- function(MTC, MTP, p.value = 0.07, threshold = 0.3, ntop = NULL) {
  #MTC module trait correlation
  #MTP module trait p.value
  #threshold is the correlation threshold
  # Check that the p.value is minor and that the absolute value of the
  # correlation is >= threshold or that ntop modules are get.
  significant <- MTP <= p.value
  vclin.names <- colnames(MTC)
  modules.names <- rownames(MTC)
  if (is.null(ntop)) {
     out <- significant & abs(MTC) >= threshold
     sapply(vclin.names, function(x, y, z) {
       a <- z[y[, x]]
       a[!sapply(a, is.na)]
     }, y = out, z = modules.names)
  } else {
    sapply(vclin.names, function(a, x, y, z, k) {
      cor.r <- abs(x[y[, a], a])
      a <- names(cor.r)[rank(cor.r) <= z]
      a[!sapply(a, is.na)]
    }, x = MTC, y = significant, z = ntop)
  }
}

IM <- select.modules(moduleTraitCor, moduleTraitPvalue,
                     p.value = 0.05, ntop = 3)

GGMMfun <- function(x, var, MM, GS, GSP, MMP, moduleColors, modNames,
                    disease){
  # Function to explore the module relationship with a trait
  module <- x
  if (is.na(module)) {
    return(NA)
  } else if (substring(module, 1, nchar("ME")) == "ME") {
    module <- substring(module, 3)
  }
  column <- match(module, modNames)

  moduleGenes <- moduleColors == module
  varc <- match(var, colnames(disease))

  data <- cbind("MM" = MM[moduleGenes, column],
                "GS" = GS[moduleGenes, varc],
                "GSP" = GSP[moduleGenes, varc],
                "MMP" = MMP[moduleGenes, column])

  if (ncol(data) == 2 | nrow(data) == 0) {
    return(NA)
  }
  # Calculates the weighted mean of genes correlation with the trait
  wgenecor <- weighted.mean(data[,"GS"], (1 - data[,"GSP"]), na.rm = TRUE)
  wmmcor <- weighted.mean(data[,"MM"], (1 - data[,"MMP"]), na.rm = TRUE)

  # Weights of the correlation between genes-trait and module-membership
  w <- (1 - data[,"GSP"]) * (1 - data[,"MMP"])

  # Taking into account if there are empty values
  gene <- !as.logical(apply(data[,c("MM", "GS")], 1,
                            function(x){sum(is.na(x))}))
  w.cor <- corr(data[gene, c("MM", "GS")], w[gene])
  if (length(data[gene, "MM"]) == 0) {
    return(NA)
  }
  u.cor <- cor(x = data[gene, "MM"], y = data[gene, "GS"])
  png(file = paste("MM_GS", var, module, ".png", sep = "_"),
       width = 700, height = 700)

  verboseScatterplot(data[gene, "MM"], data[gene, "GS"],
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", var),
       main = paste0("Module membership vs. gene significance\nWeighted cor=",
                                   signif(w.cor, digits = 2), ", unweighted"),
                     abline = 1,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                     col = ifelse(module %in% c("white", "floralwhite"),
                                  "black", module),
       sub = paste("Correlation of genes with trait: Weighted mean",
                   signif(wgenecor, 2), "normal",
                   signif(mean(data[, "GS"], na.rm = TRUE), 2)))
  dev.off()


  png(file = paste("MM_GS_ggplot", var, module, ".png", sep = "_"),
      width = 700, height = 700)
  # With ggplot with size for the weights
  ab <- lm(data[gene, "GS"] ~ data[gene, "MM"])
  plot.g <- ggplot() +
    geom_point(aes(x = data[gene, "MM"], y = data[gene, "GS"],
                   size = w[gene]),
               col = ifelse(module %in% c("white", "floralwhite"),
                            "black", module)) + theme_bw() +
  geom_abline(slope = ab$coefficients[2], intercept = ab$coefficients[1]) +
    labs(title = paste0("Module membership vs. gene significance",
                        "\nWeighted cor=", signif(w.cor, digits = 2),
                        ", unweighted cor=", signif(u.cor, digits = 2),
                        " p=",
                        signif(corPvalueStudent(u.cor, nrow(data[gene, ])),
                        digits = 2))) +
         xlab(paste("Module Membership in", module, "module")) +
         ylab(paste("Gene significance for", var))
  # Point of the weighted mean
    # geom_point(aes(wmmcor, wgenecor),
    #            col = ifelse(module == "red", "green", "red")) +
  # Point of the unweighted mean
    # geom_point(aes(mean(data[gene, "MM"], na.rm = TRUE),
    #                mean(data[gene, "GS"]), na.rm = TRUE),
    #            col = ifelse(module == "green", "orange", "green"))

  ggsave(filename = paste("MM_GS", var, module, "ggplot.png", sep = "_"),
         plot = plot.g)
  dev.off()
  }

# Explore for all the variables of trait the selected modules
f.results <- "shared_signed_unsigned"
dir.create(file.path(f.results))
orig <- setwd(file.path(f.results))

a <- sapply(names(IM), function(y, d){
  sapply(d[[y]],
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease)
}, d = IM)
setwd(orig)
# ==============================================================================
#
#  Code chunk 8: Annotate the probes with a gene name
#
# ==============================================================================

# annots <- select(hgu133plus2.db, keys = rownames(exprs),
#                  columns = c("GO", "SYMBOL", "GENENAME", "ENTREZID"),
#                  keytype = "PROBEID")
# save(annots, file = "annots_study.RData")
# load(file = "annots_study.RData", verbose = TRUE)
# dim(annots)
# names(annots)
# probes <- colnames(data.wgcna)
# probes2annot <- match(probes, annots$PROBEID)
# # The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# Should return 0.


# ==============================================================================
#
#  Code chunk 9: Store the results of the WGCNA
#
# ==============================================================================
stop("Done")

# Create a dataframe with gene symbol, modul, MM and MMPvalue
geneInfo <- data.frame(probes = probes,
                        moduleColor = moduleColors)
geneInfo <- merge(geneInfo,
                   unique(annots[,c("PROBEID", "SYMBOL")]),
                   by.x = "probes", by.y = "PROBEID",
                   all.x = TRUE, all.y = FALSE)
# Append all the data of MM and MMP
geneInfo0 <- merge(geneInfo, geneModuleMembership, by.x = "probes",
                   by.y = 0, all = TRUE)
geneInfo0 <- merge(geneInfo0, MMPvalue, by.x = "probes",
                   by.y = 0, all = TRUE)

# Extract them for the specific probes
geneInfo1 <- lapply(unique(geneInfo$moduleColor), function(x, gI){
    mm.positions <- grep(paste0("MM", x), colnames(gI))
    info.re <- gI[gI$moduleColor == x, c(1, mm.positions)]
    MM <- data.frame("probes" = info.re[, 1], "MM" = info.re[, 2])
    MMP <- data.frame("probes" = info.re[, 1], "MMP.value" = info.re[, 3])
    return(list(unique(MM), unique(MMP)))
  }, gI = geneInfo0)

# Concatenate them
MM <- data.frame()
MMP <- data.frame()
for (i in 1:length(geneInfo1)) {
  MM <- rbind(MM, geneInfo1[[i]][[1]])
  MMP <- rbind(MMP, geneInfo1[[i]][[2]])
}
# Add them to the original information
geneInfo2 <- merge(geneInfo, MM, by = "probes", all.x = TRUE)
geneInfo2 <- merge(geneInfo2, MMP, by = "probes", all.x = TRUE)
write.csv(geneInfo2, "genes_modules.csv", row.names = FALSE, na = "")

# Reading genes currently looked up in the laboratory with other experiments
int.genes <- read.csv("genes_int.csv")
int.genes.modules <- geneInfo2[geneInfo2$SYMBOL %in% int.genes$Genes_human, ]
matrx <- table(int.genes.modules$moduleColor, int.genes.modules$SYMBOL)
matrx <- matrx[order(table(int.genes.modules$moduleColor), decreasing = TRUE), ]

genes <- sapply(rownames(matrx), function(x, a){
  paste(colnames(a)[a[x, ] != 0], collapse = ", ")
}, a = matrx)
write.csv(as.data.frame(genes), file = "int_genes_module.csv")

# Foreach module create a table in a file with genes, GS GS-P.values
geneInfo1 <- lapply(unique(geneInfo$moduleColor),
                    function(x, gI, GS, GSP, d, ...){
  gT <- gI[gI$moduleColor == x, ]
  gT <- merge(gT, GS, by.x = "probes", by.y = 0, all.x = TRUE, all.y = FALSE)
  gT <- merge(gT, GSP, by.x = "probes", by.y = 0, all.x = TRUE, all.y = FALSE)
  ord <- sapply(colnames(d), function(x){grep(x, colnames(gT))})
  ord <- c(1, 2, 3, unique(unlist(ord)))
  gT <- gT[, ord]
  write.csv(gT, paste(x, "trait.csv", sep = "_"), row.names = FALSE, na = "")
}, gI = geneInfo, GS = geneTraitSignificance, GSP = GSPvalue, d = disease)
