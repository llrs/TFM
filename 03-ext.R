# ==============================================================================
#
#  Code chunk 1: Read the saved data from previous steps
#
# ==============================================================================

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)
# Load the expression and trait data saved in the first part
load(file = "Input.RData", verbose = TRUE)

#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
load(file = "modules_ME.RData", verbose = TRUE)
# rfiles <- list.files(pattern = ".RData")
# sapply(rfiles, load, verbose = TRUE)

# ==============================================================================
#
#  Code chunk 2: Correlate eigenvalues of modules with Variables
#
# ==============================================================================


# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
nSamples <- nrow(vclin)

disease.rm <- apply(vclin, 2, function(x){length(unique(x[!is.na(x)]))}) == 1
vclin <- vclin[, !disease.rm]

keepSamples <- rownames(data.wgcna) %in% vclin$files # Samples
disease <- vclin[vclin$files %in% rownames(data.wgcna), 3:ncol(vclin)]
names.disease <- colnames(disease)
names.samples <- vclin$Samples[keepSamples]
if (sum(keepSamples) < 3) {
  disease <- vclin[vclin$Sample %in% rownames(data.wgcna), 3:ncol(vclin)]
  names.disease <- colnames(disease)
  keepSamples <- rownames(data.wgcna) %in% vclin$Sample
  names.samples <- vclin$Samples[keepSamples]
} else if (sum(keepSamples) == 0) {
  stop("Subset correctly the samples with clinical data")
}

if (!all(rownames(data.wgcna[keepSamples]) == names.samples)) {
  stop("Order of samples in clinical variable and expression is not the same!")
}

moduleTraitCor <- cor(MEs[keepSamples, ],
                      disease,
                      use = "p")

keep.variables <- apply(moduleTraitCor, 2, function(x){!all(is.na(x))})
moduleTraitCor <- moduleTraitCor[, keep.variables]
# Calculating the adjusted p-value
# moduleTraitPvalue <- p.adjust(corPvalueStudent(moduleTraitCor, nSamples), "fdr")
# dim(moduleTraitPvalue) <- dim(moduleTraitCor)
# dimnames(moduleTraitPvalue) <- dimnames(moduleTraitCor)

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


# TODO #########################################################################
#
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

# 3 ============================================================================
#
#  Code chunk 3: Display the correlations of modules and variables in a heatmap
#
# ==============================================================================

# Will display correlations and their p-values as text
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)

colors_mo <- coloring(moduleTraitCor, moduleTraitPvalue)
# Calculate the number of samples used for the correlation
n <- apply(disease, 2, function(x){sum(!is.na(x))})
t.colors <- table(moduleColors)
colors <- substring(names(MEs), 3)
ylabels <- paste0("ME", orderby(t.colors, colors, names.x = TRUE),
                  " (", orderby(t.colors, colors), ")")
pdf(file = "heatmap_ME.pdf", width = 10, height = 6,
    onefile = TRUE)
par(mar = c(7, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap.multiPage(Matrix = colors_mo,
               xLabels = paste0(names.disease, " (", n, ")"),
               yLabels = ylabels,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               12,
               addPageNumberToMain = FALSE,
               main = "Module Eigengene-trait relationships")
dev.off()
save(moduleTraitCor, moduleTraitPvalue, file = "Module_info.RData")

# ==============================================================================
#
#  Code chunk 4: Calculates the correlation between gene significance of a
#                variable and the module membership
#
# ==============================================================================

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

names(geneTraitSignificance) <- paste0("GS.", names.disease)
names(GSPvalue) <- paste0("p.GS.", names.disease)


# ==============================================================================
#
#  Code chunk 5: Study the relationship between GS and MM of all modules
#
# ==============================================================================

# IM <- select.modules(moduleTraitCor, moduleTraitPvalue,
#                                  p.value = 0.05)
# IM
# if (length(IM) == 0) {
#   stop("Not significant modules")
# }

IM0 <- select.modules(moduleTraitCor, moduleTraitPvalue,
                            p.value = 1, threshold = 0)
# IM0
# Explore for all the variables of trait the selected modules
GS.MM.cor <- sapply(names(IM0), function(y, d){
  sapply(d[[y]],
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease, cor.out = TRUE)
}, d = IM0)

GS.MM.p.value <- sapply(names(IM0), function(y, d){
  sapply(d[[y]],
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease, p.value = TRUE)
}, d = IM0)

# Plot for all the variables of trait the selected modules
# a <- sapply(names(IM0), function(y, d){
#   sapply(d[[y]],
#          GGMMfun, var = y, MM = geneModuleMembership,
#          GS = geneTraitSignificance,
#          GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
#          modNames = modNames, disease = disease)
# }, d = IM0)

# # Plot the graphs of the interesting modules according to IM.
# a <- sapply(names(IM), function(y, d){
#   sapply(d[[y]],
#          GGMMfun, var = y, MM = geneModuleMembership,
#          GS = geneTraitSignificance,
#          GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
#          modNames = modNames, disease = disease)
# }, d = IM)

GS.MM.cor <- orderby(GS.MM.cor, colors, names.x = TRUE)
GS.MM.p.value <- orderby(GS.MM.p.value, colors, names.x = TRUE)


# Will display correlations and their p-values as text
textMatrix <- paste0(signif(GS.MM.cor, digits = 2), "\n(",
                     signif(GS.MM.p.value, digits = 2), ")")
dim(textMatrix) <- dim(GS.MM.cor)

colors_mo <- coloring(GS.MM.cor, GS.MM.p.value)

# Calculate the number of samples used for the correlation
ylabels <- paste0("(", orderby(t.colors, rownames(GS.MM.cor)), ") ",
                  orderby(t.colors, rownames(GS.MM.cor), names.x = TRUE))

pdf("heatmap_GS_MM.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap.multiPage(Matrix = colors_mo,
                         xLabels = colnames(GS.MM.cor),
                         yLabels = ylabels,
                         colors = greenWhiteRed(50),
                         textMatrix = textMatrix,
                         colorLabels = FALSE,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         12,
                         addPageNumberToMain = FALSE,
                         main = "ModuleMembership-GeneSignificance relationships")
dev.off()

IM2 <- select.modules(GS.MM.cor, GS.MM.p.value, p.value = 0.05)
IM2
# ==============================================================================
#
#  Code chunk 5b: Plots the relationship between GS and connectivity
#
# ==============================================================================
load(file = "kIM.RData", verbose = TRUE)
load(file = "sft.RData", verbose = TRUE)


# Explore the connectivity of all modules for a variables

# MM vs kWithin
# MM_kWithin(geneModuleMembership, connect, moduleColors,
#            power = sft$powerEstimate)

# GS vs kWithin
# connectivity.plot(moduleColors, connect,
#                   geneTraitSignificance, "meld")
# connectivity.plot(moduleColors, connect,
#                   geneTraitSignificance, "ggt")


# Furhter Screening ####
#
autoScreen <- apply(disease, 2, automaticNetworkScreening,
                    datExpr = data.wgcna,
                    datME = MEs,
                    minimumSampleSize = 4,
                    power = sft$powerEstimate)


# ==============================================================================
#
#  Code chunk 6: Explore the genes top related to each clinical variable
#
# ==============================================================================

genes.interes <- select.genes(geneTraitSignificance, GSPvalue,
                              p.value = 0.05, ntop = 100)
fnlist(genes.interes, "significant_genes_variables.csv")

# ==============================================================================
#
#  Code chunk 7: Explore the top genes related to each module in a clinical var
#
# ==============================================================================
# genes.modules <- function(disease, GTS, modules, threshold = 0.3) {
#   genes <- GTS[, match(disease, colnames(GTS))]
#   genes.i <- abs(genes) >= threshold
#   colnames(genes)[genes.i]
# }


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

# Create a dataframe with gene symbol, modul, MM and MMPvalue
geneInfo <- data.frame(genes = colnames(data.wgcna),
                        moduleColor = moduleColors)
# geneInfo <- merge(geneInfo,
#                    unique(annots[,c("PROBEID", "SYMBOL")]),
#                    by.x = "genes", by.y = "PROBEID",
#                    all.x = TRUE, all.y = FALSE)
# Append all the data of MM and MMP
geneInfo0 <- merge(geneInfo, geneModuleMembership, by.x = "genes",
                   by.y = 0, all = TRUE)
geneInfo0 <- merge(geneInfo0, MMPvalue, by.x = "genes",
                   by.y = 0, all = TRUE)

# Extract them for the specific genes
geneInfo1 <- lapply(unique(geneInfo$moduleColor), function(x, gI){
    mm.positions <- grep(paste0("MM", x), colnames(gI))
    info.re <- gI[gI$moduleColor == x, c(1, mm.positions)]
    MM <- data.frame("genes" = info.re[, 1], "MM" = info.re[, 2])
    MMP <- data.frame("genes" = info.re[, 1], "MMP.value" = info.re[, 3])
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
geneInfo2 <- merge(geneInfo, MM, by = "genes", all.x = TRUE)
geneInfo2 <- merge(geneInfo2, MMP, by = "genes", all.x = TRUE)
write.csv(geneInfo2, "genes_modules.csv", row.names = FALSE, na = "")

# Reading genes currently looked up in the laboratory with other experiments
int.genes <- read.csv(file.path(study.dir, "genes_int.csv"))
int.genes.modules <- geneInfo2[geneInfo2$genes %in% int.genes$Genes_human, ]
matrx <- table(int.genes.modules$moduleColor, int.genes.modules$genes)
matrx <- matrx[order(table(int.genes.modules$moduleColor), decreasing = TRUE), ]

genes <- sapply(rownames(matrx), function(x, a){
  paste(colnames(a)[a[x, ] != 0], collapse = ", ")
}, a = matrx)
write.csv(as.data.frame(genes), file = "int_genes_module.csv")

# Foreach module create a table in a file with genes, GS GS-P.values
geneInfo1 <- lapply(unique(geneInfo$moduleColor),
                    function(x, gI, GS, GSP, d, ...){
  gT <- gI[as.character(gI$moduleColor) == x, ]
  gT <- merge(gT, GS, by.x = "genes", by.y = 0, all.x = TRUE, all.y = FALSE)
  gT <- merge(gT, GSP, by.x = "genes", by.y = 0, all.x = TRUE, all.y = FALSE)
  ord <- sapply(colnames(d), function(x){grep(x, colnames(gT))})
  ord <- c(1, 2, unique(unlist(ord)))
  gT <- gT[, ord]
  keep.Trait <- apply(gT, 2, function(x){all(is.na(x))})
  gT <- gT[, !keep.Trait]
  write.csv(gT, name.file(x, "trait.csv"), row.names = FALSE, na = "")
}, gI = geneInfo, GS = geneTraitSignificance, GSP = GSPvalue, d = disease)
