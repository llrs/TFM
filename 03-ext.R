# Load ####

source("~/Documents/TFM/00-general.R", echo = TRUE)
warning("Done")
setwd(data.files.out)
# Load the expression and trait data saved in the first part
load(file = "Input.RData", verbose = TRUE)

#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
# rfiles <- list.files(pattern = ".RData")
# sapply(rfiles, load, verbose = TRUE)

load(file = "modules_ME_orig.RData", verbose = TRUE)
data.wgcna <- data.wgcna[, moduleColors %in% c("grey60", "darkgrey",
"plum1", "tan")]
load(file = "modules_ME.RData", verbose = TRUE)
# MEs <- MEs$eigengenes
# ME var ####
# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
nSamples <- nrow(data.wgcna)

# Remove if there isnt' any variability
disease.rm <- apply(vclin, 2, function(x){length(unique(x[!is.na(x)]))}) == 1
n <- apply(vclin, 2, function(x){sum(!is.na(x))})
keep <- n != 0
disease <- vclin[, !disease.rm & keep]
# Calculate the non empty clinical data of each variable:
# used for labels and P-value calculation!
n <- apply(disease, 2, function(x){sum(!is.na(x))})
names.disease <- colnames(disease)

# Use just the samples with their clinical data
keep.samples <- rownames(data.wgcna) %in% rownames(vclin)
data.wgcna <- data.wgcna[keep.samples, ]
MEs <- MEs[keep.samples, ]

moduleTraitCor <- cor(MEs, disease, use = "p")

keep.variables <- apply(moduleTraitCor, 2, function(x){!all(is.na(x))})
moduleTraitCor <- moduleTraitCor[, keep.variables]
# Calculating the adjusted p-value
# moduleTraitPvalue <- p.adjust(corPvalueStudent(moduleTraitCor, nSamples), "fdr")
# dim(moduleTraitPvalue) <- dim(moduleTraitCor)
# dimnames(moduleTraitPvalue) <- dimnames(moduleTraitCor)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,
                                      # Filling as many rows as moduleTrait cor
                                      # with the right amount of samples of
                                      # disease used:
                                      t(replicate(nrow(moduleTraitCor), n)))

rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)
moduleTraitPvalue <- orderby(moduleTraitPvalue, rownames(moduleTraitCor),
                             names.x = TRUE)

# Modules variables ####

# Will display correlations and their p-values as text
textMatrix <- paste0(sprintf("%.2f", moduleTraitCor), "\n(",
                     sprintf("%.2e", moduleTraitPvalue), ")")
dim(textMatrix) <- dim(moduleTraitCor)
dimnames(textMatrix) <- dimnames(moduleTraitCor)

# Colors of the
colors_mo <- coloring(moduleTraitCor, moduleTraitPvalue)

# Calculate the number of samples used for the correlation
# X labels
xlabels <- paste0(names.disease, " (", n, ")")

# Y labels
t.colors <- table(moduleColors)
colors.modules <- substring(names(MEs), 3)
# Number of genes in each group by order of colors.modules
y <- orderby(t.colors, colors.modules)
ylabels <- paste0("ME", orderby(t.colors, colors.modules, names.x = TRUE),
                  " (", y, ")")
# Heatmap ME ####
pdf(file = "heatmap_ME.pdf", width = 10, height = 6,
    onefile = TRUE)
par(mar = c(7, 8.5, 3, 3))
labeledHeatmap.multiPage(Matrix = colors_mo,
               xLabels = xlabels,
               yLabels = ylabels,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               12,
               addPageNumberToMain = FALSE,
               main = "Module Eigengene-trait relationships")
dev.off()


# calculate GS MM ###

modNames <- substring(names(MEs), 3)

geneModuleMembership <- cor(data.wgcna, MEs, use = "p")
MMPvalue <- corPvalueStudent(geneModuleMembership,
                            # filling the right n
                            t(replicate(nrow(geneModuleMembership), y)))

colnames(geneModuleMembership) <- paste0("MM", modNames)
colnames(MMPvalue) <- paste0("p.MM", modNames)

geneTraitSignificance <- cor(data.wgcna, disease, use = "p")

GSPvalue <- corPvalueStudent(geneTraitSignificance,
                             t(replicate(nrow(geneTraitSignificance), n)))

colnames(geneTraitSignificance) <- paste0("GS.", names.disease)
colnames(GSPvalue) <- paste0("p.GS.", names.disease)

# Convert to data.frame ####
geneTraitSignificance <- as.matrix(geneTraitSignificance)
geneModuleMembership <- as.matrix(geneModuleMembership)
GSPvalue <- as.matrix(GSPvalue)
MMPvalue <- as.matrix(MMPvalue)
save(geneTraitSignificance,  GSPvalue,
     moduleTraitCor, moduleTraitPvalue,
     geneModuleMembership, MMPvalue, file = "heatmaps.RData")
# Heatmap GS MM ####

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
# calculates the correlation of each module between the GS and the MM
GS.MM.cor <- sapply(names(IM0), function(y, d){
  sapply(d[[y]],
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease, cor.out = TRUE)
}, d = IM0)

# Calculates the p.value of each correlation between GS and MM
GS.MM.p.value <- sapply(names(IM0), function(y, d){
  sapply(d[[y]],
         GGMMfun, var = y, MM = geneModuleMembership,
         GS = geneTraitSignificance,
         GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
         modNames = modNames, disease = disease, p.value = TRUE)
}, d = IM0)

rownames(GS.MM.cor) <- substring(rownames(GS.MM.cor), 3)
GS.MM.cor <- orderby(GS.MM.cor, colors.modules, names.x = TRUE)
rownames(GS.MM.p.value) <- substring(rownames(GS.MM.p.value), 3)
GS.MM.p.value <- orderby(GS.MM.p.value, colors.modules, names.x = TRUE)

# Will display correlations and their p-values as text
textMatrix <- paste0(sprintf("%.2f", GS.MM.cor), "\n(",
                     sprintf("%.2e", GS.MM.p.value), ")")
dim(textMatrix) <- dim(GS.MM.cor)
dimnames(textMatrix) <- dimnames(GS.MM.cor)

colors_mo <- coloring(GS.MM.cor, GS.MM.p.value)

# Calculate the number of samples used for the correlation
ylabels <- paste0("(", orderby(t.colors, rownames(GS.MM.cor)), ") ",
                  orderby(t.colors, rownames(GS.MM.cor), names.x = TRUE))

save(colors_mo, xlabels, ylabels, geneTraitSignificance, geneModuleMembership,
     MEs, file = "heatmap_GS_MM.RData")
pdf("heatmap_GS_MM.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap.multiPage(Matrix = colors_mo,
                         xLabels = xlabels,
                         yLabels = ylabels,
                         ySymbols = substring(names(MEs), 3),
                         colors = greenWhiteRed(50),
                         textMatrix = textMatrix,
                         colorLabels = FALSE,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         zlim = c(-1, 1),
                         12,
                         addPageNumberToMain = FALSE,
                         main = "ModuleMembership-GeneSignificance relationships")
dev.off()

# Heatmap mean ####

w.mean <- sapply(names(IM0), function(y, d) {
  sapply(d[[y]], function(x, var){
    data <- extract(x, var, geneModuleMembership,
                    geneTraitSignificance, GSPvalue,
                    MMPvalue, moduleColors)
    weighted.mean(abs(data[, "GS"]), w = (1 - data[, "GSP"]))
  }, var = y)
}, d = IM0)
text.mean <- paste(sprintf("%.2f", w.mean))
dim(text.mean) <- dim(w.mean)
save(w.mean, text.mean, file = "heatmap_GS_mean.RData")
pdf("heatmap_GS_mean.pdf")
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap.multiPage(Matrix = w.mean,
                         xLabels = colnames(w.mean),
                         # xSymbols = colnames(w.mean),
                         yLabels = rownames(w.mean),
                         # ySymbols = rownames(w.mean),
                         colors = greenWhiteRed(50)[25:50],
                         textMatrix = text.mean,
                         colorLabels = FALSE,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         zlim = c(0, 1),
                         signed = FALSE,
                         12,
                         addPageNumberToMain = FALSE,
                         maxRowsPerPage = 20,
                         main = "Weighted mean gene significance of modules")
dev.off()

# Plotting modules ####
# IM2 <- select.modules(GS.MM.cor, GS.MM.p.value, p.value = 0.05, ntop = 3)
# Set manually the name of the modules to plot for all the variables
# man.int <- c("MEthistle1", "MEfloralwhite", "MEpink", "MEgreen", "MEbrown4",
#              "MElightpink4", "MEbisque4", "MEpaleturquoise", "MEdarkslateblue",
#              "MEthistle2", "MEblue")
# IM2 <- lapply(IM0, function(x){x[x %in% man.int]})
save(IM0, file = "selected_modules.RData")
# fnlist(IM2, "modules_variables.csv")

# Plot the graphs GS_MM of the interesting modules according to IM2.
# a <- sapply(names(IM2), function(y, d){
#   sapply(d[[y]],
#          GGMMfun, var = y, MM = geneModuleMembership,
#          GS = geneTraitSignificance,
#          GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
#          modNames = modNames, disease = disease)
# }, d = IM2)


# GS connectivity ####
# load(file = "kIM.RData", verbose = TRUE)
# load(file = "sft.RData", verbose = TRUE)


# Explore the connectivity of all modules for a variables

# MM vs kWithin
# MM_kWithin(geneModuleMembership, connect, moduleColors,
#            power = sft$powerEstimate)

# GS vs kWithin
# connectivity.plot(moduleColors, connect,
#                   geneTraitSignificance, "meld")
# connectivity.plot(moduleColors, connect,
#                   geneTraitSignificance, "ggt")

# Automatic Screening ####
# # Automatic screening with weighted? screening not know how it works.
# autoScreen <- apply(disease, 2, automaticNetworkScreening,
#                     datExpr = data.wgcna,
#                     datME = MEs,
#                     minimumSampleSize = 4,
#                     power = sft$powerEstimate)
#
# save(autoScreen, file = "autoScreen.RData")

# Tops related genes ####
#
#  Code chunk 6: Explore the genes top related to each clinical variable
#

# genes.interes <- select.genes(geneTraitSignificance, GSPvalue,
#                               p.value = 0.05, ntop = 100,
#                               addMEy = FALSE)
# fnlist(genes.interes, "significant_genes_variables.csv")

# Top related genes per module ####

# genes.modules <- function(disease, GTS, modules, threshold = 0.3) {
#   genes <- GTS[, match(disease, colnames(GTS))]
#   genes.i <- abs(genes) >= threshold
#   colnames(genes)[genes.i]
# }

# Store results ####

# Create a dataframe with gene symbol, modul, MM and MMPvalue
geneInfo <- data.frame(genes = colnames(data.wgcna),
                        moduleColor = moduleColors)
# geneInfo <- merge(geneInfo,
#                    unique(annots[,c("PROBEID", "SYMBOL")]),
#                    by.x = "genes", by.y = "PROBEID",
#                    all.x = TRUE, all.y = FALSE)
#
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
int.genes <- read.csv(file.path(data.dir, "genes_int.csv"))
int.genes.modules <- geneInfo2[geneInfo2$genes %in% int.genes$Genes_human, ]
matrx <- table(int.genes.modules$moduleColor, int.genes.modules$genes)
matrx <- matrx[order(table(int.genes.modules$moduleColor), decreasing = TRUE), ]

genes <- sapply(rownames(matrx), function(x, a){
  paste(colnames(a)[a[x, ] != 0], collapse = ", ")
}, a = matrx)
if (length(genes) == 0) {
  warning("Genes under study were not found. Maybe it is a miRNA study?")
} else {
  write.csv(as.data.frame(genes), file = "int_genes_module.csv")
}

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
