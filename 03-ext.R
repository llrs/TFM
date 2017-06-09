# Load ####

source("~/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)
# Load the expression and trait data saved in the first part
load(file = "~/Documents/RNA-seq/Whole_Network.RData", verbose = TRUE)

#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
# rfiles <- list.files(pattern = ".RData")
# sapply(rfiles, load, verbose = TRUE)

#load(file = "modules_ME_orig.RData", verbose = TRUE)
#data.wgcna <- data.wgcna[, moduleColors %in% c("grey60", "darkgrey",
#"plum1", "tan")]
# load(file = "modules_ME.RData", verbose = TRUE)
load(file = "../msimilarity_06.RData", verbose = TRUE)
# load(file = "../msimilarity_mean.RData", verbose = TRUE)
moduleColors <- modules
MEs <- moduleEigengenes(data.wgcna, moduleColors)
MEs <- MEs$eigengenes
# ME var ####
# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
nSamples <- nrow(data.wgcna)

# To include responders
phenoData <- read.csv("../../../data/ALD_ramon/NGS_AHsteps.design.csv", row.names = 1)
phenoData$G <- as.character(phenoData$G)
responders <- ifelse(phenoData$G == "AH.Severe.-.non-responders", 1,
                     ifelse(phenoData$G == "AH.Severe.-.responders", 0, NA))
responders <- as.data.frame(responders)
rownames(responders) <- rownames(phenoData)
vclin <- merge(vclin, responders, all.x = TRUE, all.y = FALSE, by = "row.names")
rownames(vclin) <- vclin[, "Row.names"]
vclin <- vclin[, -1]

# Remove if there isnt' any variability
disease.rm <- apply(vclin, 2, function(x){length(unique(x[!is.na(x)]))}) == 1
n <- apply(vclin, 2, function(x){sum(!is.na(x))})
keep <- n != 0
disease <- vclin[, !disease.rm & keep, drop = FALSE]
# Calculate the non empty clinical data of each variable:
# used for labels and P-value calculation!
disease <- cbind(disease, "random" = runif(nrow(disease))) #Adding dummy variable to plot heatmap
n <- apply(disease, 2, function(x){sum(!is.na(x))})
# n <- ncol(disease)
names.disease <- colnames(disease)
n

# Use just the samples with their clinical data
keep.samples <- rownames(data.wgcna) %in% rownames(vclin)
data.wgcna <- data.wgcna[keep.samples, ]
MEs <- MEs[keep.samples, ]

moduleTrait <- corAndPvalue(MEs, disease, use = "p")
moduleTraitCor <- moduleTrait$cor
keep.variables <- apply(moduleTraitCor, 2, function(x){!all(is.na(x))})
moduleTraitCor <- moduleTraitCor[, keep.variables]
# Calculating the adjusted p-value
# moduleTraitPvalue <- p.adjust(corPvalueStudent(moduleTraitCor, nSamples), "fdr")
# dim(moduleTraitPvalue) <- dim(moduleTraitCor)
# dimnames(moduleTraitPvalue) <- dimnames(moduleTraitCor)
moduleTraitPvalue <- moduleTrait$p[, keep.variables]

# rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)
# moduleTraitPvalue <- orderby(moduleTraitPvalue, rownames(moduleTraitCor),
#                              names.x = TRUE)

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
# Remove the numbers of the names of the colors...
# names(t.colors) <- sapply(names(t.colors), function(x){
#   if (grepl("[0-9]", x)) {
#     nc <- nchar(x)
#     substring(x, 1, nc - 1 )
#   } else {
#     x
#   }
# })
# names(t.colors) <- sapply(names(t.colors), function(x){
#   if (grepl("[0-9]", x)) {
#     nc <- nchar(x)
#     substring(x, 1, nc - 1 )
#   } else {
#     x
#   }
# })

# Number of genes in each group by order of colors.modules
y <- orderby(t.colors, colors.modules)
ylabels <- paste0(orderby(t.colors, colors.modules, names.x = TRUE),
                  " (", y, ")")
yk <- grep("brown |orangered|sienna|skyblue |red|black|turquoise ",
           ylabels, ignore.case = TRUE)
yk.invert <- grep("saddlebrown |lightcyan |darkturquoise |paleturquoise |darkred",
                  ylabels[yk], invert = TRUE)
yk <- yk[yk.invert]
# ylabels <- paste0("Module ", seq_along(y), " (", y, ")")
xk <- grep("Child|MELD|ABIC|Status_90|responders", xlabels)
xk <- xk[-1] # Remove child_clase
# colors_mo <- t(MEs)
# Order the matrix according to the labels
yk2 <- grep("MEbrown|orangered4|sienna3|skyblue3|black|MEturquoise|MEred",
           rownames(colors_mo), ignore.case = TRUE)
xk2 <- grep("Child|MELD|ABIC|Status_90|responders", colnames(textMatrix))
xk2 <- xk2[-1] # Remove child_clase

# Match ones with the others
l <- strsplit(ylabels[yk], " ")
l <- sapply(l, "[", element = 1)
l <- sapply(l, function(x){paste0("ME", x)})
ord <- sapply(l, grep, x = rownames(textMatrix[yk2, ]))
ord2 <- c(2, 1, 3, 7, 5, 4, 6)
# Remove the numbers at the end
ylabels <- sub("(.*)[0-9]+ (.*)","\\1 \\2", ylabels)

# Heatmap ME ####
pdf(file = "heatmap_ME.pdf", width = 10, height = 6,
    onefile = TRUE)
marHeatmap <- c(7, 10, 3, 3)
par(mar = marHeatmap)
labeledHeatmap.multiPage(Matrix = colors_mo,
                         xLabels = xlabels,
                         ySymbols = ylabels,
                         yLabels = names(MEs),
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

labeledHeatmap(Matrix = colors_mo[yk2, xk2][as.numeric(ord), ][ord2, ],
               xLabels = xlabels[xk],
               yLabels = ylabels[yk][ord2],
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix[yk2, xk2][as.numeric(ord), ],
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               12,
               main = "Module Eigengene-trait relationships")
# calculate GS MM ###

modNames <- substring(names(MEs), 3)

geneModule <- corAndPvalue(data.wgcna, MEs, use = "p")
geneModuleMembership <- geneModule$cor
MMPvalue <- geneModule$p

colnames(geneModuleMembership) <- paste0("MM", modNames)
colnames(MMPvalue) <- paste0("p.MM", modNames)

geneTrait <- corAndPvalue(data.wgcna, disease, use = "p")
geneTraitSignificance <- geneTrait$cor
GSPvalue <- geneTrait$p

colnames(geneTraitSignificance) <- paste0("GS.", names.disease)
colnames(GSPvalue) <- paste0("p.GS.", names.disease)

# Convert to matri ####
save(geneModule, geneTrait, file = "heatmaps.RData")
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
par(mar = marHeatmap)
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
par(mar = marHeatmap)
labeledHeatmap.multiPage(Matrix = w.mean,
                         xLabels = xlabels,
                         # xSymbols = colnames(w.mean),
                         yLabels = ylabels,
                         ySymbols = substring(names(MEs), 3),
                         colors = greenWhiteRed(50)[25:50],
                         textMatrix = text.mean,
                         colorLabels = FALSE,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         zlim = c(0, 1),
                         # signed = FALSE,
                         12,
                         addPageNumberToMain = FALSE,
                         maxRowsPerPage = 20,
                         main = "Weighted mean gene significance of modules")
dev.off()

# Plotting modules ####
IM2 <- select.modules(GS.MM.cor, GS.MM.p.value, p.value = 1, threshold = 0)
# Set manually the name of the modules to plot for all the variables
# man.int <- c("MEbrown", "MEfloralwhite", "MEgreen", "MEdarkred",
#              "MEdarkgreen", "MEblack")
# IM2 <- lapply(IM0, function(x){x[x %in% man.int]})
save(IM2, file = "selected_modules.RData")
# fnlist(IM2, "modules_variables.csv")

# e_shs.KRT7    l_shs.KRT7     l_ss.KRT7     e_ss.KRT7     e_us.KRT7     l_us.KRT7
# "darkred"       "brown"       "green" "floralwhite"   "darkgreen"       "black"

# Plot the graphs GS_MM of the interesting modules according to IM2.
 a <- sapply(names(IM2), function(y, d){
   sapply(d[[y]],
          GGMMfun, var = y, MM = geneModuleMembership,
          GS = geneTraitSignificance,
          GSP = GSPvalue, MMP = MMPvalue, moduleColors = moduleColors,
          modNames = modNames, disease = disease)
 }, d = IM2)


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

# Top related genes ####
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
