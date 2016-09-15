# ==============================================================================
#
#  Code chunk 1: Read the saved data from previous steps
#
# ==============================================================================

source("00-general.R")

# Load the expression and trait data saved in the first part
load(file = "Input.RData", verbose = TRUE)
#The variable lnames contains the names of loaded variables.
# Load network data saved in the second part.
load(file = "miRNA-network.RData", verbose = TRUE)

data.wgcna <- t(exprs[, pheno$group  == "AH"])

# ==============================================================================
#
#  Code chunk 2: Correlate eigenvalues of modules with Variables
#
# ==============================================================================


# Define numbers of genes and samples
nGene <- ncol(data.wgcna)

# Order pheno to match data.wgcna
pheno <- pheno[match(rownames(data.wgcna), pheno$samplename), ]

vclin <- pheno[, !colnames(pheno) %in% c("patientid", "sampleid",
                                            "samplename", "group")]
disease.r <- apply(vclin, 2, as.numeric)
nam <- c("status_90", "infection_hospitalization", "aki", "mir_182_dico",
         "gender", "follow_up_90", "status_180", "status_1year")
for (n in nam) {
  a <- as.factor(vclin[,n])
  levels(a)[levels(a) == ""] <- NA
  disease.r[, n] <- a
}
disease <- disease.r

nSamples <- nrow(disease)
keepSamples <- rownames(data.wgcna) %in% pheno$samplename
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


pdf(file = "variables_heatmap.pdf", width = 10, height = 6, onefile = TRUE)
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
                         colorLabels = FALSE,
                         colors = greenWhiteRed(50),
                         textMatrix = textMatrix,
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         ycolorLabels = TRUE,
                         12,
                         addPageNumberToMain = FALSE,
                         main = "Module-trait relationships")
dev.off()
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
warning("Manually made: cor >= abs(0.3) and p.value <= 0.05, or close in p.value
        ")
IM <- list( "mir_146a" = c("black", "tan", "cyan"), # n = 3...
            "mir_155" = c("tan"),
            "mir_422" = c("red", "black"),
            "mir_182_dico" = c("yellow", "salmon"),
            "mir_21" = c("midnightblue", "blue"),
            "status_90" = c("cyan", "salmon"), # Not significant
            # "gender" = c(), # Not significant
            "age" = c("turquoise", "yellow"), # Almost yellow significancy
            "mir_182" = c( "midnightblue", "yellow"),
            "meld" = c("yellow", "red"), # Almost red significancy
            "maddrey" = c("yellow", "red"),
            # "lille" = c(), # Not significant low correlations
            # "follow_up_90" = c(),  # Not significant low correlations
            # "status_180" = c(),  # Not significant low correlations
            # "status_1year" = c(),  # Not significant low correlations
            # "glucose" = c(),  # Not significant low correlations
            "ast" = c("cyan", "black"),
            "alt" = c("pink"),
            # "creatinine" = c(),  # Not significant low correlations
            "ggt" = c("midnightblue", "magenta"), # Not significant
            "ap" = c("magenta"),
            "hb_g.dl" = c("cyan"),  # Not significant low correlations
            "tp" = c("salmon"),
            "hvpg" = c("salmon", "tan"), # Not significant
            "aki" = c("cyan")
            # "infection_hospitalization" = c() # Not significant
            )


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
#  Code chunk 8: Store the results of the WGCNA
#
# ==============================================================================


# Create a dataframe with gene symbol, modul, MM and MMPvalue
geneInfo <- data.frame(miRNA = colnames(data.wgcna),
                       moduleColor = moduleColors)
# Append all the data of MM and MMP
geneInfo0 <- merge(geneInfo, geneModuleMembership, by.x = "miRNA",
                   by.y = 0, all = TRUE)
geneInfo0 <- merge(geneInfo0, MMPvalue, by.x = "miRNA",
                   by.y = 0, all = TRUE)

# Extract them for the specific probes
geneInfo1 <- lapply(unique(geneInfo$moduleColor), function(x, gI){
  mm.positions <- grep(paste0("MM", x), colnames(gI))
  info.re <- gI[gI$moduleColor == x, c(1, mm.positions)]
  MM <- data.frame("miRNA" = info.re[, 1], "MM" = info.re[, 2])
  MMP <- data.frame("miRNA" = info.re[, 1], "MMP.value" = info.re[, 3])
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
geneInfo2 <- merge(geneInfo, MM, by = "miRNA", all.x = TRUE)
geneInfo2 <- merge(geneInfo2, MMP, by = "miRNA", all.x = TRUE)
write.csv(geneInfo2, "genes_modules.csv", row.names = FALSE, na = "")


# Foreach module create a table in a file with genes, GS GS-P.values
geneInfo1 <- lapply(
  unique(geneInfo$moduleColor),
  function(x, gI, GS, GSP, d, ...){
    gT <- gI[gI$moduleColor == x, ]
    gT <- merge(gT, GS, by.x = "miRNA", by.y = 0, all.x = TRUE, all.y = FALSE)
    gT <- merge(gT, GSP, by.x = "miRNA", by.y = 0, all.x = TRUE, all.y = FALSE)
    ord <- sapply(colnames(d), function(x){grep(x, colnames(gT))})
    ord <- c(1, 2, 3, unique(unlist(ord)))
    gT <- gT[, ord]
    write.csv(gT, paste(x, "trait.csv", sep = "_"), row.names = FALSE, na = "")
  },
  gI = geneInfo, GS = geneTraitSignificance, GSP = GSPvalue, d = disease)
