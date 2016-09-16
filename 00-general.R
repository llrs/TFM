# General options and libraries to be called for each files of the workflow
#
#
# Packages ####
library("GEOquery")
library("affy")
library("affyPLM")
library("simpleaffy")
library("Heatplus")
library("corpcor")
library("sva")
library("data.table")
library("ggplot2")
library("ggbiplot")
library("corpcor")
library("topGO")
library("WGCNA")
library("edgeR")
library("limma")
library("RColorBrewer")
library("plyr")
library("dplyr")
library("MergeMaid")
library("annotate")
library("hgu219.db")
library("hgu133plus2.db")
library("Affyhgu133Plus2Expr")
library("hgu133plus2probe")
library("hgu133plus2cdf")
library("boot")
library("Rgraphviz")
library("ReactomePA")
library("clusterProfiler")
library("igraph")
library("org.Hs.eg.db")

# Options and configurations ####

enableWGCNAThreads(6) # Speeding up certain calculations with multi-threading.
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Negative correlations are as much valued as positive cor.
adj.opt <- "unsigned"
# Reduce the impact in genes when correlations are both positive and negative
TOM.opt <- "signed"
# Powers to test with
powers <- c(1:30)
base.dir <- "/home/lrevilla/Documents"
data.dir <- file.path(base.dir, "data")
code.dir <- file.path(base.dir, "TFM")
bio.corFnc <- FALSE


# Study's options ####
study <- "comparison"
pheno1 <- "pheno.isa.txt"
pheno2 <- "pheno.silvia.txt"

study.dir <- file.path(data.dir, "hepatitis")
orig.dir <- setwd(study.dir)
gse.number <- "GSE28619"
path.files <- file.path(study.dir, paste0(gse.number, "_RAW"))
raw.tar <- paste0(gse.number, "_RAW.tar")
path.raw <- file.path(data.dir, raw.tar)
data.out <- file.path(base.dir, study)
dir.create(data.out)
run.dir <- paste(adj.opt, TOM.opt, sep = "_")
data.files.out <- file.path(data.out, run.dir)
dir.create(data.files.out)

# Development ####

source(file.path(code.dir, "bio_cor.R"))

# Functions ####

# Normalize data and plots PCA of samples
pca.graph <- function(celfiles=NULL, data=NULL, file, outcome = NULL,
                      col = NULL, ...){
  # Data is the normalized, celfiles are the raw files
  if (is.null(data)) {
    if (!is.null(celfiles)) {
      data <- rma(celfiles)
    } else {
      stop("celfiles or data must be provided")
    }
  }

  if ("ExpressionSet" %in% is(data)) {
    if (is.null(outcome)) {
      outcome <- as.character(phenoData(data)$Type)
    }
    dists <- as.dist(1 - cor(exprs(data), method = "spearman"))
  } else {
    if (is.null(outcome)) {
      outcome <- rep("AH", ncol(data))
    }
    dists <- as.dist(1 - WGCNA::cor(data, method = "spearman",
                                    use = "pairwise.complete.obs"))
  }

  cmd <- cmdscale(dists, eig = TRUE, add = TRUE)
  perc.v <- round(cmd$eig/sum(cmd$eig)*100, 1)
  pdf(file)
  pca.plo <- ggbiplot(prcomp(dists, scale. = TRUE),
                      obs.scale = 1, var.scale = 1, ellipse = TRUE,
                      group = outcome, var.axes = FALSE,
                      # circle = TRUE
  )
  plot(pca.plo)
  plot(cmd$points[, 1], cmd$points[, 2], type = "n", main = "MDS samples",
       xlab = paste0("PC1 (", perc.v[1], "% explained var.)"),
       ylab = paste0("PC2 (", perc.v[2], "% explained var.)"),
       ...)
  text(cmd$points[,1 ], cmd$points[, 2], col = col, cex = 0.9,
       labels = outcome)
  dev.off()
  invisible(data)
}

# Given a expression set transforms it to the gene symbols and averages the expr
sum.e <- function(eset){

  ann <- annotation(eset)
  pkg <- paste0(ann, ".db")
  pkg.ann <- eval(parse(text = pkg))
  # print(class(pkg.ann))
  symbols <- AnnotationDbi::select(pkg.ann, keys(pkg.ann), "SYMBOL")
  expr <- exprs(eset)
  probes <- rownames(expr)
  rowGroup <- symbols[match(probes, symbols$PROBE),"SYMBOL"]
  # the length of rowGroup and probes should be the same
  # even if there is a warning we cannot omit it
  # length(unique(rowGroup)) == length(probes.isa)
  corr.e <- collapseRows(expr, rowGroup = rowGroup, rowID = probes,
                         method = "Average")
  corr.sy <- corr.e$datETcollapsed
  return(corr.sy)
}

# Close any device and open a pdf. Allow the same options as pdf
pdfn <- function(...){

  if (length(dev.list()) > 1) {
    dev.off()
  }
  pdf(...)
}

# Calculate the percentatge above per of data
count.p <- function(data, per){
  sum(data >= per)/length(data)
}

# Function to explore the module relationship with a trait
GGMMfun <- function(x, var, MM, GS, GSP, MMP, moduleColors, modNames,
                    disease){
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

  verboseScatterplot(data[gene, "MM"],
                     data[gene, "GS"],
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", var),
                     main = paste0("Module membership vs. gene ",
                                   "significance\nWeighted cor=",
                                   signif(w.cor, digits = 2), ", unweighted"),
                     abline = 1,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                     col = ifelse(module %in% c("white", "floralwhite"),
                                  "black", module),
                     sub = paste("Correlation of genes with trait:",
                                 "Weighted mean",
                                 signif(wgenecor, 2), "normal",
                                 signif(mean(data[, "GS"],
                                             na.rm = TRUE), 2)))
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

# Select genes above a threshold of correlation or the top ntop
select.genes <- function(GTS, GSP, p.value = 0.05, threshold = 0.3,
                         ntop = NULL) {
  # GTS is the geneTraitSignificance
  # GSP is the p-values
  #
  vclin <- substring(colnames(GTS), 4)
  colnames(GTS) <- vclin
  colnames(GSP) <- vclin
  genes <- rownames(GTS)

  sign <- GSP <= p.value
  if (is.null(ntop)) {
    out <- sign & abs(GTS) >= threshold
    sapply(colnames(GSP), function(x, y, z) {
      a <- z[y[, x]]
      a[!sapply(a, is.na)]
    }, y = out, z = rownames(GTS))
  } else {
    sapply(vclin, function(a, x, y, z, k) {
      cor.r <- abs(x[y[, a], a])
      names(cor.r) <- k[y[, a]]
      b <- names(cor.r)[rank(cor.r) <= z]
      b[!sapply(b, is.na)]
    }, x = GTS, y = sign, z = ntop, k = genes)
  }
}

# Select modules by ME above a threshold of correlation or the top ntop
select.modules <- function(MTC, MTP, p.value = 0.07,
                           threshold = 0.3, ntop = NULL) {
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

# Coloring taking into account both the correlation value and the p-value
coloring <- function(MTC, MTP) {
  colors <- sapply(colnames(MTC), function(x){
    MTC[, x]/(1 + MTP[, x])})
  colors.value <- sapply(colnames(colors), function(x){
    y <- coloring[,x]
    2*(y - min(y))/(max(y) - min(y)) - 1
  })
  colors.value
}

# Function to generate function to select the module
moduleSel <- function(modul, a){
  selFun <- function(genes){
    # Function to select those genees of the same group
    # return(a[x])
    return(genes == a[modul])
  }
  return(selFun)
}

# Change the name of the column
rename.col <- function(df, colname, new.colname) {
  colnames(df)[colnames(df) == colname] <- new.colname
  df
}

# Convert factors using eq level names as factor and eq values as number
fact2num <- function(vec, level, new.level){
  vec <- as.factor(vec)
  levels(vec) <- c(levels(vec), new.level)
  vec[grep(level, vec, ignore.case = TRUE)] <- new.level
  vec <- droplevels(vec)
  vec
}

# If any level which is a number is left, then is converted to NA
level.na <- function(vec){
  vec <- as.factor(vec)
  levels(vec) <- as.numeric(levels(vec))
  vec
}

# returns string w/o leading whitespace
trim.leading <- function(x) {
  sub("^\\s+", "", x)
}

# returns string w/o trailing whitespace
trim.trailing <- function(x) {
  sub("\\s+$", "", x)
}

# returns string w/o leading or trailing whitespace
trim <- function(x) {
  gsub("^\\s+|\\s+$", "", x)
}


