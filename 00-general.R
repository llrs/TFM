# General options and libraries to be called for each files of the workflow
#
#
# Packages ####
library("GEOquery") # Accessing NCBI GEO datasets
library("affy") # Affymetrix utilities
library("affyPLM") # Affymetrix utilities
library("simpleaffy") # Affymetrix utilities
library("Heatplus") # Heatmaps
library("corpcor") # pcas
library("sva") # combat, batch removing
library("data.table") # data manipulation
library("ggplot2") # Nice plots
library("ggbiplot") # nice PCA
library("corpcor") #
library("topGO") # GO analysis
library("edgeR") # Preprocessing of microarrays
library("limma") # Standard microarrays analysis
library("RColorBrewer") # Colors
library("plyr") # Data manipulation
library("dplyr") # Data manipulation
library("MergeMaid") # Compare two datasets by how they correlate in each
library("annotate") # Annotate things
# library("hgu219.db")
# library("hgu133plus2.db")
# library("Affyhgu133Plus2Expr")
# library("hgu133plus2probe")
# library("hgu133plus2cdf")
library("boot") # Weighted mean, bootstrap...
library("ReactomePA") # Enrichment analysis on reactome
library("clusterProfiler") # Enrichment analysis
library("igraph") # Plotting graphs
library("biomaRt") # Accessing NCBI online data
# library("hgu133plus2.db") # Chip information
library("testthat") # Testing
library("GOstats") # Calculate with GO,
library("graphite")
library("WGCNA") # Weighted gene correlation networks
library("KEGGgraph") # Kegg online accessor
library("KEGG.db") # KEGG database, old release
library("RBGL") # Coloring
library("org.Hs.eg.db") # Translation between id of humans
library("graph") # Graphs representation and handling
library("Rgraphviz") # graphic
# library("lumi") # Analyse illumina chips
library("reshape2")
library("mirbase.db")
library("targetscan.Hs.eg.db")
library("STRINGdb")

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

if (bio.corFnc) {
  source(file.path(code.dir, "bio_cor.R"))
}

# Study's options ####
study <- "miRNA_circulant"

study.dir <- file.path(data.dir, "miRNA")
orig.dir <- setwd(study.dir)
gse.number <- "GSE25609"
path.files <- file.path(study.dir, gse.number)
raw.tar <- paste0(gse.number, "_RAW.tar")
path.raw <- file.path(path.files, raw.tar)
data.out <- file.path(base.dir, study)
dir.create(data.out)
# run.dir <- paste(adj.opt, TOM.opt, sep = "_")
run.dir <- "bicor"
data.files.out <- file.path(data.out, run.dir)
dir.create(data.files.out)

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

# Subset in a list for each module the data
extract <- function(module, clinvar, MM, GS, GSP, MMP, moduleColors) {

  if (is.na(module)) {
    return(NA)
  } else if (substring(module, 1, nchar("ME")) == "ME") {
    module <- substring(module, 3)
  }

  column.mm <- match(paste0("MM", module), colnames(MM))
  column.mmp <- match(paste0("p.MM", module), colnames(MMP))
  moduleGenes <- moduleColors == module
  varc.gs <- match(paste0("GS.", clinvar), colnames(GS))
  varc.gsp <- match(paste0("p.GS.", clinvar), colnames(GSP))

  if (length(varc.gs) > 1 | length(column.mm) > 1) {
    if (length(varc.gs) > 1) {
      warning("Choose between ", paste(colnames(GS)[varc.gs],
                                       collapse = " "))
    } else {
      warning("Choose between ", paste(colnames(MM)[column.mm],
                                      collapse = " "))
    }
  }
  data <- cbind(MM[moduleGenes, column.mm, drop = FALSE],
                GS[moduleGenes, varc.gs, drop = FALSE],
                GSP[moduleGenes, varc.gsp, drop = FALSE],
                MMP[moduleGenes, column.mmp, drop = FALSE])

  colnames(data) <- c("MM", "GS", "GSP", "MMP")
  keep <- apply(data, 1, function(x){!all(is.na(x))})
  data <- unique(data[keep, ])
  if (nrow(data) == 0) {
    return(NA)
  }
  return(data)

}

# Function to calculate the wight
# cor is correlation value between -1 and 1
# p.value1 the p.value of the correlation or a p.value, between 0 and 1
# p.value2 the p.value of other measure, between 0 and 1
# return the values of the weight the higher the better
weight <- function(cor = NULL, p.value1, p.value2 = NULL) {
  if (is.null(cor) & is.null(p.value2)) {
    stop("Which weight you want to calculate?\n",
         "Correlation with p.value or two p.values?")
  }
  if (any(p.value1 > 1) | any(p.value1 < 0)) {
    stop("p.value is not in between 0 and 1")
  }

  if (!is.null(p.value2)) {
    out <- (1 - p.value1) * (1 - p.value2)
  } else {
    out <- cor * (1 - p.value1)
  }
  return(out)
}

# Given data calculates the weight, the weighte correlation and its p.value
w.cor <- function(data) {

  weights <- weight(p.value1 = data$GSP, p.value2 = data$MMP)
  w.cor <- corr(data[, c("MM", "GS")], weights) # requires the boot package
  p.value.w <- corPvalueStudent(w.cor, nrow(data))
  list(weights, w.cor, p.value = p.value.w)
}

plot.GGMM <- function(data) {
  # Weights of the correlation between genes-trait and module-membership
  corr.res <- w.cor(data)
  weight <- corr.res[1]
  w.cor <- corr.res[2]
  p.value.w <- corr.res[3]

  u.cor <- cor(x = data$MM, y = data$GS)
  p.value.u <- corPvalueStudent(u.cor, nrow(data))
  # With ggplot with size for the weights
  ab <- lm(data$GS ~ data$MM)
  plot.g <- ggplot() +
    geom_point(aes(x = data$MM, y = data$GS,
                   size = weight),
               col = ifelse(module %in% c("white", "floralwhite"),
                            "black", module)) + theme_bw() +
    geom_abline(slope = ab$coefficients[2], intercept = ab$coefficients[1]) +
    labs(title = paste0("Module membership vs. gene significance",
                        "\nWeighted cor=", signif(w.cor, digits = 2),
                        ", unweighted cor=", signif(u.cor, digits = 2),
                        " p=",
                        signif(p.value.u, digits = 2), ", ",
                        signif(p.value.w, digits = 2))) +
    xlab(paste("Module Membership in", module, "module")) +
    ylab(paste("Gene significance for", var))
  warning(paste("Plotting", module, "in", var, "."))
  ggsave(filename = name.file("MM_GS", var, module, ".png"),
         plot = plot.g)
}

# Function to explore the module relationship with a trait
GGMMfun <- function(x, var, MM, GS, GSP, MMP, moduleColors, modNames,
                    disease, cor.out = FALSE, p.value = FALSE){

  data <- extract(x, var, MM, GS, GSP, MMP, moduleColors)

  if (is.na(x)) {
    return(NA)
  } else if (substring(x, 1, nchar("ME")) == "ME") {
    module <- substring(x, 3)
  } else {
    module <- x
  }

  # Calculates the weighted mean of genes correlation with the trait
  # wgenecor <- weighted.mean(data[,"GS"], (1 - data[,"GSP"]), na.rm = TRUE)
  # wmmcor <- weighted.mean(data[,"MM"], (1 - data[,"MMP"]), na.rm = TRUE)

  # Weights of the correlation between genes-trait and module-membership
  weight <- (1 - data$GSP) * (1 - data$MMP)
  w.cor <- corr(data[, c("MM", "GS")], weight)
  p.value.w <- corPvalueStudent(w.cor, nrow(data))
  if (cor.out) {
    return(w.cor)
  } else if (p.value) {
    return(p.value.w)
  }
  if (length(data$MM) == 0) {
    return(NA)
  }
  # w.m <- weighted.mean(data[gene, "GS"], w = 1 - data[gene, "GSP"])
  u.cor <- cor(x = data$MM, y = data$GS)
  p.value.u <- corPvalueStudent(u.cor, nrow(data))
  # png(file = paste("MM_GS", var, module, ".png", sep = "_"),
  #     width = 700, height = 700)
  #
  # verboseScatterplot(data[gene, "MM"],
  #                    data[gene, "GS"],
  #                    xlab = paste("Module Membership in", module, "module"),
  #                    ylab = paste("Gene significance for", var),
  #                    main = paste0("Module membership vs. gene ",
  #                                  "significance\nWeighted cor=",
  #                                  signif(w.cor, digits = 2), ", unweighted"),
  #                    abline = 1,
  #                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
  #                    col = ifelse(module %in% c("white", "floralwhite"),
  #                                 "black", module),
  #                    sub = paste("Correlation of genes with trait:",
  #                                "Weighted mean",
  #                                signif(wgenecor, 2), "normal",
  #                                signif(mean(data[, "GS"],
  #                                            na.rm = TRUE), 2)))
  # dev.off()

  # With ggplot with size for the weights
  ab <- lm(data$GS ~ data$MM)
  plot.g <- ggplot() +
    geom_point(aes(x = data$MM, y = data$GS,
                   size = weight),
               col = ifelse(module %in% c("white", "floralwhite"),
                            "black", module)) + theme_bw() +
    geom_abline(slope = ab$coefficients[2], intercept = ab$coefficients[1]) +
    labs(title = paste0("Module membership vs. gene significance",
                        "\nWeighted cor=", signif(w.cor, digits = 2),
                        ", unweighted cor=", signif(u.cor, digits = 2),
                        " p=",
                        signif(p.value.u, digits = 2), ", ",
                        signif(p.value.w, digits = 2))) +
    xlab(paste("Module Membership in", module, "module")) +
    ylab(paste("Gene significance for", var))
  # Point of the weighted mean
  # geom_point(aes(wmmcor, wgenecor),
  #            col = ifelse(module == "red", "green", "red")) +
  # Point of the unweighted mean
  # geom_point(aes(mean(data[gene, "MM"], na.rm = TRUE),
  #                mean(data[gene, "GS"]), na.rm = TRUE),
  #            col = ifelse(module == "green", "orange", "green"))
  warning(paste("Plotting", module, "in", var, "."))
  ggsave(filename = name.file("MM_GS", var, module, ".png"),
         plot = plot.g)
  # dev.off()
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
  # Selecting all the ones that pass the threshold
  if (is.null(ntop)) {
    out <- significant & abs(MTC) >= threshold
    sapply(vclin.names, function(x, y, z) {
      a <- z[y[, x]]
      a[!sapply(a, is.na)]
    }, y = out, z = modules.names, simplify = FALSE)

  } else if (is.numeric(ntop)) {
    # Selecting the top ntop modules based on highest correlation
    sapply(vclin.names, function(a, x, y, z, k) {
      cor.r <- abs(x[y[, a], a])
      # Ordering by the highest correaltion
      a <- names(cor.r)[order(cor.r, decreasing = TRUE)][1:z]
      a[!sapply(a, is.na)]
    }, x = MTC, y = significant, z = ntop, simplify = FALSE)
  } else {
    stop("Please select a number of modules to be selected")
  }
}

# Coloring taking into account both the correlation value and the p-value
# It is a wrapper to the weight function
coloring <- function(MTC, MTP) {
  weight(cor = MTC, p.value1 = MTP)
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

# Print all the elements of a list in columns in the file fil
fnlist <- function(x, fil) {
  nams <- names(x)
  if (file.exists(fil)) {
    warning("Removing existing file ", fil)
    file.remove(fil)
  }
  for (i in seq_along(x)) {
    cat(nams[i], "\t",  x[[i]], "\n",
        file = fil, append = TRUE)
  }
}

select.interesting.modules <- function(MTC, MTP, p.value = 0.07,
                                       threshold = 0.3, ntop = NULL,
                                       var = c("meld", "ggt", "bilirubin"),
                                       var.sign = c(1, 1, -1)) {
  #MTC module trait correlation
  #MTP module trait p.value
  #threshold is the correlation threshold
  # Check that the p.value is minor and that the absolute value of the
  # correlation is >= threshold or that ntop modules are get.
  # TODO: implement this
  # var indicates the clinical variables of interest.
  # var.sign indicate the sign of correlation each module should have with the
  #    clinical variable
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

# Plots GS vs kWithin
connectivity.plot <- function(modules, con, GS, var){
  colorlevels <- unique(modules)
  colorh1 <- modules
  column <- grep(var, colnames(GS))
  if (length(column) >= 2) {
    stop("Variable matches with two columns")
  }
  a <- lapply(colorlevels, function(x){
    if (x == "grey") {
      return(NULL)
    }
    png(name.file("K", var, x, ".png"))
    restrict1 <- (colorh1 == x)
    verboseScatterplot(con$kWithin[restrict1],
         abs(GS[restrict1, column]),
         col = colorh1[restrict1],
         main = x,
         xlab = "Intramodular connectivity (kWithin)",
         ylab = paste("abs(Gene Significance) with", var))
    dev.off()
  })
}

# Plots MM vs kWithin on all modules or return the correlation and p.value
MM_kWithin <- function(MM, con, col, power, cor.out = FALSE, p.value = FALSE) {
  out <- sapply(unique(col), function(x){
    if (x == "grey") {
      return(NULL)
    }
    keepGenes <- col == x
    cor.o <- cor(con$kWithin[keepGenes],
                 abs(MM[keepGenes, paste0("MM", x)] ^ power))
    if (cor.out) {
      return(cor.o)
    } else if (p.value) {
      corPvalueStudent(cor.o, sum(keepGenes))
    }
    png(name.file("MM_Kwithin", x, ".png"))
    verboseScatterplot(con$kWithin[keepGenes],
                       abs(MM[keepGenes, paste0("MM", x)] ^ power),
                       xlab = "Intramodular Connectivity (kWithin)",
                       ylab = paste("Module Membership ^", power),
                       main = paste("Module", x),
                       col = x,
                       abline = TRUE)
    dev.off()

  })
}

# Join with sep, except the last one
# ... the characters to join
# return a character of length 1
name.file <- function(..., sep = "_"){
  arg <- c(...)
  if (length(arg) > 2) {
    paste0(paste0(arg[-length(arg)], collapse = sep), arg[length(arg)])
  } else {
    paste0(arg, collapse = sep)
  }
}

# Helper function to order x by by/y
orderby <- function(x, by, names.x = FALSE) {
  # Both names.x and arrays of 2 dimensions must be provided to order a df
  if (length(dim(x)) == 2 & names.x) {
    out <- x[order(match(rownames(x), by)), ]
  } else {
    out <- x[order(match(names(x), by))]
    if (names.x) {
      out <- names(x)[order(match(names(x), by))]
    }
  }

  return(out)
}

#Remove the 0 of the second position
convert <- function(x){

  ifelse(substring(x, 2, 2) == "0", paste0(substring(x, 1, 1),
                                           substring(x, 3, 3)),
         x)
}

# Check if the power is enough to avoid noise over a normal distribution
checkPower <- function(power, nGenes) {
  out <- 1/sqrt(nGenes) ^ power * nGenes
  if (out >= 0.1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# Extract the data of each hub and represent it for each sample
module.expr <- function(data.wgcna, modules, color) {
  out <- data.wgcna[, modules == color]
  out <- scale(out)
  df <- melt(out)
  ggplot(df, aes(Var1, value, group = factor(Var2))) +
    geom_line(color = color) + xlab("Samples") + ylab("Expression") +
    ggtitle(paste("Expression scaled on module", color)) + theme_bw()
}

# Given miRNA look for targets in the miRBase.db
#
# Uses targetscan to map entrez gene identifers to miRNA
# miRNA a list of character of miRNA
# mature logical, are the miRNA matures? The ids for miRNA in miRBase are if
# mature hsa-miR-xxx if not mature hsa-mir-xxxx
miRNA.target <- function(miRNA, mature = TRUE) {

  if (!mature) {
    myMature <- matureName(get(miRNA, mirbaseMATURE))
  } else {
    myMature <- miRNA
  }
  myMature <- myMature[!is.na(unlist(myMature))]
  myMirFam <- unlist(unique(mget(myMature, targetscan.Hs.egMIRBASE2FAMILY,
                          ifnotfound = NA)))
  if (is.null(myMirFam)) {
    return(NA)
  }
  if (all(is.na(myMirFam))) {
    stop("Couldn't find the miRNA in the dabase")
  }
  myMirFam <- unlist(myMirFam[!is.na(myMirFam)])
  myMirTargets <- unique(unlist(mget(myMirFam,
                                     revmap(targetscan.Hs.egTARGETS))))
  myMirTargets[!is.na(myMirTargets)]
}

# Function to search the targeted genes in the kegg and reactome database
# targets is the name of the targets as the output from miRNA.target
# kegg logical, should the KEGG pathways be used? Caution the kegg package has
# data from 2011
# reactome logical, should the Reactome pathway be used?
target.database <- function(targets, kegg = TRUE, reactome = FALSE) {

  if (kegg) {
    myMirTargetsDB <- unlist(mget(targets, org.Hs.egPATH, ifnotfound = NA))
    myMirTargetsDB <- myMirTargetsDB[!is.na(myMirTargetsDB)]
    myMirTargetsDBNames <- unlist(lapply(myMirTargetsDB,
                                         function(i) {mget(i, KEGGPATHID2NAME)}))
    myMirTargetsDBNames[!is.na(myMirTargetsDBNames)]

  } else if (reactome) {
    myMirTargetsDB <- unlist(mget(targets, reactomeEXTID2PATHID,
                                  ifnotfound = NA))
    myMirTargetsDB <- myMirTargetsDB[!is.na(myMirTargetsDB)]
    myMirTargetsDBNames <- lapply(myMirTargetsDB,
                                  function(i) {mget(i, reactomePATHID2NAME)})
    myMirTargetsDBNames <- unlist(myMirTargetsDBNames)
    myMirTargetsDBNames[!is.na(myMirTargetsDBNames)]
  }
}

# Given the metabolic targets of the groups and the total calculates the
# participation of each group on the total
# group count of times a pathway is targeted by miRNAs of a module/group
# total count of times a pathway is targeted by miRNAs of the univers/chip/study
# absolute logical; should the absolute counts be used or the relative frequency
compare.targets <- function(group, total, absolute = TRUE) {
  if (!is.table(group)) {
    t.group <- table(group)
  } else {
    t.group <- group
  }
  if (!is.table(total)) {
    t.total <- table(total)
  } else {
    t.total <- total
  }
  keep <- names(t.total) %in% names(t.group)
  if (absolute) {
    t.group/t.total[keep]
  } else {
    (t.group/sum(t.group))/(t.total[keep]/sum(t.total))
  }
}
