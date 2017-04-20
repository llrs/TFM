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
library("hgu219.db")
library("Affyhgu133Plus2Expr")
library("hgu133plus2probe")
library("hgu133plus2cdf")
library("boot") # Weighted mean, bootstrap...
library("ReactomePA") # Enrichment analysis on reactome
library("clusterProfiler") # Enrichment analysis
library("igraph") # Plotting graphs
library("biomaRt") # Accessing NCBI online data
library("hgu133plus2.db") # Chip information
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
# library("biomaRt")
library("reactome.db")
library("AnnotationDbi")
library("STRINGdb")
library("foreach")
library("doParallel")
library("BiocParallel")
library("parallel")
library("GSVA")
library("GSVAdata")
library("GSEABase")
library("snowfall")
library("BioCor")
# Options and configurations ####

enableWGCNAThreads(4) # Speeding up certain calculations with multi-threading.
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Negative correlations are as much valued as positive cor.
adj.opt <- "signed hybrid"
# Reduce the impact in genes when correlations are both positive and negative
TOM.opt <- "signed"
# Powers to test with
powers <- c(1:30)
base.dir <- "~/Documents"
data.dir <- file.path(base.dir, "data")
code.dir <- file.path(base.dir, "TFM")
bio.corFnc <- FALSE

# Study's options ####
study <- "RNA-seq"
pheno1 <- "pheno.isa.txt"
pheno2 <- "pheno.silvia.txt"
rd <- "POS_NEG_TOTAL_16SAMPLES.csv"
data.folder <- "ductular_reaction"
study.dir <- file.path(data.dir, data.folder)
orig.dir <- setwd(study.dir)
gse.number <- "GSE28619"
path.files <- file.path(study.dir, paste0(gse.number, "_RAW"))
raw.tar <- paste0(gse.number, "_RAW.tar")
path.raw <- file.path(data.dir, raw.tar)
data.out <- file.path(base.dir, study)
dir.create(data.out)
# subdirectory_opt <- paste(adj.opt, TOM.opt, sep = "_")
subdirectory_opt <- "late"
subdirectory <- "Consensus"
data.files.out <- file.path(data.out, subdirectory)#, subdirectory_opt) #08_01_01
dir.create(data.files.out)

# Functions ####

# Normalize data and plots PCA of samples
pca.graph <- function(celfiles=NULL, data=NULL, file, outcome = NULL,
                      col = NULL, ...) {
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
                      labels = outcome
                      # circle = TRUE
  ) + theme_bw() + ggtitle("PCA samples")
  plot(pca.plo)
  plot(cmd$points[, 1], cmd$points[, 2], type = "n", main = "MDS samples",
       xlab = paste0("PC1 (", perc.v[1], "% explained var.)"),
       ylab = paste0("PC2 (", perc.v[2], "% explained var.)"),
       ...)
  text(cmd$points[,1 ], cmd$points[, 2], col = col, cex = 0.9,
       labels = outcome)
  dev.off()
  invisible(cmd)
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
  if (sum(is.na(c(column.mmp, column.mmp, varc.gs, varc.gsp))) >= 1) {
    stop("Unable to find the apropiate data")
  }
  if (length(varc.gs) > 1 | length(column.mm) > 1) {
    if (length(varc.gs) > 1) {
      warning("Choose between ", paste(colnames(GS)[varc.gs],
                                       collapse = " "))
    } else {
      warning("Choose between ", paste(colnames(MM)[column.mm],
                                      collapse = " "))
    }
  }
  mod_data <- cbind(MM[moduleGenes, column.mm, drop = FALSE],
                GS[moduleGenes, varc.gs, drop = FALSE],
                GSP[moduleGenes, varc.gsp, drop = FALSE],
                MMP[moduleGenes, column.mmp, drop = FALSE])

  colnames(mod_data) <- c("MM", "GS", "GSP", "MMP")
  keep <- apply(mod_data, 1, function(x){!all(is.na(x))})
  mod_data <- unique(mod_data[keep, ])
  if (nrow(mod_data) == 0) {
    warning("Data for ", module, " isn't available")
    return(NA)
  }
  return(mod_data)

}

# Given data calculates the weight, the weighte correlation and its p.value
w.cor <- function(mod_data) {
  weight <- (1 - mod_data[, "GSP"]) * (1 - mod_data[, "MMP"])
  w.cor <- corr(mod_data[, c("MM", "GS")], weight) # requires the boot package
  p.value.w <- corPvalueStudent(w.cor, nrow(mod_data))
  list(weight, w.cor, p.value = p.value.w)
}

plot.GGMM <- function(mod_data, module) {
  # Weights of the correlation between genes-trait and module-membership
  if (sum(is.na(mod_data)) >= 1) {
    warning("data for ", module, "wasn't available")
  }
  corr.res <- w.cor(mod_data)
  weight <- corr.res[1]
  w.cor <- corr.res[2]
  p.value.w <- corr.res[3]

  u.cor <- cor(x = mod_data[, "MM"], y = mod_data[, "GS"])
  p.value.u <- corPvalueStudent(u.cor, nrow(mod_data))
  # With ggplot with size for the weights
  ab <- lm(mod_data[, "GS"] ~ mod_data[, "MM"])
  plot.g <- ggplot() +
    geom_point(aes(x = mod_data[, "MM"], y = mod_data[, "GS"],
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
  message(paste("Plotting", module, "in", var, "."))
  ggsave(filename = name.file("MM_GS", var, module, ".png"),
         plot = plot.g)
}

# Function to explore the module relationship with a trait
GGMMfun <- function(x, var, MM, GS, GSP, MMP, moduleColors, modNames,
                    disease, cor.out = FALSE, p.value = FALSE){

  mod_data <- extract(x, var, MM, GS, GSP, MMP, moduleColors)
  if (sum(is.na(mod_data)) >= 1) {
    return(NA)
  }

  if (is.na(x)) {
    return(NA)
  } else if (substring(x, 1, nchar("ME")) == "ME") {
    module <- substring(x, 3)
  } else {
    module <- x
  }

  # Calculates the weighted mean of genes correlation with the trait
  # wgenecor <- weighted.mean(mod_data[,"GS"], (1 - mod_data[,"GSP"]), na.rm = TRUE)
  # wmmcor <- weighted.mean(mod_data[,"MM"], (1 - mod_data[,"MMP"]), na.rm = TRUE)

  # Weights of the correlation between genes-trait and module-membership
  weight <- (1 - mod_data[, "GSP"]) * (1 - mod_data[, "MMP"])
  w.cor <- corr(mod_data[, c("MM", "GS")], weight)
  p.value.w <- corPvalueStudent(w.cor, nrow(mod_data))
  if (cor.out) {
    return(w.cor)
  } else if (p.value) {
    return(p.value.w)
  }
  if (length(mod_data[, "MM"]) == 0) {
    return(NA)
  }
  # browser()
  # w.m <- weighted.mean(mod_data[gene, "GS"], w = 1 - mod_data[gene, "GSP"])
  u.cor <- cor(x = mod_data[, "MM"], y = mod_data[, "GS"])
  p.value.u <- corPvalueStudent(u.cor, nrow(mod_data))

  # With ggplot with size for the weights
  ab <- lm(mod_data[, "GS"] ~ mod_data[, "MM"])
  plot.g <- ggplot() +
    geom_point(aes(x = mod_data[, "MM"], y = mod_data[, "GS"],
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
# MTC the value of correlation
# MTP the p-value of such correlations
# p.value the threshold of p.values
# threshold the threshold of the correlation
# ntop if present the top ntop number of modules is selected for each variables
# return a list of variables and modules selected.
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
coloring <- function(MTC, MTP) {
  colors <- sapply(colnames(MTC), function(x){
    MTC[, x]*(1 - MTP[, x])})
  colors
}

# Function to generate function to select the module
moduleSel <- function(modul, a){
  b <- a[modul]
  selFun <- function(genes){
    # Function to select those genees of the same group
    # return(a[x])
    # message("Using ", b)
    return(genes == b)
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
    cat(nams[i], "\t",  x[[i]], "\n", file = fil, append = TRUE)
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

#Remove the 0 of the second position of a string
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
# data.wgcna is the data with rows each sample and columns the genes
# modules is the assignment of modules of each gene
# color is the color to plot
# amplitude is the b-a difference of the normalization
# centered logical; it should have the 0 inside
# returns a ggplot object
module.expr <- function(data.wgcna, modules, color, amplitude = 2,
                        centered = FALSE) {
  out <- data.wgcna[, modules == color]
  if (centered) {
    max.x <- amplitude/2
    min.x <- -max.x
  } else {
    max.x <- amplitude
    min.x <- 0
  }
  out <- apply(out, 1, scaling, min.x = min.x, max.x = max.x)
  df <- melt(t(out))
  ggplot(df, aes(Var1, value, group = factor(Var2))) + theme_bw() +
    geom_line(color = color) + xlab("Samples") + ylab("Expression")
}

# x numbers to be scaled
# min.x minimum value of values scaled
# max.x maximum value of values scaled
scaling <- function(x, min.x = -1 , max.x = 1) {
  return((max.x - min.x)*(x - min(x, na.rm = T))/(
    max(x, na.rm = T) - min(x, na.rm = T)) + min.x)
}

# Function which computes the enrichement in GOdata objects and plot them
go.enrich <- function(GOdata, moduleName, ont) {
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(GOdata, algorithm = 'weight01', statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  avgResult <- combineResults(resultFisher, resultKS, resultKS.elim,
                              method = "mean")

  allRes <- GenTable(GOdata, weight01 = resultKS,
                     elim = resultKS.elim, classic = resultFisher,
                     topNodes = 50, numChar = 1000)
  write.csv(allRes, file = name.file("GO_table", ont, moduleName, ".csv"),
            row.names = FALSE)

  pdf(name.file("GO", ont, moduleName, ".pdf"), onefile = TRUE)

  showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 2,
                 useInfo = 'all')
  title(main = "GO analysis using Fisher algorithm")
  showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 2,
                 useInfo = 'all')
  title(main = "GO analysis using Weight01 algorithm")
  showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 2,
                 useInfo = 'all')
  title(main = "GO analysis using KS elim algorithm")
  showSigOfNodes(GOdata, score(avgResult), firstSigNodes = 2, useInfo = 'all')
  title(main = "GO analysis using average")
  dev.off()
}

# Plots in a pdf the object in a module
# path_object the object to be used for plotting
# moduleName name of the modulename
pathway.enrich <- function(path_object, moduleName, prefix = "kegg_") {
    pdf(paste0(prefix, moduleName, ".pdf"), onefile = TRUE)
    tryCatch({dotplot(path_object)},
             error = function(e) {
               message("Couldn't draw the dotplot")
               message(e, "\n")
             })
    tryCatch({cnetplot(path_object, showCategory = 15, categorySize = "geneNum",
                       layout = igraph::layout_nicely)},
             error = function(e) {
               message("Couldn't draw the cnetplot")
               message(e, "\n")
             })
    tryCatch({enrichMap(path_object, layout = igraph::layout_nicely,
                        vertex.label.cex = 1, n = 15)},
             error = function(e) {
               message("Couldn't draw the enrichMap")
               message(e, "\n")
             })
    dev.off()
}

# Given many datasets for WGCNA it create a multiExpr data set
# Works either with traits and expression data
multiSet <- function(...) {
  sets <- list(...)
  multiExpr <- vector(mode = "list", length = length(sets))
  names(multiExpr) <- names(sets)
  for (i in 1:length(sets)) {
    setExpr <- sets[[i]]
    multiExpr[[i]] <- list(data = as.data.frame(setExpr))
  }

  keep.columns <- Reduce(intersect, lapply(multiExpr,
                                           function(x){colnames(x$data)}))
  message("Keeping ", length(keep.columns), " columns of the sets.")
  for (i in 1:length(sets)) {
    setExpr <- sets[[i]]
    multiExpr[[i]]$data <- multiExpr[[i]]$data[, colnames(setExpr) %in% keep.columns]
  }
  out <- checkSets(multiExpr, checkStructure = TRUE)
  if (out$structureOK & out$nGenes == length(keep.columns)) {
    return(multiExpr)
  } else {
    stop("Unable to create a multiSet/multiExpr, with the sets given.")
  }
}

# Given a powerTable finds the minimum power to reach the
# min value of threshold for SFT.R.sq fitting.
multiple.softThreshold <- function(multiPower, min = 0.85){
  unlist(lapply(multiPower, function(x){
    y <- x$data[x$data$SFT.R.sq > min, "Power"]
    return(y[1])}))
}

# Function to advise on the weight to use
# data.wgcna expression data in WGCNA format
# power  power used to build the network
# adj.opt options of the adjacency construction
# bio_mat the previously calculated biological information
# Tom.opt how to build the TOM matrix
weight.bio.cor <- function(data.wgcna, power, adj.opt, bio_mat, TOM.opt){
  g <- seq(0, 1, by = 0.1)
  combin.weights <- expand.grid(g, g, g, g)
  colnames(combin.weights) <- c("Expr",  names(bio_mat))
  combin.weights <- unique(combin.weights[, !is.na(colnames(combin.weights))])
  keep.weights <- apply(combin.weights, 1, function(x)(sum(x) == 1))
  combin.weights <- combin.weights[keep.weights, ]
  # Filter to those whose the expression is at least 50%
  combin.weights <- combin.weights[combin.weights[, "Expr"] >= 0.5, ]

  adj <- adjacency(data.wgcna, type = adj.opt, power = power)
  out <- apply(combin.weights, 1, function(we){
    adj.bio <- cor.all(adj, bio_mat, weights = we)
    TOM <- TOMsimilarity(adj.bio, TOMType = TOM.opt)
    dissTOM <- 1 - TOM
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = 30)
    moduleColors <- labels2colors(dynamicMods$labels)
    # prop.table(table(moduleColors))}
    table(moduleColors)
    }
  )
  # Select column names
  nam <- unique(unlist(sapply(out, names)))
  makeDF <- function(List, Names) { # Set the element to the right column
    m <- t(vapply(List,
                  FUN = function(X) unlist(X)[Names],
                  FUN.VALUE = numeric(length(Names))))
    as.data.frame(m)
  }

  tables <- makeDF(out, nam)
  colnames(tables) <- nam # In some cases it doesn't apply correctly the names
  full <- cbind(combin.weights, tables)
  rownames(full) <- 1:nrow(full)
  return(full)
}


# Plot the combinations of parameters given by weight.bio.cor
# full the weights and size of the resulting modules
# labels.bio_mat the name of the data introduced as bio_mat in weight.bio.cor
weight.plot <- function(full, labels.bio_mat) {
  new.df <- melt(full, id.vars = labels.bio_mat)
  new.df <- new.df[!is.na(new.df$value), ]
  ggplot(new.df, aes(x = value, fill = variable)) +
    facet_wrap(labels.bio_mat,
               labeller = labeller(.rows = label_both, .multi_line = FALSE)) +
    geom_histogram(bins = 15) +
    theme_bw() +
    scale_x_log10() +
    # scale_y_log10() +
    scale_fill_manual(values = unique(as.character(new.df$variable))) +
    xlab("# of genes by module") +
    ylab("# of modules") +
    guides(fill = FALSE)

  ggplot(new.df) +
    geom_violin(aes(x = type, y = value)) +
    theme_bw() +
    scale_y_log10() +
    ylab("# of genes by module") +
    xlab("Method to build the network") +
    geom_point(aes(x = type, y = value, color = variable, size = value)) +
    scale_color_manual(values = unique(as.character(new.df$variable))) +
    theme(legend.position = "none") +
    ggtitle("Distribution and size of the modules in each network")
}

# exprMat The expression matrice
# gctFn gct file name.
# Create a file with the specificities of GSE java application
gsea.write.gct <- function(exprMat, gctFn) {
  nGenes <- nrow(exprMat)
  nConds <- ncol(exprMat)
  write("#1.2", file = gctFn, append = FALSE) # All .gct files require this dummy header line.
  write(paste(nGenes, nConds, sep = "\t"), file = gctFn, append = TRUE)
  write(paste("Name", "Description", paste(colnames(exprMat), collapse = "\t"), sep = "\t"), file = gctFn, append = T)
  # The second column of the .gct file, "Description", is filled out with "na"'s.
  rownames(exprMat) <- paste(rownames(exprMat), "na", sep = "\t") # Append "\tna" to every gene name.
  write.table(exprMat, file = gctFn, append = TRUE, quote = FALSE, sep = "\t",
              na = "", row.names = TRUE, col.names = FALSE)
}

# Write the cls file of the phenotype to be used
# pheno the phenotype of the data
# column the column of the pheno
# file_name the name of the output file
# returns a file
gsea.write.cls <- function(pheno, column, file_name) {
  # (number of samples) (space) (number of classes) (space) 1
  clases <- unique(pheno[, column])
  write(paste(nrow(pheno), length(clases), "1"),
        file = file_name)
  write(paste("#", paste(clases, collapse = " ")), file = file_name, append = TRUE)
  write(paste(pheno[, column], collapse = " "), file = file_name, append = TRUE)
}
