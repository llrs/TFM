# File to store the common code used in the branch
# Packages ####

library("WGCNA")
library("ggplot2")


# Options ####

enableWGCNAThreads(6)
options(stringsAsFactors = FALSE)

# Functions ####

pdfn <- function(...){
  # Close any device and open a pdfn with the same options
  if (length(dev.list()) > 1) {
    dev.off()
  }
  pdf(...)
}

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
  # Calculates the weighted mean of genes correlation with the trait
  wgenecor <- weighted.mean(data[,"GS"], (1 - data[,"GSP"]), na.rm = TRUE)
  wmmcor <- weighted.mean(data[,"MM"], (1 - data[,"MMP"]), na.rm = TRUE)

  # Taking into account if there are empty values
  gene <- !as.logical(apply(data[,c("MM", "GS")], 1,
                            function(x){sum(is.na(x))}))
  w.cor <- corr(data[gene, c("MM", "GS")], w[gene])
  u.cor <- cor(x = data[gene, "MM"], y = data[gene, "GS"])

  png(file = paste("MM_GS", var, module, ".png", sep = "_"),
      width = 700, height = 700)
  verboseScatterplot(data[gene, "MM"], data[gene, "GS"],
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = paste("Gene significance for", var),
                     main = paste0(
                       "Module membership vs. gene significance\nWeighted cor=",
                       signif(w.cor, digits = 2), ", unweighted"),
                     abline = 1,
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                     col = ifelse(module %in% c("white", "floralwhite"),
                                  "black", module),
                     sub = paste(
                       "Correlation of genes with trait: Weighted mean",
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
}

