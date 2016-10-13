#  Analyse the modules ####

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
library("STRINGdb")
setwd(data.files.out)

topGO <- FALSE
Reactome <- FALSE
Kegg <- FALSE
GSEA <- FALSE
STRING <- TRUE

# Load previously work done ####
load(file = "Input.RData", verbose = TRUE)
load(file = "modules_ME.RData", verbose = TRUE)
load(file = "selected_modules.RData", verbose = TRUE)

keepSamples <- rownames(data.wgcna) %in% rownames(vclin)

# Define numbers of genes and samples
nGene <- ncol(data.wgcna)
nSamples <- nrow(vclin)
disease.rm <- apply(vclin, 2, function(x){length(unique(x[!is.na(x)]))}) == 1

vclin <- vclin[, !disease.rm]

if (!all(rownames(data.wgcna) == rownames(vclin))) {
  stop("Order of samples in clinical variable and expression is not the same!")
}


# Reconvert the data to the "normal" format, of each column a sample.
exprs <- t(data.wgcna)

# Obtain the annotation of the data
# load(file = "annots_study.RData", verbose = TRUE)

# I assume I keep the same order of genes (Which I do, as I don't change them)
m <- length(unique(moduleColors))
genes <- as.factor(moduleColors)
numb.col <- 0:(m - 1)
names(numb.col) <- c("grey", standardColors(m - 1))
match.colors <- sum(unique(moduleColors) %in% c("grey", standardColors(203)))
if (match.colors != length(unique(moduleColors))) {
  stop("Colors not correctly assigned.")
}
# converts the name to the right number
lg <- levels(genes)
for (x in names(numb.col)) {
  levels(genes)[lg == x] <- numb.col[x]
}
genes <- as.numeric(levels(genes))[genes]
names(genes) <- rownames(exprs)

# Grup all genes of the same group
clusters <- sapply(unique(moduleColors), function(x, genes, nc){
  names(genes[genes == nc[x]])
}, genes = genes, nc = numb.col)


# Compare modules ==============================================================
#
#  Code chunk 0: , general overview of the functions
#


#Function to translate from symbols to entrezid
clustersEntrez <- sapply(clusters, function(x){
  # a <- unique(annots[annots$PROBEID %in% x,
  #                    "ENTREZID"])
  a <- select(org.Hs.eg.db, keys = x, columns = "ENTREZID", keytype = "SYMBOL")
  a <- unique(a)
  a[!is.na(a)]
})

# Compare clusters on a single plot ####
# pdf("clusters_.pdf", onefile = TRUE, width = 20, height = 20)
# x.axis <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
# eGO <- compareCluster(clustersEntrez, fun = "enrichGO")
# save(eGO, file = "eGO.RData")
# plot(eGO) + ggtitle("Enrich GO") + x.axis
# cGO <- compareCluster(clustersEntrez, fun = "enrichGO", ont = "CC")
# save(cGO, file = "eGO.RData")
# plot(cGO) + ggtitle("Enrich cc GO") + x.axis
# # gGO <- compareCluster(clustersEntrez, fun = "groupGO")
# # save(gGO, file = "gGO.RData")
# # plot(gGO) + ggtitle("Group GO") + x.axis
# eP <- compareCluster(clustersEntrez, fun = "enrichPathway")
# save(eP, file = "eP.RData")
# plot(eP) + ggtitle("Enrich Pathways") + x.axis
# eK <- compareCluster(clustersEntrez, fun = "enrichKEGG")
# save(eK, file = "eK.RData")
# plot(eK) + ggtitle("Enrich KEGG") + x.axis
# dev.off()

string_db <- STRINGdb$new(version = "10", species = 9606,
                           score_threshold = 0, input_directory = "" )

imodules <- unique(unlist(IM2))

universeGenesEntrez <- unique(AnnotationDbi::keys(org.Hs.eg.db))
universeGenesEntrez <- universeGenesEntrez[!is.na(universeGenesEntrez)]

# Study each module ####
out <- sapply(imodules, function(x) {

  if (substr(x, 0, 2) %in% c("ME", "MM")) {
    moduleName <- substring(x, 3)
  } else {
    moduleName <- x
  }

  message(paste("Analyzing", moduleName, "module!"))
  selFun <- moduleSel(moduleName, numb.col)

  # Preparing the objects with Entrezid for the reactome and kegg analysis
  moduleGenes <- clusters[moduleName][[1]]
  moduleGenesEntrez <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = moduleGenes,
                                     keytype = "SYMBOL",
                                     columns = "ENTREZID"))
  moduleGenesEntrez <- moduleGenesEntrez[!is.na(moduleGenesEntrez)]


  # topGO ======================================================================
  #
  #  Code chunk 1: topGO analysis of each module

  # topGOdata object ####
  if (topGO) {
    GOdata.bp <- new("topGOdata",
                  ontology = "BP",
                  description = paste("Biological process of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = "alias",
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)

    GOdata.mp <- new("topGOdata",
                  ontology = "MP",
                  description = paste("Molecular process of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = "alias",
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)
    GOdata.cc <- new("topGOdata",
                  ontology = "CC",
                  description = paste("Cellular component of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = "alias",
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)

    # save(GOdata, file = "array_BP.RData")
    # load(file = "array_BP.RData", verbose = TRUE)
    # GOdata
    # geneSelectionFun(GOdata) <- selFun
    # GOdata
    # description(GOdata) <- paste("Molecular function of the",
    #                              moduleName, "module.")
    # topGO tests ####
    go.enrich(GOdata.cc, moduleName, "CC")
    go.enrich(GOdata.bp, moduleName, "BP")
    go.enrich(GOdata.mp, moduleName, "MP")
  }
  # Reactome ===================================================================
  #
  #  Code chunk 2: Reactome analysis of the module
  #
  if (Reactome) {
    reactome_enrich <- enrichPathway(gene = moduleGenesEntrez,
                                     universe = universeGenesEntrez,
                                     pvalueCutoff = 0.05, readable = TRUE,
                                     minGSSize = 2)
    if (nrow(summary(reactome_enrich)) != 0) {
      write.csv(summary(reactome_enrich),
                file = paste0("reactome_", moduleName, ".csv"))
      pdf(paste0("reactome_", moduleName, ".pdf"), onefile = TRUE)
      tryCatch({dotplot(reactome_enrich)},
               error = function(e) {
                 message("Couldn't plot the dotplot for reactome")
                 message(e)
               })

      # One can use the fold change to visualize how are the genes expressed
      # with a foldChange = vector
      tryCatch({cnetplot(reactome_enrich, showCategory = 15,
                         categorySize = "geneRatio",
                         layout = layout_nicely)},
               error = function(e) {
                 message("Couldn't plot the cnetplot for reactome")
                 message(e)
               })
      # summary(reactome_enrich)
      # dput(summary(reactome_enrich))
      # Can't have titles
      tryCatch({enrichMap(reactome_enrich, layout = layout_nicely,
                          vertex.label.cex = 1, n = 15)},
               error = function(e) {
                 message("Couldn't map the enrichMap for reactome")
                 message(e)
               })
      dev.off()
    }
  }

  #  Kegg ======================================================================
  #
  #  Code chunk 3: Kegg analysis of each module
  #
  if (Kegg) {
    kegg_enrich <- enrichKEGG(moduleGenesEntrez,
                              universe = universeGenesEntrez,
                              use_internal_data = TRUE,
                              minGSSize = 2)
    if (nrow(summary(kegg_enrich)) != 0) {
      write.csv(summary(kegg_enrich),
                file = paste0("kegg_", moduleName, ".csv"))
      pdf(paste0("kegg_", moduleName, ".pdf"), onefile = TRUE)
      tryCatch({dotplot(kegg_enrich)},
               error = function(e) {
                 message("Couldn't map the dotplot for KEGG")
                 message(e)
               })
      tryCatch({cnetplot(kegg_enrich, showCategory = 15, categorySize = "geneNum",
                         layout = igraph::layout_nicely)},
               error = function(e) {
                 message("Couldn't map the cnetplot for KEGG")
                 message(e)
               })
      tryCatch({enrichMap(kegg_enrich, layout = igraph::layout_nicely,
                          vertex.label.cex = 1, n = 15)},
               error = function(e) {
                 message("Couldn't map the enrichMap for KEGG")
                 message(e)
               })
      dev.off()
    }
  }
  # GSEA =======================================================================
  #
  #  Code chunk 4: GSEA
  #
  if (GSEA) {
    gse <- gsePathway(as.vector(moduleGenesEntrez)[order(moduleGenesEntrez,
                                                         decreasing = TRUE)],
                      nPerm = 1000, pvalueCutoff = 0.2,
                      pAdjustMethod = "BH", verbose = TRUE)

    # General plotting
    if (!is.null(gse)) {
      summary(gse)
      message(paste("Plotting GSE for module", moduleName))
      pdfn(paste0("gsea_", moduleName, ".pdf"), onefile = TRUE)
      enrichMap(gse)
      dev.off()
    }
    # # Individual gene plot
    # gseaplot(gse, geneSetID = moduleGenesEntrez[1])
  }

  # STRING =====================================================================
  if (STRING) {
    if (length(moduleGenes) < 375) {
      string_id <- string_db$map(data.frame(gene = moduleGenes), "gene",
                                 takeFirst = FALSE)
      # For each variable it could paint the GS of each gene
      # payload <- string_db$post_payload(string_id$STRING_id, )
      png(name.file("STRING_module", moduleName, ".png"),
          width = 2000, height = 2000, bg = "transparent")
      string_db$plot_network(string_id$STRING_id, # payload_id = payload,
                             add_link = FALSE)
      dev.off()
      # png(name.file("STRING_enrichment", moduleName, ".png"))
      # tryCatch({
      #   tryCatch({string_db$plot_ppi_enrichment(string_id$STRING_id,
      #                                           title = paste("Interaction enrichment of", moduleName))},
      #            error = function(e) {
      #              dput(string_id$STRING_id)
      #              message("Couldn't map the enrichment of STRING")
      #              # sessionInfo()
      #              # message(e)
      #            })
      # }, error = function(e) {
      #   message(length(string_id$STRING_id))
      #   message("Another error")
      #   # message(e)
      # })
      # dev.off()
    } else {
      warning(paste("Moduel", moduleName, "is bigger!"))
    }
  }
})

