#  Analyse the modules ####

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

compare <- TRUE
topGO <- FALSE
Reactome <- FALSE
Kegg <- FALSE
GSEA <- FALSE
# load("../clusters.RData")
# load(design, )
STRING <- FALSE

# Initial format of input all will be converted to entrez
keytype <- "REFSEQ" # c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS",
# "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
# "GO", "GOALL", "IPI", "MAP", "OMIM", "ONTOLOGY", "ONTOLOGYALL",
# "PATH", "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UCSCKG",
# "UNIGENE", "UNIPROT")
# Entrez for GO
GO.ID <- "entrez" # c("entrez", "genbank", "alias", "ensembl", "symbol",
                  # "genename", "unigene")

# Load previously work done ####
load(file = "Input.RData", verbose = TRUE)
# load(file = "modules_ME.RData", verbose = TRUE)
# load(file = "modules_ME_orig.RData", verbose = TRUE)
# data.wgcna <- data.wgcna[, moduleColors %in% c("grey60", "darkgrey",
#                                                "plum1", "tan")]
load(file = "modules_ME.RData", verbose = TRUE)
# MEs <- MEs$eigengenes
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

# Convert all into Entrezid. # And only keep those
names.genes <- unique(AnnotationDbi::select(org.Hs.eg.db,
                                             keys = rownames(exprs),
                                             keytype = keytype,
                                             columns = "ENTREZID"))
name <- data.frame(keytype = rownames(exprs), mod = genes)
ids <- merge(name, names.genes, by.y = keytype, by.x = "keytype")
genes <- ids$mod
names(genes) <- ids$ENTREZID

# Grup all genes of the same group
clusters <- sapply(unique(moduleColors), function(x, genes, nc){
  ng <- names(genes[genes == nc[x]])
  ng[!is.na(ng)]
}, genes = genes, nc = numb.col)

save(clusters, file = "modules_entrezid.RData")
# Compare modules ####

 if (compare) {
   pdf("clusters_.pdf", onefile = TRUE, width = 20, height = 20)
   x.axis <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
   eGO <- compareCluster(clusters, fun = "enrichGO")
   save(eGO, file = "eGO.RData")
   plot(eGO) + ggtitle("Enrich GO") + x.axis
   cGO <- compareCluster(clusters, fun = "enrichGO", ont = "CC")
   save(cGO, file = "cGO.RData")
   plot(cGO) + ggtitle("Enrich cc GO") + x.axis
   # gGO <- compareCluster(clusters, fun = "groupGO")
   # save(gGO, file = "gGO.RData")
   # plot(gGO) + ggtitle("Group GO") + x.axis
   eP <- compareCluster(clusters, fun = "enrichPathway")
   save(eP, file = "eP.RData")
   plot(eP) + ggtitle("Enrich Pathways") + x.axis
   eK <- compareCluster(clusters, fun = "enrichKEGG")
   save(eK, file = "eK.RData")
   plot(eK) + ggtitle("Enrich KEGG") + x.axis
   dev.off()
 }

string_db <- STRINGdb$new(version = "10", species = 9606,
                           score_threshold = 0, input_directory = "" )

imodules <- unique(unlist(IM2))

 if (Reactome | Kegg) {
   universeGenesEntrez <- unique(AnnotationDbi::keys(org.Hs.eg.db))
   universeGenesEntrez <- universeGenesEntrez[!is.na(universeGenesEntrez)]
 }

if (GSEA) {
  # TODO
  # Requires exprs with controls and other conditions and design
  # to be able to calculate the contrast
  mc <- makeContrasts(, levels = design)
  keep <- ids2indices(gene.sets, rownames(exprs))
  camera(exprs, keep, design, contrast = 2)
  roast(exprs, keep, design, contrast = 2)
  romer(exprs, keep, design, contrast = 2)
}

# Study each module ####
out <- sapply(imodules, function(x) {

  if (substr(x, 0, 2) %in% c("ME", "MM")) {
    moduleName <- substring(x, 3)
  } else {
    moduleName <- x
  }

  message("Analyzing ", moduleName, " module!")
  selFun <- moduleSel(moduleName, numb.col)

  # Preparing the objects with Entrezid for the reactome and kegg analysis
  moduleGenes <- clusters[moduleName][[1]]

  # topGO ######
  if (topGO) {
    tryCatch({GOdata.bp <- new("topGOdata",
                  ontology = "BP",
                  description = paste("Biological process of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = GO.ID,
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)
    go.enrich(GOdata.bp, moduleName, "BP")}, error = function(x){
      message("\nUnable to calculate BP.\n")
      message(x)
    })

    tryCatch({GOdata.mp <- new("topGOdata",
                  ontology = "MP",
                  description = paste("Molecular process of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = GO.ID,
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)
    go.enrich(GOdata.mp, moduleName, "MP")}, error = function(x){
      message("\nUnable to calculate MP.\n")
      message(x)
    })
    tryCatch({GOdata.cc <- new("topGOdata",
                  ontology = "CC",
                  description = paste("Cellular component of the",
                                      moduleName, "module."),
                  allGenes = genes,
                  # annot = annFUN.gene2GO, ## the new annotation function
                  # affyLib = "org.Hs.eg.db",
                  annot = annFUN.org,
                  ID = GO.ID,
                  mapping = "org.Hs.eg",
                  geneSelectionFun = selFun)
    go.enrich(GOdata.cc, moduleName, "CC")}, error = function(x){
      message("\nUnable to calculate CC.\n")
      message(x)
    })
  }
  # Reactome ####
  if (Reactome) {
    reactome_enrich <- enrichPathway(gene = moduleGenes,
                                     universe = universeGenesEntrez,
                                     pvalueCutoff = 0.05, readable = TRUE,
                                     minGSSize = 2)
    if (length(summary(reactome_enrich)) != 0) {
      write.csv(summary(reactome_enrich),
                file = paste0("reactome_", moduleName, ".csv"),
                row.names = FALSE)
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
    } else {
      message("Not enough data in reactome for module ", moduleName)
    }
  }

  #  Kegg ####
  if (Kegg) {
    kegg_enrich <- enrichKEGG(moduleGenes,
                              universe = universeGenesEntrez,
                              use_internal_data = TRUE,
                              minGSSize = 2)
    if (length(summary(kegg_enrich)) != 0) {
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
    } else {
      message("Not enough data in KEGG for module ", moduleName)
    }
  }

  # GSEA ####
  if (GSEA) {
    message("Calculating GSEA in KEGG in module ", moduleName)
    # TODO correct the way moduleGenes are ordered, to a more meaningful
    # Use some of the ideas of R. Castelo
    gse <- gsePathway(as.vector(moduleGenes)[order(moduleGenes,
                                                   decreasing = TRUE)],
                      nPerm = 1000, pvalueCutoff = 0.2,
                      pAdjustMethod = "BH", verbose = TRUE)

    # General plotting
    if (!is.null(gse)) {
      summary(gse)
      message("Plotting GSE for module ", moduleName)
      pdfn(paste0("gsea_", moduleName, ".pdf"), onefile = TRUE)
      enrichMap(gse)
      dev.off()
    }
    # # Individual gene plot
    # gseaplot(gse, geneSetID = moduleGenes[1])
  }

  # STRING ####
  if (STRING) {
    if (length(moduleGenes) < 375) {
      string_id <- string_db$map(data.frame(gene = moduleGenes), "gene",
                                 takeFirst = FALSE)
      # For each variable it could paint the GS of each gene
      # payload <- string_db$post_payload(string_id$STRING_id, )
      png(name.file("STRING_module", moduleName, ".png"),
          width = 2000, height = 2000)
      string_db$plot_network(string_id$STRING_id, # payload_id = payload,
                             add_link = FALSE)
      dev.off()
    } else {
      warning(paste("Moduel", moduleName, "is bigger!"))
    }
  }
})

