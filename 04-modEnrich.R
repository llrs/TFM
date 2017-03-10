#  Analyse the modules ####

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)
setwd(data.files.out)

compare <- FALSE
topGO <- FALSE
Reactome <- TRUE
Kegg <- TRUE
GSEA <- FALSE
# load("../clusters.RData")
# load(design, )
STRING <- TRUE

# Initial format of input all will be converted to entrez
keytype <- "SYMBOL" # c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS",
# "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME",
# "GO", "GOALL", "IPI", "MAP", "OMIM", "ONTOLOGY", "ONTOLOGYALL",
# "PATH", "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UCSCKG",
# "UNIGENE", "UNIPROT")
# Entrez for GO
GO.ID <- "entrez" # c("entrez", "genbank", "alias", "ensembl", "symbol",
                  # "genename", "unigene")

# Load previously work done ####
load(file = "../../Late_Network.RData", verbose = TRUE)
#load(file = "modules_ME.RData", verbose = TRUE)
#load(file = "modules_ME_orig.RData", verbose = TRUE)
#data.wgcna <- data.wgcna[, moduleColors %in% c("grey60", "darkgrey",
#                                               "plum1", "tan")]
load(file = "modules_ME.RData", verbose = TRUE)
# MEs <- MEs$eigengenes
load(file = "selected_modules.RData", verbose = TRUE)

# Reconvert the data to the "normal" format, of each column a sample.
exprs <- t(data.wgcna)

# Obtain the annotation of the data
# load(file = "annots_study.RData", verbose = TRUE)

# I assume I keep the same order of genes (Which I do, as I don't change them)
m <- length(unique(moduleColors))
genes <- as.factor(moduleColors)
numb.col <- 0:(m - 1)
names(numb.col) <- labels2colors(numb.col)
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
  x.axis <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
  eGO <- compareCluster(clusters, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                        ont = "MF")
  save(eGO, file = "eGO.RData")
  plot.eGO <- dotplot(eGO) + ggtitle("Enrich MF GO") + x.axis
  bpGO <- compareCluster(clusters, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                         ont = "BP")
  save(bpGO, file = "bpGO.RData")
  plot.bpGO <- dotplot(bpGO) + ggtitle("Enrich BP GO") + x.axis
  cGO <- compareCluster(clusters, fun = "enrichGO", ont = "CC",
                        OrgDb = org.Hs.eg.db)
  save(cGO, file = "cGO.RData")
  plot.cGO <- dotplot(cGO) + ggtitle("Enrich CC GO") + x.axis
  # gGO <- compareCluster(clusters, fun = "groupGO")
  # save(gGO, file = "gGO.RData")
  # plot(gGO) + ggtitle("Group GO") + x.axis
  eP <- compareCluster(clusters, fun = "enrichPathway")
  save(eP, file = "eP.RData")
  plot.eP <- dotplot(eP) + ggtitle("Enrich Pathways") + x.axis
  eK <- compareCluster(clusters, fun = "enrichKEGG")
  save(eK, file = "eK.RData")
  plot.eK <- dotplot(eK) + ggtitle("Enrich KEGG") + x.axis
  pdf("clusters_.pdf", onefile = TRUE, width = 20, height = 20)
  plots <- list(plot.eGO, plot.bpGO, plot.eP, plot.eK)
  invisible(lapply(plots, print))
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
message("Analysing modules ", paste(imodules, collapse = ", "))
out <- sapply(imodules, function(x) {

  if (substr(x, 0, 2) %in% c("ME", "MM")) {
    moduleName <- substring(x, 3)
  } else {
    moduleName <- x
  }

  message("Analyzing ", moduleName, " module!")

  # Preparing the objects with Entrezid for the reactome and kegg analysis
  moduleGenes <- clusters[[moduleName]]

  # topGO ######
  if (topGO) {
    # selFun <- moduleSel("grey")
    GOdata.bp <- new("topGOdata",
                     ontology = "BP",
                     description = "Biological process of the genes in the study.",
                     allGenes = genes,
                     # annot = annFUN.gene2GO, ## the new annotation function
                     # affyLib = "org.Hs.eg.db",
                     annot = annFUN.org,
                     ID = GO.ID,
                     mapping = "org.Hs.eg",
                     geneSelectionFun = function(x) {
                       x == numb.col[moduleName]}
                     )
    # save(GOdata.bp, file = "BP_GO.RData")
    print(GOdata.bp)
    go.enrich(GOdata.bp, moduleName, "BP")

    # GOdata.mp <- new("topGOdata",
    #                  ontology = "MP",
    #                  description = "Molecular process of the genes in the study.",
    #                  allGenes = genes,
    #                  # annot = annFUN.gene2GO, ## the new annotation function
    #                  # affyLib = "org.Hs.eg.db",
    #                  annot = annFUN.org,
    #                  ID = GO.ID,
    #                  mapping = "org.Hs.eg",
    #                  geneSelectionFun = function(x){x == numb.col[moduleName]})
    # save(GOdata.mp, file = "MP_GO.RData")
    # tryCatch({
    #   geneSelectionFun(GOdata.mp) <- function(x){x == numb.col[moduleName]}
    #   print(GOdata.mp)
    #   go.enrich(GOdata.mp, moduleName, "MP")},
    #   error = function(x) {
    #   message("\nUnable to calculate MP.")
    #   message(x)
    # })
    GOdata.cc <- new("topGOdata",
                     ontology = "CC",
                     description = "Cellular component of the genes in study",
                     allGenes = genes,
                     # annot = annFUN.gene2GO, ## the new annotation function
                     # affyLib = "org.Hs.eg.db",
                     annot = annFUN.org,
                     ID = GO.ID,
                     mapping = "org.Hs.eg",
                     geneSelectionFun = function(x) {
                       x == numb.col[moduleName]}
                     )
    # save(GOdata.cc, file = "CC_GO.RData")
    print(GOdata.cc)
    go.enrich(GOdata.cc, moduleName, "CC")
  }
  # Reactome ####
  if (Reactome) {
    message("\tDoing the reactome analysis")
    reactome_enrich <- enrichPathway(gene = moduleGenes,
                                     universe = universeGenesEntrez,
                                     pvalueCutoff = 0.05, readable = TRUE,
                                     minGSSize = 2, maxGSSize = 2000)
    if (is.null(reactome_enrich)) {
      message("Module ", moduleName, " is not enriched in a Reactome pathway.")
    } else if (nrow(reactome_enrich) >= 1) {
      write.csv(as.data.frame(reactome_enrich),
                file = paste0("reactome_", moduleName, ".csv"),
                row.names = FALSE, na = "")
      pathway.enrich(reactome_enrich, moduleName)
    } else {
      message("\tWithout significant pathways")
    }
    message("\tFinishing the reactome analysis")
  }

  #  Kegg ####
  if (Kegg) {
    message("\tDoing the Kegg analysis")
    kegg_enrich <- enrichKEGG(moduleGenes,
                              universe = universeGenesEntrez,
                              minGSSize = 2, maxGSSize = 2000)
    if (is.null(kegg_enrich)) {
      message("Module ", moduleName, " is not enriched in a KEGG pathway.")
    } else if (nrow(kegg_enrich) >= 1) {
      write.csv(as.data.frame(kegg_enrich),
                file = paste0("kegg_", moduleName, ".csv"),
                row.names = FALSE, na = "")
      pathway.enrich(kegg_enrich, moduleName)
    } else {
      message("\tWithout significant pathways")
    }
    message("\tFinishing the kegg analysis")
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
      message("\tDoing the STRING analysis")
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
      warning(paste("Module", moduleName, "is bigger!"))
    }
  }
})

