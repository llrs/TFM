#  Analyse the modules with globaltest o GOstats, SPIA, hopach
# First steps similar to the DE analysis

lnames <- load(file = "TNF_AH-network-auto.RData") 
lnames
lnames <- load(file = "Input.RData") 
lnames
lnames <- load(file = "ME.RData")
lnames
lnames <- load(file = "Module_info.RData")
lnames

library("globaltest")
library("limma")
# library("SPIA")
# library("GOstats")
# library("topGO")
# library("ggbio")
# library("affy")
# library("RColorBrewer")
# library("gcrma")
library("sva")
library("svd")
library("hgu133plus2.db")

moduleName <- "lightcoral"
exprs <- t(data.wgcna)
colnames(exprs) <- samples

# annots <- select(hgu133plus2.db, keys = rownames(exprs),
#                  columns = c("GO", "SYMBOL", "GENENAME"), keytype = "PROBEID")
# save(annots, file = "annots_study.RData")
load(filee = "annots_study.RData")
myGenes <- rownames(exprs)[net$colors  ==  moduleName]
myGenes2GO <- sapply(myGenes, function(x){
  annots[annots$PROBEID == x & annots$ONTOLOGY == "MF", "GO"]
})

genes <- ifelse(rownames(exprs) %in% myGenes, 1, 0)
names(genes) <- rownames(exprs)

selFun <- function(x){
  if (x == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

GOdata <- new("topGOdata",
              ontology = "MF",
              description = paste("Molecular function of the", 
                                  moduleName, "module."),
              allGenes = genes,
              annot = annFUN.gene2GO, ## the new annotation function 
              gene2GO = myGenes2GO,
              geneSelectionFun = selFun) 

resultFisher <- runTest(GOdata, 
                        algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

showSigOfNodes(GOdata, 
               score(resultKS.elim), firstSigNodes = 2, useInfo = 'all')
