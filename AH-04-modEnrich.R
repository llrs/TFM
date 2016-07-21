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
library("topGO")
library("Rgraphviz")
# library("ggbio")
# library("affy")
# library("RColorBrewer")
# library("gcrma")
library("sva")
library("svd")
library("hgu133plus2.db")
library("WGCNA")
enableWGCNAThreads(6)

moduleName <- "lightcoral"
exprs <- t(data.wgcna)
colnames(exprs) <- samples

# # Not needed since we can use it directly on the topGOdata creation
# annots <- select(hgu133plus2.db, keys = rownames(exprs),
#                  columns = c("GO", "SYMBOL", "GENENAME"), keytype = "PROBEID")
# save(annots, file = "annots_study.RData")
# load(file = "annots_study.RData")
# myGenes <- rownames(exprs)[moduleColors  ==  moduleName]
# genesGO <- function(probeid, annots, ontology = c("MF", "BP", "CC")){
#   myGenes2GO <- sapply(probeid, function(x){
#     annots[annots$PROBEID == x & annots$ONTOLOGY == "MF", "GO"]
#   })
#   return(myGenes2GO)
# }

# myGenes2GO <- genesGO(myGenes, annots, "BP")
# allGenes2GO <- genesGO(rownames(exprs), annots, "BP")
# save(allGenes2GO, file = "genesGO.RData")
# load(file = "genesGO.RData")

# I assume I keep the same order of genes
genes <- as.numeric(as.factor(moduleColors))
names(genes) <- rownames(exprs)
n <- length(unique(moduleColors))

moduleSel <- function(modul, n){
  selFun <- function(x){
    number <- grep(paste0("^", modul, "$"),
                   labels2colors(0:n))
    x == number
  } 
}
selFun <- moduleSel(moduleName, n)

# GOdata <- new("topGOdata",
#               ontology = "BP",
#               description = paste("Molecular function of the", 
#                                   moduleName, "module."),
#               allGenes = genes,
#               annot = annFUN.db , ## the new annotation function 
#               affyLib = "hgu133plus2.db",
#               geneSelectionFun = selFun)
# 
# save(GOdata, file = "array_BP.RData")
load(file = "array_BP.RData")

resultFisher <- runTest(GOdata, 
                        algorithm = "classic", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata, classic = resultFisher, 
                  orderBy = "classic", ranksOf = "classic")
write.csv(allRes, file = paste0("table_", moduleName, ".csv"), 
          row.names = FALSE)

pdf(paste0("BP_GO_fisher_", moduleName, ".pdf"))
showSigOfNodes(GOdata, 
               score(resultFisher), firstSigNodes = 2, useInfo = 'all')
dev.off()


moduleName <- "plum"
geneSelectionFun(GOdata) <- moduleSel(moduleName, n)
resultFisher <- runTest(GOdata, 
                        algorithm = "classic", statistic = "fisher")

pdf(paste0("BP_GO_fisher_", moduleName, ".pdf"))
showSigOfNodes(GOdata, 
               score(resultFisher), firstSigNodes = 2, useInfo = 'all')
dev.off()