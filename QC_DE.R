# Analyse data from article

# Load required libraries
library("GEOquery")
library("affy")
library("simpleaffy")
library("edgeR")
library("limma")
library("gcrma")


# Prepare some variables
gse_number <- "GSE28619"
path_files <- "../Documents/GSE28619_RAW/"
test_file <- "GSM709348.CEL.gz"
path_file <- paste0(path_files, test_file)
setwd("~/Documents")
raw.tar <- paste0(gse_number, "_RAW.tar")
path_raw <- paste(gse_number, raw.tar, sep="/")
#By exploring the data on GEO I found that the platform is :
# GPL570 	[HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
# So we need to load the bioconductor library with the information about it
library("Affyhgu133Plus2Expr")

# Tutorial on https://www.biostars.org/p/53870/
# Based on this tutorial:
# http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
# Download set
getGEOSuppFiles(gse_number)

#Unpack the CEL files
untar(path_raw, exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
cels  <-  list.files("data/", pattern = "CEL")
setwd("data/")

# Phenodata contains the information of the experiment
celfiles <- read.affy("phenodata.txt")

#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# From https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE28619
gset <- getGEO("GSE28619", GSEMatrix =TRUE)
if (length(gset) > 1) tdx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
expres <- exprs(gset)
