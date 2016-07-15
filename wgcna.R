# Analyse data with WGCNA from doi:10.1136/gutjnl-2011-301146.

library("GEOquery", quietly = T)

# Prepare some variables
gse_number <- "GSE28619"

# Download files
files.info <- getGEO(gse_number, destdir = "~/Documents")
print(files.info)
info <- files.info$GSE28619_series_matrix.txt.gz
disease <- as.numeric(pData(info)$description) - 1
dim(phenoData(info))

data.wgcna <- t(exprs(info))

save(data.wgcna, disease, file = "Input.RData")
