# Analyse data with WGCNA from doi:10.1136/gutjnl-2011-301146.

library("GEOquery", quietly = T)
library("WGCNA")

# Prepare some variables
gse_number <- "GSE28619"
#
# # Download files
# files.info <- getGEO(gse_number, destdir = "~/Documents")
# print(files.info)
# info <- files.info$GSE28619_series_matrix.txt.gz
#
#
# disease <- as.numeric(pData(info)$description) - 1
# dim(phenoData(info))
load("Input.RData")

samples <- c("N57", "C9", "C7", "C12", "N54", "N17", "N53", "CA64", "CA45",
             "CA62", "CA51", "CA79", "CA92", "CA47", "CA20", "CA22", "CA82",
             "CA39", "CA80", "CA28", "CA77", "CA5")
names(samples) <- rownames(data.wgcna)

disease <- read.csv("clean_variables.csv")
ids <- disease$id
disease <- disease[, -c(1, 5, 6, 9, 42)]
disease.r <- apply(disease, 2, as.numeric)
nam <- c("status_90", "gender", "infection_hospitalization", "aki", "hvpg_corte20",
  "hvpg_corte20", "lille_corte", "ch")
for (n in nam) {
  disease.r[,n] <- as.factor(disease[,n])
}
disease <- disease.r
disease <- disease[ids %in% samples, ]
# data.wgcna <- t(exprs(info))

pdf("dendro_traits.pdf")
# Re-cluster samples
sampleTree2 = hclust(dist(data.wgcna[samples %in% ids, ]), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(disease, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(disease),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
pdf("dendro.pdf")

 # Plot the sample dendrogram
plotClusterTreeSamples(sampleTree2, 
                    main = "Sample dendrogram and trait heatmap")
dev,off()
save(data.wgcna, disease, samples, ids, file = "Input.RData")


