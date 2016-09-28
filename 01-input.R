# Read the files

source("/home/lrevilla/Documents/TFM/00-general.R", echo = TRUE)

pheno <- read.csv("miRNA_phenotype.csv")
exprs <- read.csv("miRNA_normQ.csv", row.names = 1)
exprs2 <- read.table("normalised_ASH_mirna.csv", sep = "\t")
vclin <- read.csv(file.path("..", "DB_MIR-Alcoholic_Hepatitis_Sept2015.csv"))
vclin2 <- read.csv(file.path("..", "DB_integrated.CCL20-Alcoholic_Hepatitis_19_06_2013_(2)_1_.csv"))
setwd(data.files.out)
# ==============================================================================
#
#  Code chunk 1: Preparing the clinical data
#
# ==============================================================================


AH <- paste0("HA-", c(23, 24, 30, 36, 39, 51, 54, 57, 60, 77, 91, 84, 10, 90,
                      18, 21, 25, 34, 40, 50, 2, 15, 26, 28, 33, 37, 52, 86,
                      87))
normal <- paste0("NQ-", c(2:5))
ids <- c(AH, normal)
samplename <- substring(pheno$sampleid, 8)

samplename <- convert(samplename)
pheno2 <- cbind(pheno, "samplename" = samplename,
               "group" = as.factor(substring(pheno$GROUP, 1, 2)),
               "patientid" = ids)

vclin <- merge(pheno2, vclin, by.x = "patientid", by.y = "id", all.x = TRUE)
inter_v <- c("patientid", "sampleid", "samplename", "group", "age",
             "gender", "meld", "lille", "glucose", "status_180",
             "ast", "alt")
am <- c("status_90", "infection_hospitalization", "aki", "hvpg_corte20",
        "hvpg", "maddrey", "status_1year", "creatinine", "ap", "ggt",
        "hb_g.dl", "tp")
vars <- c(am, inter_v)
pheno <- vclin[, colnames(vclin) %in% vars]

pheno$status_90 <- ifelse(pheno$status_90 == "alive", 1, 0)
pheno$status_180 <- ifelse(pheno$status_180 == "alive", 1, 0)
pheno$status_1year <- ifelse(pheno$status_1year == "alive", 1, 0)
pheno$gender <- fact2num(pheno$gender, "Male", 1)
pheno$gender <- fact2num(pheno$gender, "Female", 0)

ids_2 <- sub("_\\..*", "", colnames(exprs2))
ids_2 <- sub("\\.", "", ids_2)
vclin2 <- vclin2[vclin2$id %in% ids_2, colnames(vclin2) %in% vars]

v <- rbind.fill(pheno, vclin2)
# ==============================================================================
#
#  Code chunk 2: Preparing the expression data
#
# ==============================================================================

data.wgcna <- t(exprs[, pheno$group  == "AH"])
gsg <- goodSamplesGenes(data.wgcna, verbose = 3)
if (!gsg$allOK)
{ # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:",
                     paste(names(data.wgcna)[!gsg$goodGenes],
                           collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:",
                     paste(rownames(data.wgcna)[!gsg$goodSamples],
                           collapse = ", ")))
  # Remove the offending genes and samples from the data:
  data.wgcna <- data.wgcna[gsg$goodSamples, gsg$goodGenes]
}


pdf("samples_dendrogram.pdf")
sampleTree <- hclust(dist(data.wgcna), method = "average")
pars <- par(mar = c(0, 4, 2, 0), cex = 0.6)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

save(pheno, exprs, file = "Input.RData")
