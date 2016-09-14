

library("plyr")
library("org.Hs.eg.db")
read.genes <- function(file, trait, th = 0){
  dat <- read.csv(file)
  c.dat <- colnames(dat)
  ncol <- grep(trait, c.dat)
  if (length(ncol) != 2) {
    stop("Too much traits")
  }
  g.names <- dat[, "SYMBOL"][dat[, ncol[2]] <= 0.05 & abs(dat[, ncol[2]]) >= th]
  g.names <- as.character(g.names)
  g.out <- AnnotationDbi::select(org.Hs.eg.db, keys = g.names,
                                 columns = "GENENAME", keytype = "SYMBOL")
  g.out <- unique(g.out[, "GENENAME"])
  g.out[!is.na(g.out)]
}

info <- list("meld" = c("skyblue", "midnightblue", "darkolivegreen"),
             "maddrey" = c("darkolivegreen", "skyblue", "steelblue"),
             "status_90" = c("lightcyan", "salmon"))

skyblue_meld <- read.genes("skyblue_trait.csv", "meld")
midnightblue_meld <- read.genes("midnightblue_trait.csv", "meld")
darkolivegreen_meld <- read.genes("darkolivegreen_trait.csv", "meld")
meld <- rbind.fill(as.data.frame(t(skyblue_meld)),
                   as.data.frame(t(midnightblue_meld)),
                   as.data.frame(t(darkolivegreen_meld)))
meld <- t(meld)
colnames(meld) <- c("skyblue_meld", "midnightblue_meld", "darkolivegreen_meld")

skyblue_maddrey <- read.genes("skyblue_trait.csv", "maddrey")
steelblue_maddrey <- read.genes("steelblue_trait.csv", "maddrey")
darkolivegreen_maddrey <- read.genes("darkolivegreen_trait.csv", "maddrey")
maddrey <- rbind.fill(as.data.frame(t(skyblue_maddrey)),
                      as.data.frame(t(steelblue_maddrey)),
                      as.data.frame(t(darkolivegreen_maddrey)))
maddrey <- t(maddrey)
colnames(maddrey) <- c("skyblue_maddrey", "steelblue_maddrey",
                       "darkolivegreen_maddrey")

lightcyan_status_90 <- read.genes("lightcyan_trait.csv", "status_90")
salmon_status_90 <- read.genes("salmon_trait.csv", "status_90")
status_90 <- rbind.fill(as.data.frame(t(lightcyan_status_90)),
                        as.data.frame(t(salmon_status_90)))
status_90 <- t(status_90)
colnames(status_90) <- c("lightcyan_status_90", "salmon_status_90")
dat.complete <- rbind.fill(as.data.frame(t(status_90)),
                           as.data.frame(t(maddrey)),
                           as.data.frame(t(meld)))
dat.all <- t(dat.complete)
colnames(dat.all) <- c(colnames(status_90),
                       colnames(maddrey),
                       colnames(meld))
write.csv(dat.all, "genes_modules_trait.csv", row.names = FALSE, na = "")
