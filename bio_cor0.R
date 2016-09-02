library("biomaRt")
library("hgu133plus2.db")
library("testthat")
library("GOstats")
library("graphite")
library("WGCNA")
library("KEGGgraph")
library("KEGG.db")
library("RBGL")

".combinadic" <- function(n, r, i) {
  
  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic
  n0 <- length(n)
  if (i < 1 | i > choose(n0,r)) stop("'i' must be 0 < i <= n!/(n-r)!")
  largestV <- function(n, r, i) {
    #v <- n-1
    v <- n                                  # Adjusted for one-based indexing
    #while(choose(v,r) > i) v <- v-1
    while (choose(v,r) >= i) v <- v - 1        # Adjusted for one-based indexing
    return(v)
  }
  
  res <- rep(NA,r)
  for (j in 1:r) {
    res[j] <- largestV(n0,r,i)
    i <- i - choose(res[j],r)
    n0 <- res[j]
    r <- r - 1
  }
  res <- res + 1
  res <- n[res]
  return(res)
}

compare_graphs <- function(g1, g2){
  # Function to estimate how much two graphs overlap by looking if the nodes
  # are the same
  prot1 <- nodes(g1)
  prot2 <- nodes(g2)
  if (length(prot1) == 0 | length(prot2) == 0) {
    return(NA)
  }
  score <- (length(intersect(prot1, prot2)))*2/(
    length(prot2) + length(prot1))
}

go_cor <- function(e_a, e_b, chip = "hgu133plus2.db", ...){
  # Calculates the degree of overlap of the GO BP ontologies of entrez ids.
  # https://support.bioconductor.org/p/85702/#85732
  LP <- simLL(e_a, e_b, "BP", measure = "LP", chip = chip, ...)
  UI <- simLL(e_a, e_b, "BP", measure = "UI", chip = chip, ...)
  
  if (is.na(LP) | is.na(UI)) {
    return(NA)
  }
  s.path <- function(ig){
    # The longest of the shortest path of a graph
    lfi <- leaves(ig, "in")
    degs <- degree(ig)
    root <- names(degs$outDegree)[degs$outDegree == 0]
    paths <- sp.between(ig, lfi, root)
    plens <- subListExtract(paths, "length", simplify = TRUE)
    max(plens)
  }
  # Calculates the score taking into account the size and the middnest path
  # Taking advantage of the fact that in GO there is a root and leaves
  (UI$sim/LP$sim)*max(s.path(LP$g1), s.path(LP$g2))
}

kegg_cor <- function(react_a, react_b){
  # Function that correlates based on reactome ids
  # Basically calculates how many nodes do overlap between pathways
  
  if (is.na(react_a) | is.na(react_b)) {
    return(NA)
  } else if (react_a == react_b) {
    return(1)
  }
  
  # Retrieve the results from internet https site
  tmp <- tempfile("hsa")
  retrieveKGML(paste0("path:", react_a), "hsa", tmp, method = "wget", 
               quiet = TRUE)
  g1 <- parseKGML2Graph(tmp, expandGenes = TRUE)
  
  tmp2 <- tempfile("hsa")
  retrieveKGML(paste0("path:", react_b), "hsa", tmp2, method = "wget", 
               quiet = TRUE)
  g2 <- parseKGML2Graph(tmp2, expandGenes = TRUE)
  
  score <- compare_graphs(g1, g2)
  score
}

# Not useful because it doesn't have any real causation 02/08/2016
dist_cor <- function(a, b, info){
  use_info <- c("chromosome_name", "strand", "start_position", "end_position",
                "gene_biotype")
  info_a <- unique(info[info$affy_hg_u133_plus_2 == a, use_info])
  info_b <- unique(info[info$affy_hg_u133_plus_2 == b, use_info])
  if (nrow(info_a) != 1 | nrow(info_b) != 1) {
    return(NA)
  }
  # Using the position and the type of genes output a cor
  if (info_a["chromosome_name"] != info_b["chromosome_name"]) {
    score <- 0
  }
  score <- 0.5 # If in the same chromosome at least 0.5 distance correlation
  
  # Maybe the default for the same chromosome could be substituted by the
  # CM distance
  # 5*10^5 is the region of upstream/downstream where regulation usually occurs
  if (info_b["strand"] == info_a["strand"]) {
    start_dist <- info_a["start_position"] - info_b["start_position"]
    end_dist <- info_a["end_position"] - info_b["end_position"]
  } else {
    start_dist <- info_a["start_position"] - info_a["end_position"]
    end_dist <-  info_b["start_position"] - info_b["end_position"]
  }
  
  if (abs(start_dist) < 5*10 ^ 5) {
    score <- score + 0.25
  } else if (abs(end_dist) < 5*10 ^ 5) {
    score <- score + 0.25
  }
  
  # Using the cathegory of biotypes to score them
  biotypes <- c("miRNA", "snRNA", "snoRNA", "scaRNA", "lincRNA", "lncRNA")
  if (((info_a["gene_biotype"] %in% biotypes)  &
       (info_b["gene_biotype"] == "protein_coding")) |
      ((info_b["gene_biotype"] %in% biotypes) |
       (info_a["gene_biotype"] == "protein_coding"))) {
    score <- score + 0.25
  }
  if ( (info_a["gene_biotype"] == "processed_pseudogene") |
       (info_b["gene_biotype"] == "processed_pseudogene")) {
    score <- score - 0.25
  }
  return(score)
}

react_cor <- function(react_a, react_b, hR){
  # Function that correlates based on reactome ids
  # Basically calculates how many nodes do overlap between pathways
  
  if (is.na(react_a) | is.na(react_b)) {
    return(NA)
  } else if (react_a == react_b) {
    return(1)
  }
  
  ids <- sapply(hR, function(x){x@id})
  react_name <- names(hR[ids %in% react_a])
  react_name2 <- names(hR[ids %in% react_b])
  
  if (length(react_name) != 1 | length(react_name2) != 1) {
    return(NA)
  }
  
  # Obtain each pathway
  g1 <- hR[[react_name]]
  g2 <- hR[[react_name2]]
  
  score <- compare_graphs(g1, g2)
  score
}

comb2mat <- function(input, func, ...){
  # Funcion to perform efficiently the conversion from combinations
  #  to symmetric matrix
  # Perform all the combinations of 2 from the input
  cobs <- list()
  for (i in 1:length(input)) {
    cobs[[i]] <- .combinadic(input, 2, i)
  }
  # cobs <- combn(input, 2)
  
  N <- sapply(cobs, function(x, ...){func(x[1], x[2], ...)})
  # Function that performs the calculus
  # N <- seq_len(ncol(combs))
  out <- matrix(ncol = length(input), nrow = length(input))
  out[lower.tri(out)] <- N
  out <- t(out)
  out[lower.tri(out)] <- N
  out <- t(out)
  diag(out) <- 1
  rownames(out) <- colnames(out) <- input
  return(out)
}

seq2mat <- function(x, dat) {
  out <- matrix(ncol = length(x), nrow = length(x))
  out[lower.tri(out)] <- unlist(dat)
  out <- t(out)
  out[lower.tri(out)] <- unlist(dat)
  out <- t(out)
  diag(out) <- 1
  rownames(out) <- colnames(out) <- x
  return(out)
}

# Extract all the id of reactome for each combination and compare them all
comb_biopath <- function(comb, info, by, biopath){
  # react_path <- apply(comb, 2, function(y){
  a <- unique(info[info[by] == comb[1], biopath])
  a <- a[a != ""]
  
  b <- unique(info[info[by] == comb[2], biopath])
  b <- b[b != ""]
  if (all(sapply(a, is.na)) | all(sapply(b, is.na))) {
    return(NA)
  }
  expand.grid(a, b)
  # })
}

# Check if something is a matrix (internal use only)
check_na <- function(x){
  if (class(x) != "matrix") {
    if (length(x) == 0) {
      return(TRUE)
    } else if (length(x) == 1) {
      if (is.na(x)) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

bio.cor <- function(x, ... ){
  # Using data correlates biologically two genes or probes
  # From the graphite package
  names_probes <- x
  humanReactome <- pathways("hsapiens", "reactome")
  # humanKegg <- pathways("hsapiens", "kegg")
  
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  attri <- c("affy_hg_u133_plus_2","entrezgene",
             "gene_biotype", "start_position", "end_position", "chromosome_name",
             "strand", "reactome",
             "hgnc_symbol")
  info <- getBM(attributes = attri,
                filters = "affy_hg_u133_plus_2", values = names_probes,
                mart = mart)
  
  affy.id <- select(hgu133plus2.db, keys = names_probes,
                    columns = c("PROBEID", "ENTREZID", "SYMBOL", "PATH"))
  
  # go_mat <- comb2mat(names_probes, go_cor)
  # dist_mat <- comb2mat(x, dist_cor, info) # Not useful
  
  # TODO: extract a lab notebook
  # comb <- combn(names_probes, 2)
  react.bio <- rep(NA, choose(length(names_probes), 2))
  kegg.bio <- react.bio
  for (i in 1:choose(length(names_probes), 2)) {
    comb <- .combinadic(names_probes, 2, i)
    react_path <- comb_biopath(comb, info, "affy_hg_u133_plus_2","reactome")
    # Calculate the max of the correlation of both ids
    if (check_na(react_path)){
      a <- NA
    } else {
      a <- apply(react_path, 1, function(y){
        react_cor(y[1], y[2], hR = humanReactome)
        }
      )
    }
    if (length(a) != 0) {
      react.bio[i] <- max(a, na.rm = TRUE)
    } else {
      react.bio[i] <- NA
    }
    # # })
    # print(head(comb))
    # print(head(affy.id))
    kegg_path <- comb_biopath(comb, affy.id, "PROBEID","PATH")
    # Calculate the max of the correlation of both ids
    # N <- lapply(kegg_path, function(x){
    #   if (check_na(x)) {
    #     return(NA)
    #   }
    if (check_na(kegg_path)) {
      a <- NA
    } else {
      # print(head(kegg_path))
      a <- apply(kegg_path, 1, function(y) {
        result <- tryCatch({
          kegg_cor(y[1], y[2])
          }, warning = function(w) {
            message(paste("Warning downloading the data for ", 
                          y[1], "or", y[2]))
            message(w)
          }, error = function(e) {
            message(paste("Couldn't download the data for ", y[1], "or", y[2]))
            # Choose a return value in case of error
            return(NA)
          })
        }
      )
    }
    
    if (length(a) != 0) {
      kegg.bio[i] <- max(a, na.rm = TRUE)
    } else {
      kegg.bio[i] <- NA
    }
    # }})
    
  }
  react_mat <- seq2mat(names_probes, react.bio)
  kegg_mat <- seq2mat(names_probes, kegg.bio)
  cor_mat <- list(reactome = react_mat, kegg = kegg_mat, go = go_mat)
  
}

weighted <- function(x, w){
  if (length(x) != length(w)) {
    stop(paste("Weights and data don't match the length.\n",
               length(x), length(w)))
  }
  sum(x*w, na.rm = TRUE)
}

# Function that used the previously calculated biological correlation to
# calculate the total correlation
cor.all <- function(datExpr, bio_mat, ...){
  cor_mat <- cor(datExpr, use = "p")
  cors <- c(cor_mat, bio_mat)
  apply(simplify2array(cors), c(1,2), weighted, w = c(0.5, 0.18, 0.10, 0.22))
}