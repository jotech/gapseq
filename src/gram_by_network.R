#
# Use similarity to metabolic networks of reference strains to infer gram staining
#
library(data.table)

find_gram_by_network <- function(file, script.dir){
  dcast.gapseq <- readRDS(paste0(script.dir, "/../dat/", "dist_reference-networks.RDS"))
  ref.info <- fread(paste0(script.dir, "/../dat/", "reference_genomes_edited.tbl"))
  
  dt.new <- fread(file)
  org    <- gsub("-all-Pathways.tbl","",basename(file))
  dt.new$org <- org
  
  dcast.new<- dcast.data.table(dt.new, formula=org~ID, value.var = "Prediction")
  idx.new <- which(colnames(dcast.new) %in% colnames(dcast.gapseq))
  idx.old <- which(colnames(dcast.gapseq) %in% colnames(dcast.new))
  dcast.all <- rbind(dcast.gapseq[,..idx.old],dcast.new[,..idx.new])
  mat.gapseq   <- as.matrix(dcast.all[,-1])
  rownames(mat.gapseq) <- dcast.all$org
  dist.gapseq <- dist(mat.gapseq)
  
  dist.mat <- as.matrix(dist.gapseq)
  idx <- match(org, rownames(mat.gapseq))
  dist.mat[idx,idx] <- Inf # tmp hack
  closest <- which.min(dist.mat[,idx])
  
  hit <- ref.info[make.names(name)==rownames(dist.mat)[closest],.(name,gram)]
  
  cat("Closest metabolic network:", hit$name, "\tgram:", hit$gram, "\n")  
  return(hit$gram)
}
