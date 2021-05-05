#
# Use similarity to metabolic networks of reference strains to infer gram staining
#
library(data.table)

find_gram_by_network <- function(file, script.dir){
  dcast.gapseq <- readRDS(paste0(script.dir, "/../dat/", "dist_reference-networks.RDS"))
  ref.info <- fread(paste0(script.dir, "/../dat/", "reference_genomes_edited.tbl"))
  
  dt.new <- fread(file)
  dt.new[,Prediction:=Prediction*1]
  org    <- gsub("-all-Pathways.tbl","",basename(file))
  
  dcast.new <- matrix(nrow=1, ncol=nrow(dt.new), data=dt.new$Prediction); colnames(dcast.new) <- dt.new$ID
  #dcast.new<- dcast.data.table(unique(dt.new), formula=~ID, value.var = "Prediction")
  comm.pwy <- intersect(colnames(dcast.new), colnames(dcast.gapseq))
  idx.old <- match(comm.pwy, colnames(dcast.gapseq))
  idx.new <- match(comm.pwy, colnames(dcast.new))
  mat.gapseq <- as.matrix(rbind(dcast.gapseq[,c(idx.old)],dcast.new[,c(idx.new)]))
  rownames(mat.gapseq)[nrow(mat.gapseq)] <- org
  
  dist.gapseq <- dist(mat.gapseq)
  
  dist.mat <- as.matrix(dist.gapseq)
  idx <- match(org, rownames(mat.gapseq))
  dist.mat[idx,idx] <- Inf # tmp hack
  closest <- which.min(dist.mat[,idx])
  
  hit <- ref.info[name==rownames(dist.mat)[closest],.(name,gram)]
  
  cat("Closest metabolic network:", hit$name, "\tgram:", hit$gram, "\n")  
  return(hit$gram)
}
