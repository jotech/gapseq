#!/usr/bin/Rscript

# first argument: genome file
# second argument: annotation file

args = commandArgs(trailingOnly=TRUE)
if ( length(args)!=2 & length(args)!=3 ) {
  stop("First argument should be a genome file and the second one an gapseq reaction file.")
}

suppressMessages(library(data.table))
suppressMessages(library(Biostrings))
suppressMessages(library(IRanges))

coverage <- function(genome.filestr, anno.filestr, anno.filestr2=NA){
  #genome <- readDNAStringSet("~/uni/gapseq/toy/ecoli.fna.gz")
  #anno   <- fread("~/uni/gapseq/toy/ecoli-all-Reactions.tbl", fill = T, skip = "rxn	")
  genome   <- readDNAStringSet(genome.filestr)
  anno.raw <- fread(anno.filestr, fill = T, skip = "rxn	")
  anno     <- unique(anno.raw[,.(sstart, send, stitle)])
  if( !is.na(anno.filestr2) ){
    anno.raw2<- fread(anno.filestr2, fill = T, skip = "id	")  
    anno2    <- unique(anno.raw2[,.(sstart, send, stitle)])
    anno     <- rbind(anno, anno2)
  }
  anno <- anno[!is.na(sstart) & !is.na(send)]
  anno[,start:=ifelse(sstart<send, sstart, send)]
  anno[,end:=ifelse(send>sstart, send, sstart)]
  
  # reduce overlapping intervals 
  anno.red <- data.table()
  for(s in unique(anno$stitle)){
    anno[stitle==s]
    tmp <- reduce(IRanges(start=anno[stitle==s, start], end=anno[stitle==s, end]))
    anno.red <- rbind(anno.red, data.table(stitle=s, start=tmp@start, end=tmp@start+tmp@width-1))
  }

  genome.len <- sum(sapply(genome, length))
  
  genome.hit <- lapply(genome, function(x){vector("numeric", length=length(x))})
  for(i in 1:nrow(anno.red)){
    feat.start <- anno.red[i,start]
    feat.end   <- anno.red[i,end]
    contig     <- anno.red[i,stitle]
    contig.idx <- match(contig, names(genome.hit))
    
    if( !is.na(feat.start) & !is.na(feat.end) ){
      genome.hit[[contig.idx]][feat.start:feat.end] <- 1  
    }
  }
  genome.cov <- sum(sapply(genome.hit, sum)) / genome.len
  cat("Annotation/genome coverage:", genome.cov, "\n")
  
  return(invisible(genome.cov))
}
#coverage("~/uni/gapseq/toy/ecoli.fna.gz", "~/uni/gapseq/toy/ecoli-all-Reactions.tbl", "~/uni/gapseq/toy/ecoli-Transporter.tbl") # 0.342
#coverage("~/uni/gapseq/toy/myb71.fna.gz", "~/uni/gapseq/toy/myb71-all-Reactions.tbl", "~/uni/gapseq/toy/myb71-Transporter.tbl") # 0.2643

if( length(args) == 2){
  coverage(args[1], args[2])  
}else{
  coverage(args[1], args[2], args[3])
}
