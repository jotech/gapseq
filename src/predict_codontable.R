suppressMessages(library(Biostrings))
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
genome <- args[1]
# genome <- "GCF_000350285.1_OR1_genomic.fna" # 25
# genome <- "HRGMv2_2645.fna" # 4
# genome <- "HRGMv2_2676.fna" # 11
# genome <- "HRGMv2_0240.fna" # 11

codes <- c("4","11","25")

prefix <- basename(sub("\\.fna$","",genome))
fna <- readDNAStringSet(genome)
names(fna) <- sub(" .*$","",names(fna))

codestat <- lapply(codes, function(code) {
  gff <- paste0(prefix,"_",code,".gff")
  # file.copy(gff,paste0("~/tmp/",prefix,"_",code,".gff"))
  faa <- readAAStringSet(paste0(prefix,"_",code,".faa"))
  concatAA <- paste(unname(as.character(faa)), collapse = "")
  concatAA <- gsub("\\*","",concatAA)
  ctsW <- str_count(concatAA,"W")
  ctsG <- str_count(concatAA,"G")
  df <- read.csv(gff, sep = "\t", header = FALSE, comment.char = "#")
  df <- df[!grepl("^#",df$V1),]
  contigs <- unique(names(fna))
  cov <- lapply(contigs, function(contig) {
    fnawidth <- width(fna[contig])
    idx <- which(df$V1 == contig)
    if(length(idx) > 0) {
      tmpmat <- as.matrix(df[idx,c("V4","V5")], drop = FALSE)
      ranges <- IRanges(start = apply(tmpmat,1,min),
                        end = apply(tmpmat,1,max))
      ranges <- reduce(ranges)
      codinglength <- sum(width(ranges))
      return(c(codinglength,fnawidth))
    } else {
      return(c(0, fnawidth))
    }
  })
  cov <- do.call("rbind",cov)
  cov <- sum(cov[,1])/sum(cov[,2])
  return(c(Coverage = round(cov*100),
           percW = round(ctsW/str_length(concatAA), digits = 5),
           percG = round(ctsG/str_length(concatAA), digits = 5)))
})
names(codestat) <- codes

usecode <- "11"

if(codestat[["4"]]["Coverage"] - codestat[["11"]]["Coverage"] > 5) {
  wperc <- c("4"  = unname(codestat[["4"]]["percW"]),
             "25" = unname(codestat[["25"]]["percW"]))
  usecode <- names(wperc)[which.min(abs(wperc-0.01))]
}

writeLines(usecode, con = paste0(prefix,"_code"))
