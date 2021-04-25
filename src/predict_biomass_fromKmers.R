# predict_biomass_fromKmers.R
suppressMessages(library(data.table))
suppressMessages(library(Biostrings))
args = commandArgs(trailingOnly=TRUE)

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

kmer_dat <- readRDS(paste0(script.dir,"/../dat/seq/16S_db/kmer5_refseq.RDS"))
kmer_dat$kmer_freq_norm <- kmer_dat$kmer_freq / rowSums(kmer_dat$kmer_freq) * 100

n <- nrow(kmer_dat$kmer_freq_norm)
a <- readDNAStringSet(args[1])
a_kmer <- colSums(oligonucleotideFrequency(a, kmer_dat$k))

res <- unlist(lapply(1:n, FUN = function(x) cor(kmer_dat$kmer_freq_norm[x,],
                                                (a_kmer/sum(a_kmer))*100)))

tophit <- order(-res)[1]

bm_pred <- kmer_dat$phenotype[tophit]

if(bm_pred == "archaeal")
  bm_pred <- "Archaea"
if(bm_pred == "Gram+")
  bm_pred <- "Gram_pos"
if(bm_pred == "Gram-")
  bm_pred <- "Gram_neg"



cat(bm_pred)
