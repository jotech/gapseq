suppressMessages(library(data.table))

#-------------------------------------------------------------------------------
# (0) Read pre-alignment data and alignment stats
#-------------------------------------------------------------------------------
load("~/tmp/prealignment_data.RData")

alicols <- c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","sseq")
alignments <- fread("~/tmp/alignments.tsv")
setnames(alignments, alicols[1:ncol(alignments)])
alignments[, file := sub("\\|.+$","",qseqid)]
alignments[, qseqid := sub("^.+\\.fasta\\|","",qseqid)]

#-------------------------------------------------------------------------------
# (1) Parse arguments and gapseq/script path
#-------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

# command line arguments
bitcutoff <- as.numeric(args[1])
identcutoff <- as.numeric(args[2])
strictCandidates <- args[3] == "true"
output_suffix <- args[4]
aliTool <- args[5]

pwyrea
reaec
seqfiles
setnames(seqfiles,"ecs","ec")

# identify reactions probably catalysed by protein complexes
cplx <- merge(seqfiles, seq_headers, by = "file", allow.cartesian = TRUE)


#
rxndt <- merge(pwyrea, reaec, by = c("rea", "reaName", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, seqfiles[use == TRUE, .(rea, reaName, ec, file, type, src)], by = c("rea","reaName", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, alignments, by = "file", all.x = TRUE, allow.cartesian = TRUE)
