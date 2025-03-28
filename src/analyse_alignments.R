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

# identify reactions probably catalyzed by protein complexes
cplx <- merge(seqfiles[use == TRUE, .(file, rea, reaName, ec)], seq_headers, by = "file", allow.cartesian = TRUE)
cplx <- unique(cplx, by = c("rea","reaName","ec", "header"))
cplx[, is_complex := any(grepl("subunit|chain|polypeptide|component",header)), by = .(rea, reaName, ec)]
cplx[is_complex == TRUE, complex := complex_detection(header), by = .(rea, reaName, ec)]
cplx[, subunit_count := length(unique(complex[!is.na(complex)])), by = .(rea, reaName, ec)]
cplx <- cplx[subunit_count < 2, is_complex := FALSE]
cplx[, qseqid := sub(" .*$","",header)]
cplx[is_complex == TRUE & is.na(complex), complex := "Subunit undefined"]
setkeyv(cplx, cols = c("rea","reaName","ec","qseqid"))
cplx <- unique(cplx, by =  c("rea","reaName","ec","qseqid"))
cplx_bin <- cplx[,.(is_complex = any(is_complex),
                    subunit_count = subunit_count[1],
                    subunits = paste(unique(complex),collapse = ",")), by = .(rea, reaName, ec)]
cplx_bin[is_complex == FALSE, subunit_count := NA_integer_]

#
rxndt <- merge(reaec[,.(rxn = rea, name = reaName, ec, dbhit)],
               seqfiles[use == TRUE, .(rxn = rea, name = reaName, ec, file, type, src)], by = c("rxn","name", "ec"),
               all.x = TRUE)
rxndt <- merge(rxndt, cplx_bin[,.(rxn = rea, name = reaName, ec, is_complex, subunit_count)], by = c("rxn","name", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, alignments, by = "file", all.x = TRUE, allow.cartesian = TRUE)
rxndt$complex <- cplx[rxndt[,.(rxn, name, ec, qseqid)], complex]
rxndt[is.na(is_complex), is_complex := FALSE]


rxndt[!is.na(file) & is_complex == TRUE, complex := cplx[.(rxn, name, ec, qseqid),]]


# rxndt <- merge(pwyrea, reaec, by = c("rea", "reaName", "ec"), all.x = TRUE)
# rxndt <- merge(rxndt, seqfiles[use == TRUE, .(rea, reaName, ec, file, type, src)], by = c("rea","reaName", "ec"), all.x = TRUE)
# rxndt <- merge(rxndt, alignments, by = "file", all.x = TRUE, allow.cartesian = TRUE)
