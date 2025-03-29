suppressMessages(library(data.table))

#-------------------------------------------------------------------------------
# (0) Read pre-alignment data and alignment stats
#-------------------------------------------------------------------------------
load("prealignment_data.RData") # load("~/tmp/prealignment_data.RData")
setnames(seqfiles,"ecs","ec")

alicols <- c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","sseq")
alignments <- fread("alignments.tsv") # alignments <- fread("~/tmp/alignments.tsv")
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
identcutoff_exception <- as.numeric(args[4])
subunit_cutoff <- as.numeric(args[5])/100
completenessCutoffNoHints <- as.numeric(args[6])/100
completenessCutoff <- as.numeric(args[7])/100

#-------------------------------------------------------------------------------
# (2) Read additional required gapseq files
#-------------------------------------------------------------------------------
source(paste0(script.dir,"/complex_detection2.R"))
dt_exceptions <- fread(paste0(script.dir,"/../dat/exception.tbl"), skip = "enzyme/reaction")

#-------------------------------------------------------------------------------
# (3) Complex detection
#-------------------------------------------------------------------------------

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


#-------------------------------------------------------------------------------
# (4) Inferring Reaction presence/absence
#-------------------------------------------------------------------------------

rxndt <- merge(reaec[,.(rxn = rea, name = reaName, ec, dbhit, spont)],
               seqfiles[use == TRUE, .(rxn = rea, name = reaName, ec, file, type, src)], by = c("rxn","name", "ec"),
               all.x = TRUE)
rxndt <- merge(rxndt, cplx_bin[,.(rxn = rea, name = reaName, ec, is_complex, subunit_count, subunits)], by = c("rxn","name", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, alignments, by = "file", all.x = TRUE, allow.cartesian = TRUE)
rxndt$complex <- cplx[rxndt[,.(rxn, name, ec, qseqid)], complex]
rxndt[is.na(is_complex), is_complex := FALSE]

exc_ec <- dt_exceptions[grepl("^[0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+$", `enzyme/reaction`),
                        paste0("(",`enzyme/reaction`,"($|/))", collapse = "|")]
rxndt[, exception := as.numeric(grepl(exc_ec, ec))]

rxndt[bitscore >= bitcutoff & pident >= identcutoff, status := "good_blast"]
rxndt[bitscore < bitcutoff | pident < identcutoff, status := "bad_blast"]
rxndt[exception == 1 & !is.na(bitscore) & pident < identcutoff_exception, status := "bad_blast"]
rxndt[is.na(qseqid) & !is.na(file), status := "no_blast"]
rxndt[is.na(qseqid) & is.na(file), status := "no_seq_data"]
rxndt[is.na(qseqid) & is.na(file) & spont == TRUE, status := "spontaneous"]

# infere if complexes are complete (enough)
rxndt[is_complex == TRUE, subunits_found := length(unique(complex[!is.na(complex) & complex != "Subunit undefined" & status == "good_blast"])), by = .(rxn,name,ec)]
rxndt[is_complex == TRUE, subunit_undefined_found := any(complex == "Subunit undefined" & status == "good_blast"), by = .(rxn,name,ec)]
rxndt[(subunits_found / subunit_count > subunit_cutoff) | (subunits_found / subunit_count == subunit_cutoff & subunit_undefined_found == TRUE), complex.status := 1]

# remove duplicate hits (same reaction & complex & gene in target genome)
rxndt <- rxndt[order(rxn, name, ec, stitle, complex, -bitscore)]
rxndt <- unique(rxndt, by = c("rxn", "name", "ec", "stitle", "complex"))

# merge in pathway info
rxndt <- merge(pwyrea[, .(rxn = rea, name = reaName, ec, pathway = pwyID, keyrea)],
               rxndt, by = c("rxn","name", "ec"),
               allow.cartesian = TRUE)
rxndt <- rxndt[order(pathway, rxn, complex, -bitscore)]

#-------------------------------------------------------------------------------
# (4) Inferring Pathway presence/absence
#-------------------------------------------------------------------------------
pwydt <- unique(rxndt, by = c("pathway","rxn"))
pwydt <- pwydt[, .(pathway, rxn, spont, keyrea, is_complex, complex.status, status)]
pwydt <- pwydt[, .(NrReaction = .N,
                   NrSpontaneous = sum(spont),
                   NrVague = sum(status == "no_seq_data"),
                   NrKeyReaction = sum(keyrea),
                   NrReactionFound = sum((is_complex == FALSE & status == "good_blast") | (is_complex == TRUE & !is.na(complex.status))),
                   NrKeyReactionFound = sum(status == "good_blast" & keyrea == TRUE),
                   ReactionsFound = paste(rxn[(is_complex == FALSE & status == "good_blast") | (is_complex == TRUE & !is.na(complex.status))], collapse = " "),
                   SpontaneousReactions = paste(rxn[spont == TRUE], collapse = " "),
                   KeyReactions = paste(rxn[keyrea == TRUE], collapse = " ")),by = pathway]
pwydt[, Completeness := NrReactionFound / (NrReaction - NrVague - NrSpontaneous) * 100]
pwydt[, Prediction := Completeness >= completenessCutoffNoHints]

#-------------------------------------------------------------------------------
# (n) Export reaction and pathway tables
#-------------------------------------------------------------------------------
fwrite(rxndt, file = "output.tbl", sep = "\t", quote = FALSE)
fwrite(rxndt, file = "output_pwy.tbl", sep = "\t", quote = FALSE)
