suppressMessages(library(data.table))
library(parallel)
library(stringr)

#-------------------------------------------------------------------------------
# (0) Read pre-alignment data and alignment stats
#-------------------------------------------------------------------------------
load("prealignment_data.RData") # load("~/tmp/prealignment_data.RData")
setnames(seqfiles,"ecs","ec")

alicols <- c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","sseq")
if(file.size("alignments.tsv") == 0) {
  alignments <- data.table(matrix(0,nrow = 1, ncol = length(alicols)-1))
  alignments[, V1 := as.character(V1)]
  alignments[, V6 := as.character(V6)]
} else {
  alignments <- fread("alignments.tsv") # alignments <- fread("~/tmp/alignments.tsv")
}
setnames(alignments, alicols[1:ncol(alignments)])
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
n_threads <- as.integer(args[8])
vagueCutoff <- as.numeric(args[9])
verbose <- as.integer(args[10])
pwyDB <- fread(args[11], sep="\t", quote=FALSE)

#-------------------------------------------------------------------------------
# (2) Read additional required gapseq files
#-------------------------------------------------------------------------------
source(paste0(script.dir,"/complex_detection.R"))
dt_exceptions <- fread(paste0(script.dir,"/../dat/exception.tbl"), skip = "enzyme/reaction")

#-------------------------------------------------------------------------------
# (3) Complex detection
#-------------------------------------------------------------------------------

# identify reactions probably catalyzed by protein complexes
cplx <- merge(seqfiles[use == TRUE, .(file, rea, reaName, ec)], seq_headers, by = "file", allow.cartesian = TRUE)
cplx <- unique(cplx, by = c("rea","reaName","ec", "header"))
cplx[, is_complex := any(grepl("subunit|chain|polypeptide|component",header)), by = .(rea, reaName, ec)]

if(any(cplx$is_complex == TRUE)) {
  cplx_sub <- cplx[is_complex == TRUE]
  split_groups <- split(cplx_sub, by = c("rea", "reaName", "ec"))
  results_list <- mclapply(split_groups, function(dt_group) {
    dt_group[, complex := complex_detection(header)]
    return(dt_group)
  }, mc.cores = n_threads)
  cplx_result <- rbindlist(results_list)
  cplx[cplx_result, on = .(rea, reaName, ec, header), complex := i.complex]
} else {
  cplx[, complex := NA_character_]
}

#cplx[is_complex == TRUE, complex := complex_detection(header), by = .(rea, reaName, ec)]

cplx[, subunit_count := length(unique(complex[!is.na(complex)])), by = .(rea, reaName, ec)]
cplx[subunit_count < 2, is_complex := FALSE]
cplx[subunit_count < 2, complex := NA_character_]
cplx[, qseqid := sub(" .*$","",header)]
cplx[is_complex == TRUE & is.na(complex), complex := "Subunit undefined"]
setkeyv(cplx, cols = c("rea","reaName","ec","qseqid"))
cplx <- unique(cplx, by =  c("rea","reaName","ec","qseqid"))
cplx_bin <- cplx[,.(is_complex = any(is_complex),
                    subunit_count = subunit_count[1],
                    subunits = paste(unique(complex),collapse = ",")), by = .(rea, reaName, ec)]
cplx_bin[is_complex == FALSE, subunit_count := NA_integer_]
#print(cplx[, table(rea, complex, useNA = "always")])

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
                   NrKeyReaction = sum(keyrea & status != "no_seq_data"),
                   NrReactionFound = sum((is_complex == FALSE & status == "good_blast") | (is_complex == TRUE & !is.na(complex.status))),
                   NrKeyReactionFound = sum((is_complex == FALSE & status == "good_blast" & keyrea == TRUE) | (is_complex == TRUE & !is.na(complex.status) & keyrea == TRUE)),
                   ReactionsFound = paste(rxn[(is_complex == FALSE & status == "good_blast") | (is_complex == TRUE & !is.na(complex.status))], collapse = " "),
                   SpontaneousReactions = paste(rxn[spont == TRUE], collapse = " "),
                   KeyReactions = paste(rxn[keyrea == TRUE], collapse = " ")),by = pathway]

# infere completeness
pwydt[NrVague/(NrReaction-NrSpontaneous) < vagueCutoff, Completeness := NrReactionFound / (NrReaction - NrVague - NrSpontaneous) * 100] # if vague reaction are considered they should not influence the completeness threshold
pwydt[NrVague/(NrReaction-NrSpontaneous) >= vagueCutoff, Completeness := NrReactionFound / (NrReaction - NrSpontaneous) * 100]

# infere pathway presence/absence
pwydt[strictCandidates == FALSE, Prediction := Completeness >= completenessCutoffNoHints*100 & NrKeyReaction == NrKeyReactionFound]
pwydt[strictCandidates == FALSE & NrKeyReaction > 0 & NrKeyReaction == NrKeyReactionFound & Completeness >= completenessCutoff*100, Prediction := TRUE]
pwydt[NrReaction == NrVague + NrSpontaneous, Prediction := FALSE]

# add pathway status term (explaining why pathway is predicted to be present)
pwydt[, pathway.status := NA_character_]
pwydt[Completeness == 100, pathway.status := "full"]
pwydt[Prediction == TRUE & Completeness < 100 & NrKeyReactionFound == NrKeyReaction, pathway.status := "threshold"]
pwydt[Prediction == TRUE & Completeness < completenessCutoffNoHints*100, pathway.status := "keyenzyme"]

# add a column with pathway names
pwydt <- merge(pwydt, pwyDB[,.(pathway = V1, Name = V2)], by = "pathway")

# reorder columns
pwydt <- pwydt[, .(pathway, Name, Prediction, Completeness, pathway.status,
                   NrReaction, NrSpontaneous, NrVague, NrKeyReaction,
                   NrReactionFound, NrKeyReactionFound,
                   ReactionsFound, SpontaneousReactions, KeyReactions)]

#-------------------------------------------------------------------------------
# (4) Add pathway prediction info to reaction table, too
#-------------------------------------------------------------------------------

rxndt <- merge(rxndt, pwydt[, .(pathway, pathway.status)], by = "pathway")

#-------------------------------------------------------------------------------
# ORFs associated with at least one found reaction (used for coverage calc.)
#-------------------------------------------------------------------------------

nmapped <- rxndt[(is_complex == FALSE & status == "good_blast") | (is_complex == TRUE & !is.na(complex.status)),
                 length(unique(stitle))]
writeLines(as.character(nmapped), "nmappedORFs.tmp")

#-------------------------------------------------------------------------------
# (n) Export reaction and pathway tables
#-------------------------------------------------------------------------------
setnames(pwydt, c("pathway", "pathway.status"), c("ID", "Status"))

fwrite(rxndt, file = "output.tbl", sep = "\t", quote = FALSE)
fwrite(pwydt, file = "output_pwy.tbl", sep = "\t", quote = FALSE)
