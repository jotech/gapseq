suppressMessages(library(data.table))

#-------------------------------------------------------------------------------
# (0) Read pre-alignment data and alignment stats
#-------------------------------------------------------------------------------
load("~/tmp/preblast_data.RData")
alignments <- fread("~/tmp/query.blast")
setnames(alignments, c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send"))
alignments[, file := sub("\\|.+$","",qseqid)]
alignments[, qseqid := sub("^.+\\.fasta\\|","",qseqid)]

#-------------------------------------------------------------------------------
# (1) Parse arguments and gapseq/script path
#-------------------------------------------------------------------------------

pwyrea
reaec
seqfiles
setnames(seqfiles,"ecs","ec")

#
rxndt <- merge(pwyrea, reaec, by = c("rea", "reaName", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, seqfiles[use == TRUE, .(rea, reaName, ec, file, type, src)], by = c("rea","reaName", "ec"), all.x = TRUE)
rxndt <- merge(rxndt, alignments, by = "file", all.x = TRUE, allow.cartesian = TRUE)
