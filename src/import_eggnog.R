# Import eggnog-mapper annotations and create gapseq-find like output files
#
# Arguments:
# `file_annotations` - eggnog-mapper annotation output file
# `file_hits` - eggnog-mapper annotation hits file (assumes diamond output currently)
#
library(data.table)

import_eggnog <- function(file_annotations, file_hits){
  eggnog_anno <- fread(cmd=paste("grep -v '^##'", file_annotations))
  eggnog_anno[EC!="-", seed:=getDBhit(EC, rea="aa", database="seed", reaName="")$dbhit, by=`#query`]
  diamond_header <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp")
  eggnog_hits <- fread(file_hits, col.names=diamond_header)
  eggnog <- merge(eggnog_anno, eggnog_hits, by.x="#query", by.y="qseqid", all.x = T)
  # old gapseq used protein as query and genome as subject, which is different in eggnog/diamond. therefore qseqid and sseqid are interchanged
  gs_eggnog <- eggnog[EC!="-", list(rxn=`#query`, name=Preferred_name, ec=EC, bihit=NA, qseqid=sseqid, pident, 
                      evalue=evalue.x, bitscore, qcovs=scovhsp, stitle=`#query`, sstart=qstart, send=qend, pathway=NA, status="good_blast", 
                      pathway.status=NA, dbhit=seed, complex=NA, exception=0, complex.status=NA)]
  return(gs_eggnog)
}

# make sure script path is defined
if(!exists("script.dir"))
  stop("Script path not defined for 'getDBhit()' function")

# load dependencies
source(paste0(script.dir, "/getDBhit.R"))
