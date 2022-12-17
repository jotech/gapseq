#!/usr/bin/env Rscript

library(data.table)
library(stringr)

# first argument: scriptdir
# second argument: hmmsearch result table (tblout)
args <- commandArgs(trailingOnly=TRUE)


hmm.markers <- fread(paste0(args[1], "/../dat/seq/hmm/gram_markers.tsv"))
hmmwhy <- readLines(args[2]) # why this table format, hmmer?
hmmwhy <- hmmwhy[!grepl("^#", hmmwhy)] # # remove header/footer/comment lines
hmmwhy.l <- unlist(lapply(str_split(hmmwhy," +", n = 8),
                          function(x) paste(x[1:7], collapse = "\t")))

# make table parsable for fread
tmpf <- tempfile()
cat(hmmwhy.l, sep = "\n", file = tmpf)
dt <- fread(tmpf, header = FALSE)
frm <- file.remove(tmpf)
setnames(dt, c("target.name","target.acc","query.name","query.acc","eval",
               "bitscore","bias"))

# remove suboptimal hits
dt <- dt[order(query.acc,eval)]
dt <- dt[!duplicated(query.acc)]

# score domains
dt <- merge(dt, hmm.markers, all = TRUE, by = "query.acc")
dt[is.na(bitscore) | bitscore < 50, bitscore := 0]
dt[, contr := sqrt(bitscore)]

gramres <- dt[, sum(contr)/.N, by = Gram.staining][order(-V1)][1, Gram.staining]

cat(gramres)

q(save = "no")

