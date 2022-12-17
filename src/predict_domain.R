#!/usr/bin/env Rscript

library(data.table)
library(stringr)

# first argument: scriptdir
# second argument: hmmsearch result table (tblout)
args <- commandArgs(trailingOnly=TRUE)


hmm.markers <- fread(paste0(args[1], "/../dat/seq/hmm/domain_markers.tsv"))
hmmwhy <- readLines(args[2]) # why this table format, hmmer?
hmmwhy <- hmmwhy[!grepl("^#", hmmwhy)]# removing header and footer lines
hmmwhy.l <- unlist(lapply(str_split(hmmwhy," +", n = 8),
                          function(x) paste(x[1:7], collapse = "\t")))


# make table parsable for fread
tmpf <- tempfile()
cat(hmmwhy.l, sep = "\n", file = tmpf)
dt <- fread(tmpf, header = FALSE, fill = TRUE)
frm <- file.remove(tmpf)
setnames(dt, c("target.name","target.acc","query.name","query.acc","eval",
               "bitscore","bias"))

# remove suboptimal hits
dt <- dt[order(query.acc,eval)]
dt <- dt[!duplicated(query.acc)]

# score domains
dt[, query.acc := gsub("\\.[0-9]$","", query.acc)]
dt <- merge(dt, hmm.markers, all = TRUE, by = "query.acc")
dt[is.na(bitscore) | bitscore < 50, bitscore := 0]
dt[, contr := sqrt(bitscore)]

taxres <- dt[, sum(contr)/.N, by = domain][order(-V1)][1, domain]

cat(taxres)

q(save = "no")

