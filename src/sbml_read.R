#!/usr/bin/Rscript

library(sybilSBML, quietly=T)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
    stop("Needs sbml file ...")
}

mod <- readSBMLmod(args[1])
write.table(mod@react_id, file="reactions.csv", row.names=F, col.names=F, quote=F)
saveRDS(mod, "./model.RDS")
