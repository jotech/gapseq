library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
# TODO check if all required arguments are provided:
# 1. path to temporary tsv file with all pathways to be tested
# 2.

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}


pwyDB <- fread(args[1])
setkey(pwyDB, "V1")

saveRDS(pwyDB, "~/tmp/pwyDB.RDS")

pwyrea <- lapply(pwyDB$V1, FUN = function(pwyi) {
  tmpdat <- pwyDB[pwyi]
  keyrea <- strsplit(tmpdat[,V8],",")[[1]]
  reaids <- tmpdat[, str_split(V6,",")][[1]]
  namestmp <- tmpdat[, V9]
  namestmp <- gsub("\\(([^()]*?);([^()]*?)\\)", "\\(\\1,\\2\\)",namestmp) # In rare cases the reactio name separtor ";" occurs within a reaction name, when inside brackets "(...)". These semicolons need to be replaced before string splitting.
  namestmp <- gsub("\\&[a-zA-Z0-9]+;","___",namestmp) # in other rare cases, reaction names contain HTML character entities such as "&pi;", causing issues when splitting at ";"
  reanames <- str_split(namestmp,";")[[1]]
  ecs <- tmpdat[, str_split(V7,",")][[1]]
  if(length(reanames) != length(reaids)) {
    message(paste0("Mismatch in the number of reaction IDs and reaction names in pathway ",pwyi,". Replacing reaction names with arbitrary rxn[1-n]..."))
    reanames <- paste0("rxn",1:length(reaids))
  }

  data.table(rea = reaids,
             reaName = reanames,
             ec = ecs,
             keyrea = reaids %in% keyrea
  )
})
names(pwyrea) <- pwyDB$V1
pwyrea <- rbindlist(pwyrea, idcol = "pwyID")

reaec <- pwyrea[, getDBhit(rea, reaName, ec, "seed"), by = .(rea, reaName, ec)]


