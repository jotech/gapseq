library(data.table)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
# TODO check if all required arguments are provided:
# 1. path to temporary tsv file with all pathways to be tested
# 2. target reaction database (vmh/seed) ("-d" in gapseq find)
# 3. Taxonomic range for reference sequences
# 4. Reference sequence type to use for alignments

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

# get arguments
pwyDB <- fread(args[1])
setkey(pwyDB, "V1")
database <- args[2]
taxonomy <- args[3]
seqSrc <- args[4]

saveRDS(pwyDB, "~/tmp/pwyDB.RDS")
# pwyDB <- readRDS("~/tmp/pwyDB.RDS")

pwyrea <- lapply(pwyDB$V1, FUN = function(pwyi) {
  tmpdat <- pwyDB[pwyi]
  keyrea <- strsplit(tmpdat[,V8],",")[[1]]
  reaids <- tmpdat[, str_split(V6,",")][[1]]
  namestmp <- tmpdat[, V9]
  namestmp <- gsub("\\(([^()]*?);([^()]*?)\\)", "\\(\\1,\\2\\)",namestmp) # In rare cases the reaction name separator ";" occurs within a reaction name, when inside brackets "(...)". These semicolons need to be replaced before string splitting.
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

# get target database hits
source(paste0(script.dir,"/getDBhit.R"))
reaec <- pwyrea[, getDBhit(rea, reaName, ec, database), by = .(rea, reaName, ec)]

#-------------------------------------------------------------------------------
# collate reference sequences to use for later alignments
#-------------------------------------------------------------------------------

# allseqfiles <- dir(paste0(script.dir,"/../dat/seq/",taxonomy,"/"),
#                    recursive = TRUE,
#                    pattern = "\\.fasta$")
# fsz <- file.size(paste0(script.dir,"/../dat/seq/",taxonomy,"/",allseqfiles))
# names(fsz) <- allseqfiles

# for some metacyc reactions, there are already uniprot accessions assigned.
metaGenes <- fread(paste0(script.dir,"/../dat/meta_genes.csv"))
metaGenes[, rxn := gsub("^\\||\\|$","",rxn)]
setkey(metaGenes, "rxn")
metaGenes <- metaGenes[uniprot != ""]

md5sum <- function(str) {
  temp_file <- tempfile()
  writeLines(str, temp_file, useBytes = TRUE)
  hash <- tools::md5sum(temp_file)
  unlink(temp_file)  # Remove the temporary file
  return(unname(hash))
}

identifySeqFiles <- function(reaID, reaName, ecs, altecs) {
  allecs <- c(unlist(str_split(ecs, "/")),
              unlist(str_split(altecs, "/")))
  allecs <- allecs[grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$",allecs)]
  rnmd5 <- md5sum(reaName)

  # MetaCyc geneXreaction assignments
  sq_rxn <- paste0("rxn/",reaID,".fasta") # what when this one is not there?

  # geneXec assignments (reviewed and unreviewed)
  if(length(allecs) > 0) {
    sq_revEC <- paste0("rev/",allecs,".fasta")
    sq_unrevEC <- paste0("unrev/",allecs,".fasta")
  } else {
    sq_revEC <- character(0L)
    sq_unrevEC <- character(0L)
  }

  # geneXreaName assignments (reviewed and unreviewed)
  if(reaName != "") {
    sq_revName <- paste0("rev/",rnmd5,".fasta")
    sq_unrevName <- paste0("unrev/",rnmd5,".fasta")
  } else {
    sq_revName <- character(0L)
    sq_unrevName <- character(0L)
  }

  # user files
  sq_user <- paste0("user/",
                    c(allecs,rnmd5,reaID),
                    ".fasta")

  fileDT <- data.table(file = c(sq_rxn,
                                sq_revEC, sq_revName,
                                sq_unrevEC, sq_unrevName,
                                sq_user),
                       type = c("metacyc",
                                rep("EC",length(sq_revEC)), rep("reaName",length(sq_revName)),
                                rep("EC",length(sq_revEC)), rep("reaName",length(sq_revName)),
                                rep("EC",length(allecs)),"reaName","metacyc"))
  fileDT[, src := sub("/.*$","",file)]
  fileDT[, use := FALSE]
  fileDT[, fex := file.exists(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]
  fileDT <- fileDT[!(src == "user" & fex == FALSE)] # when user fasta does not exist, nothing to collect/use here
  fileDT[fex == TRUE, file_size := file.size(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]

  # TODO: Check if something needs to be downloaded here
  # Everything of type "EC" should be there
  # files for hashed reaName should be there in cases where all EC files are empty
  # files for "rxn" should be there, if the reaction ID is listed in table metaGenes

  if(seqSrc == "1") {
    # only reviewed

  } else if(seqSrc == "2") {
    # unreviewed only if reviewed not available

  } else if(seqSrc == "3") {
    # reviewed+unreviewed

  } else if(seqSrc == "4") {
    # only unreviewed

  }

}
