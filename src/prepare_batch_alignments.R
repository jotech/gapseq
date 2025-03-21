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
pwyDB <- fread(args[1], sep="\t", quote=FALSE)
setkey(pwyDB, "V1")
database <- args[2]
taxonomy <- args[3]
seqSrc <- args[4]
force_offline <- args[5] == "true"
update_manually <- args[6] == "true"

# ensure correct column types
pwyDB$V8 <- as.character(pwyDB$V8) # key reactions
pwyDB$V14 <- as.character(pwyDB$V14) # spontaneous reactions

# for debugging
saveRDS(pwyDB, "~/tmp/pwyDB.RDS")
# pwyDB <- readRDS("~/tmp/pwyDB.RDS")

pwyrea <- lapply(pwyDB$V1, FUN = function(pwyi) {
  tmpdat <- pwyDB[pwyi]
  keyrea <- strsplit(tmpdat[,V8],",")[[1]]
  spontrea <- strsplit(tmpdat[,V14],",")[[1]]
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
             keyrea = reaids %in% keyrea,
             spont = reaids %in% spontrea
  )
})
names(pwyrea) <- pwyDB$V1
pwyrea <- rbindlist(pwyrea, idcol = "pwyID")

# get target database hits
source(paste0(script.dir,"/getDBhit.R"))
reaec <- pwyrea[spont == FALSE, getDBhit(rea, reaName, ec, database), by = .(rea, reaName, ec)] # spont == FALSE because we don't need to download sequences for spontaneous reactions

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

  fileDT <- data.table(rea = reaID, reaName = reaName, ecs = ecs, altecs = altecs,
                       file = c(sq_rxn,
                                sq_revEC, sq_revName,
                                sq_unrevEC, sq_unrevName,
                                sq_user),
                       type = c("metacyc",
                                rep("EC",length(sq_revEC)), rep("reaName",length(sq_revName)),
                                rep("EC",length(sq_revEC)), rep("reaName",length(sq_revName)),
                                rep("EC",length(allecs)),"reaName","metacyc"))
  fileDT[, src := sub("/.*$","",file)]
  fileDT[, use := FALSE] # tbd later, dep. on "seqSrc"

  # If only unreviewed or only reviewed, we can directly remove some file entries
  if(seqSrc == "1") {
    # only reviewed
    fileDT <- fileDT[src %in% c("rxn","rev","user")]
  } else if(seqSrc == "4") {
    # only unreviewed
    fileDT <- fileDT[src %in% c("rxn","unrev","user")]
  }

  # if the reaction ID is not listed in metaGenes, we can remove the "rxn" entry
  if(reaID %notin% metaGenes$rxn) {
    fileDT <- fileDT[src != "rxn"]
  }

  fileDT[, fex := file.exists(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]
  fileDT <- fileDT[!(src == "user" & fex == FALSE)] # when user fasta does not exist, nothing to collect/use here
  fileDT[fex == TRUE, file_size := file.size(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]

  if(!force_offline) {
    # TODO: Check if something needs to be downloaded here
    # Everything of type "EC" should be there
    if(fileDT[type == "EC" & !fex, .N] > 0) {
      print(reaID)
      print(fileDT[type == "EC" & !fex])
    }
    # files for hashed reaName should be there in cases where all EC files are empty
    # files for "rxn" should be there, if the reaction ID is listed in table metaGenes
  }

  return(fileDT)
}

seqfiles <- lapply(1:nrow(reaec), function(i) {
  return(identifySeqFiles(reaID = reaec[i,rea],
                          reaName = reaec[i,reaName],
                          ecs = reaec[i,ec],
                          altecs = reaec[i,altec]))
})
seqfiles <- rbindlist(seqfiles)

#-------------------------------------------------------------------------------
# (n) Download sequences (if required)
#-------------------------------------------------------------------------------

seqfiles$uniprot_query <- NA_character_
if(!force_offline) {
  seqfiles[type == "EC" & (!fex | update_manually) & src == "rev",
           uniprot_query := paste0(script.dir,"/uniprot.sh -e \"",gsub("^rev/|\\.fasta$","",file),"\" -t \"",taxonomy,"\" -i 0.9 -o")]
  seqfiles[type == "EC" & (!fex | update_manually) & src == "unrev",
           uniprot_query := paste0(script.dir,"/uniprot.sh -u -e \"",gsub("^unrev/|\\.fasta$","",file),"\" -t \"",taxonomy,"\" -i 0.5 -o")]

  # perform download for ECs
  seqfiles_dlEC <- seqfiles[!is.na(uniprot_query)][!duplicated(file)]
  ndl <- nrow(seqfiles_dlEC)
  if(ndl > 0) {
    cat(ndl, "sequence files need to be downloaded from UniProt (via EC / metacyc-genes):\n")
    for(i in 1:ndl) {
      cat("\r", i, "/", ndl,"(",seqfiles_dlEC[i,file],")\n")
      system(seqfiles_dlEC[i,uniprot_query], ignore.stdout = TRUE, ignore.stderr = FALSE)
    }
    cat("\n")
  }

  # update file existence and size columns
  seqfiles[type == "EC" & !is.na(uniprot_query), fex := file.exists(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]
  seqfiles[type == "EC" & !is.na(file_size) & fex == TRUE, file_size := file.size(paste0(script.dir,"/../dat/seq/",taxonomy,"/",file))]

  # use reaction names for uniprot queries if no EC for reaction, or EC fastas are empty
  # seqfiles[, use_reanames := all(ecs == "") || sum(file_size * (type == "EC" & fex), na.rm = TRUE) == 0,by = .(rea, reaName, ecs)]

  # perform sequence download for direct uniprot links from metacyc reaction IDs
  mc_uplinks <- merge(seqfiles[src=="rxn" & (!fex | update_manually),
                               .(rea, file)],metaGenes[,.(rea = rxn, uniprot)], by = "rea")
  updlids <- unique(mc_uplinks$uniprot)
  ndl <- length(updlids)
  if(ndl > 0) {
    cat(ndl, "files need to be downloaded from UniProt via diect metacyc gene links:\n")
    for(i in 1:ndl) {
      cat("\r",i,"/",ndl)
      system(paste0(
        script.dir,"/uniprot.sh -d ",updlids[i]," -t ",taxonomy," -i 0.9 -o"
      ), ignore.stdout = TRUE, ignore.stderr = FALSE)
    }
    cat("\n")

    for(reai in unique(mc_uplinks$file)) {
      file.create(paste0(script.dir,"/../dat/seq/",taxonomy,"/",reai))
      file.append(paste0(script.dir,"/../dat/seq/",taxonomy,"/",reai),
                  paste0(script.dir,"/../dat/seq/",taxonomy,"/rxn/",mc_uplinks[file == reai,uniprot],".fasta"))
    }
  }


}
