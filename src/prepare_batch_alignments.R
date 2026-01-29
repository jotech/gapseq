suppressMessages(library(data.table))
library(stringr)
suppressMessages(library(Biostrings))
library(parallel)

#-------------------------------------------------------------------------------
# (0) Parse arguments and gapseq/script path
#-------------------------------------------------------------------------------

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
pwyDB[, V9 := gsub("  "," ",V9)] # remove double spaces in reaction names
setkey(pwyDB, "V1")
database <- args[2]
taxonomy <- args[3]
seqSrc <- args[4]
force_offline <- args[5] == "true"
update_manually <- args[6] == "true"
use_gene_seq <- args[7] == "true"
n_threads <- as.integer(args[8])
verbose <- as.integer(args[9])
onlyList <- args[10] == "true"
seqdb <- args[11]

# some hard-coded correction of syntax errors in reaction names that would break string parsing (remove when meta_pwy.tbl is corrected/updated)
pwyDB[V1 %in% c("|PWY18C3-10|","|PWY18C3-12|"), V9 := sub("sucrose-3-isobutanoyl-4-isovaleryl-3;-isovaleroyltransferase",
                                                        "sucrose-3-isobutanoyl-4-isovaleryl-3'-isovaleroyltransferase",
                                                        V9)]
pwyDB[V1 == "|PWY-6024|", V9 := sub("vitexin 2;;-O-arabinoside 7-O-galactosyltransferase",
                                    "vitexin 2''-O-arabinoside 7-O-galactosyltransferase",
                                    V9)]

# some reaction names contain HTML character entities such as "&pi;". In those cases, drop the "&" and ";"
pwyDB[, V9 := gsub("&([a-zA-Z0-9#]+);", "\\1", V9)]

# ensure correct column types
pwyDB$V8 <- as.character(pwyDB$V8) # key reactions
pwyDB$V14 <- as.character(pwyDB$V14) # spontaneous reactions

# when ec numbers are searched -> empty reaction name
pwyDB[is.na(V9), V9 := ""]

# for debugging
# saveRDS(pwyDB, "~/tmp/pwyDB.RDS")
# pwyDB <- readRDS("~/tmp/pwyDB.RDS")

if(onlyList) {
  pwyDB[, tmpstr := paste0(gsub("\\|","",V1), " – ",
                           V2, " with ",
                           str_count(V6,",")+1, " reactions\n"), by = V1]
  parseCats <- function(str) {
    out <- unlist(str_split(str,","))
    out <- gsub("\\|","",out)
    out <- out[out %notin% c("THINGS","FRAMES","Pathways","Generalized-Reactions")]
    out <- paste0("(",paste(out,collapse = ","),")")
    return(out)
  }
  pwyDB[, cats := parseCats(V4), by = V1]
  pwyDB[, tmpstr := paste0(1:.N,"/",.N," ",tmpstr,cats,"\n\n")]
  cat("\n")
  cat(pwyDB$tmpstr, sep = "")
  cat("\n")
  quit(save = "no")
}

#-------------------------------------------------------------------------------
# (1) Create pathway table
# keys: `pwyID`, `rea`
#-------------------------------------------------------------------------------

pwyrea <- mclapply(pwyDB$V1, mc.cores = n_threads, FUN = function(pwyi) {
  tmpdat <- pwyDB[pwyi]
  keyrea <- strsplit(tmpdat[,V8],",")[[1]]
  spontrea <- strsplit(tmpdat[,V14],",")[[1]]
  reaids <- tmpdat[, str_split(V6,",")][[1]]
  namestmp <- tmpdat[, V9]
  namestmp <- gsub("\\(([^()]*?);([^()]*?)\\)", "\\(\\1,\\2\\)",namestmp) # In rare cases the reaction name separator ";" occurs within a reaction name, when inside brackets "(...)". These semicolons need to be replaced before string splitting.
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

#-------------------------------------------------------------------------------
# (2) Create reaction data base hits (target database: gapseq/seed reactions)
# keys: `rea`, `reaName`, `ec`)
#-------------------------------------------------------------------------------

# get target database hits (parallel comp?)
source(paste0(script.dir,"/getDBhit.R"))
reaec <- unique(pwyrea[,.(rea, reaName, ec, spont)])
reaec <- cbind(reaec,
               reaec[, getDBhit(rea, reaName, ec, database, n_threads)])

#-------------------------------------------------------------------------------
# (3) identify fasta files with reference sequences to use for later alignments
#-------------------------------------------------------------------------------

# for some metacyc reactions, there are already uniprot accessions assigned.
metaGenes <- fread(paste0(script.dir,"/../dat/meta_genes.csv"))
metaGenes[, rxn := gsub("^\\||\\|$","",rxn)]
metaGenes <- metaGenes[uniprot != ""]
metaGenes <- metaGenes[, .(uniprot = unlist(strsplit(uniprot,";"))), by = .(rxn,genes,pwy,ncbi,location,go,org)]
metaGenes <- metaGenes[uniprot != "Q91194"] # This is a wrong uniprot ID in the metaGenes table. THe correct one is "Q9I194" and it's also listed.
setkey(metaGenes, "rxn")

# md5sums of reaction names
rnuniq <- unique(reaec$reaName)
md5_hashes <- vapply(
  rnuniq,
  function(x) tools::md5sum(bytes = charToRaw(x)),
  character(1)
)
names(md5_hashes) <- rnuniq
rm(rnuniq)

identifySeqFiles <- function(reaID, reaName, ecs, altecs, spont) {
  allecs <- c(unlist(str_split(ecs, "/")),
              unlist(str_split(altecs, "/")))
  allecs <- allecs[grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$",allecs)]
  #print(reaName)
  rnmd5 <- md5_hashes[reaName]

  # MetaCyc geneXreaction assignments
  sq_rxn <- paste0("rxn/",reaID,".fasta")

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

  # If only unreviewed or only reviewed, we can directly drop some file entries
  if(seqSrc == "1") {
    # only reviewed
    fileDT <- fileDT[src %in% c("rxn","rev","user")]
  } else if(seqSrc == "4") {
    # only unreviewed
    fileDT <- fileDT[src %in% c("rxn","unrev","user")]
  }

  # if no direct metacyc gene links should be used ('rxn'), drop those entries
  if(!use_gene_seq)
    fileDT <- fileDT[type != "metacyc"]

  # reaName fasta (md5hash file names) should not be used if reaction names are arbitrary (e.g., rxn1)
  fileDT <- fileDT[!(type == "reaName" & grepl("^rxn[0-9]+$|^Rxn[0-9]+$|^RXN",reaName))]

  # if the reaction ID is not listed in metaGenes, we can remove the "rxn" entry
  if(use_gene_seq && reaID %notin% metaGenes$rxn) {
    fileDT <- fileDT[src != "rxn"]
  }

  # set dir of sequence database
  fileDT[, sdir := ifelse(src=="user",paste0(script.dir,"/../dat/seq"),seqdb)]

  fileDT[, fex := file.exists(paste0(sdir,"/",taxonomy,"/",file))]
  fileDT <- fileDT[!(src == "user" & fex == FALSE)] # when user fasta does not exist, nothing to collect/use here
  fileDT[fex == TRUE, file_size := file.size(paste0(sdir,"/",taxonomy,"/",file))]

  return(fileDT)
}

reaec_nospont <- reaec[spont == FALSE]
seqfiles <- mclapply(1:nrow(reaec_nospont), mc.cores = n_threads,
                     FUN = function(i) {
                       return(identifySeqFiles(reaID = reaec_nospont[i,rea],
                                               reaName = reaec_nospont[i,reaName],
                                               ecs = reaec_nospont[i,ec],
                                               altecs = reaec_nospont[i,altec]))
                     })
seqfiles <- rbindlist(seqfiles)
#print(seqfiles[sample(1:.N, 50)])
rm(reaec_nospont)

#-------------------------------------------------------------------------------
# (4) Download sequences (if required/wanted)
#-------------------------------------------------------------------------------

seqfiles$uniprot_query <- NA_character_
if(!force_offline) {
  seqfiles[type == "EC" & (!fex | update_manually) & src == "rev",
           uniprot_query := paste0(script.dir,"/uniprot.sh -e \"",gsub("^rev/|\\.fasta$","",file),"\" -t \"",taxonomy,"\" -i 0.9 -o -D ",seqdb)]
  seqfiles[type == "EC" & (!fex | update_manually) & src == "unrev",
           uniprot_query := paste0(script.dir,"/uniprot.sh -u -e \"",gsub("^unrev/|\\.fasta$","",file),"\" -t \"",taxonomy,"\" -i 0.5 -o -D ",seqdb)]

  # perform download for ECs
  seqfiles_dlEC <- seqfiles[!is.na(uniprot_query)][!duplicated(file)]
  ndl <- nrow(seqfiles_dlEC)
  if(ndl > 0) {
    cat(ndl, "sequence files need to be downloaded from UniProt (via EC / metacyc-genes):\n")
    for(i in 1:ndl) {
      cat(" ", i, "/", ndl,"(",seqfiles_dlEC[i,file],")\n")
      system(seqfiles_dlEC[i,uniprot_query], ignore.stdout = TRUE, ignore.stderr = FALSE)
    }
    cat("\n")
  }

  # update file existence and size columns
  seqfiles[type == "EC" & !is.na(uniprot_query), fex := file.exists(paste0(sdir,"/",taxonomy,"/",file))]
  seqfiles[type == "EC" & !is.na(uniprot_query) & fex == TRUE, file_size := file.size(paste0(sdir,"/",taxonomy,"/",file))]

  # use reaction names for uniprot queries if no valid EC for reaction, or EC fastas are empty
  seqfiles[, use_reanames := !any(type == "EC") || sum(file_size * (type == "EC" & fex), na.rm = TRUE) == 0,
           by = .(rea, reaName, ecs)]
  seqfiles[use_reanames == TRUE & type == "reaName" & (!fex | update_manually) & src == "rev" & reaName != "",
           uniprot_query := paste0(script.dir,"/uniprot.sh -r \"",gsub("\"","\\\\\"",reaName),"\" -t \"",taxonomy,"\" -i 0.9 -o -D ",seqdb)]
  seqfiles[use_reanames == TRUE & type == "reaName" & (!fex | update_manually) & src == "unrev" & reaName != "",
           uniprot_query := paste0(script.dir,"/uniprot.sh -u -r \"",gsub("\"","\\\\\"",reaName),"\" -t \"",taxonomy,"\" -i 0.5 -o -D ",seqdb)]
  seqfiles_dlRN <- seqfiles[type == "reaName" & !is.na(uniprot_query)][!duplicated(file)]
  ndl <- nrow(seqfiles_dlRN)
  if(ndl > 0) {
    cat(ndl, "sequence files need to be downloaded from UniProt (via reaction name):\n")
    for(i in 1:ndl) {
      cat(" ", i, "/", ndl,"(",seqfiles_dlRN[i,reaName],"/",seqfiles_dlRN[i,file],")\n")
      system(seqfiles_dlRN[i,uniprot_query], ignore.stdout = TRUE, ignore.stderr = FALSE)
    }
    cat("\n")
  }

  # update file existence and size columns
  seqfiles[type == "reaName" & !is.na(uniprot_query), fex := file.exists(paste0(sdir,"/",taxonomy,"/",file))]
  seqfiles[type == "reaName" & !is.na(uniprot_query) & fex == TRUE, file_size := file.size(paste0(sdir,"/",taxonomy,"/",file))]

  # perform sequence download for direct uniprot links from metacyc reaction IDs
  mc_uplinks <- merge(seqfiles[src=="rxn" & (!fex | update_manually),
                               .(rea, file)],metaGenes[,.(rea = rxn, uniprot)], by = "rea")
  updlids <- unique(mc_uplinks$uniprot)
  ndl <- length(updlids)
  if(ndl > 0) {
    cat(ndl, "files need to be downloaded from UniProt via direct metacyc gene links:\n")
    for(i in 1:ndl) {
      cat(" ",i,"/",ndl,"(",updlids[i],")\n")
      system(paste0(
        script.dir,"/uniprot.sh -d ",updlids[i]," -t ",taxonomy," -i 0.9 -o -D ",seqdb
      ), ignore.stdout = TRUE, ignore.stderr = FALSE)
    }
    cat("\n")

    for(reai in unique(mc_uplinks$file)) {
      file.create(paste0(seqdb,"/",taxonomy,"/",reai))
      file.append(paste0(seqdb,"/",taxonomy,"/",reai),
                  paste0(seqdb,"/",taxonomy,"/rxn/",mc_uplinks[file == reai,uniprot],".fasta"))
    }
  }
  seqfiles[, use_reanames := NULL] # column not needed anymore
}
seqfiles[, uniprot_query := NULL] # column not needed anymore

#-------------------------------------------------------------------------------
# (5) Decide which reference sequences files to use
#-------------------------------------------------------------------------------

# Decision rules:
# (a) sequences in folder "rxn" are always used
# (b) sequences in folder "user" are always used
# (c) Sequences of type "EC" are always of priority, meaning that if there are sequences via the EC number (also in unrev), sequences of type "reaName" are not used, even if they are reviewed.
#     Thus the decision tree is:
#       has rev EC seqs? --no--> has unrev EC seqs? --no--> has rev reaName seqs? --no--> use unrev reaName seqs
# (d) if sequences are there in user file, files in "unrev" or "rev" are ignored (this might be changed in future versions)
# (e) sequences in folder "rev" are only used when seqSrc is in c("1","2","3") and when there are no user sequences
# (f) sequences in folder "unrev" are only used when (seqSrc is in c("3","4") OR (seqSrc is "2" AND reviewed seqs don't exist)) and when there are no user sequences
# (g) If sequence file does not exist or is empty, don't use it. Obviously.

seqfiles[, has_user_seqs := any(src == "user", na.rm = TRUE), by = .(rea, reaName, ecs)]
seqfiles[, has_revEC_seqs := sum(((type == "EC") * (src == "rev") * file_size), na.rm = TRUE) > 0, by = .(rea, reaName, ecs)]
seqfiles[, has_unrevEC_seqs := sum(((type == "EC") * (src == "unrev") * file_size), na.rm = TRUE) > 0, by = .(rea, reaName, ecs)]
seqfiles[, has_revRN_seqs := sum(((type == "reaName") * (src == "rev") * file_size), na.rm = TRUE) > 0, by = .(rea, reaName, ecs)]
seqfiles[, has_unrevRN_seqs := sum(((type == "reaName") * (src == "unrev") * file_size), na.rm = TRUE) > 0, by = .(rea, reaName, ecs)]

# (a+b)
seqfiles[src == "rxn", use := TRUE]
seqfiles[src == "user", use := TRUE]

if(seqSrc == "1") {
  # only reviewed
  # (c+d+e)
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == TRUE & src %in% c("rev") & type == "EC", use := TRUE]
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == FALSE & has_revRN_seqs == TRUE & src %in% c("rev") & type == "reaName", use := TRUE]
} else if(seqSrc == "2") {
  # only reviewed
  # (c+d+e) EC
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == TRUE & src %in% c("rev") & type == "EC", use := TRUE]
  # (c+d+f) EC
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == FALSE & has_unrevEC_seqs == TRUE & src %in% c("unrev") & type == "EC", use := TRUE]
  # (c+d+e) reaName
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == FALSE & has_unrevEC_seqs == FALSE & has_revRN_seqs == TRUE & src %in% c("rev") & type == "reaName", use := TRUE]
  # (c+d+f) reaName
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == FALSE & has_unrevEC_seqs == FALSE & has_revRN_seqs == FALSE & has_unrevRN_seqs == TRUE & src %in% c("unrev") & type == "reaName", use := TRUE]
} else if(seqSrc == "3") {
  # reviewed+unreviewed
  # (c+d+e+f)
  seqfiles[has_user_seqs == FALSE & src %in% c("unrev","rev") & type == "EC", use := TRUE]
  seqfiles[has_user_seqs == FALSE & has_revEC_seqs == FALSE & has_unrevEC_seqs == FALSE & src %in% c("unrev","rev") & type == "reaName", use := TRUE]
} else if(seqSrc == "4") {
  # only unreviewed
  # (c+d+e+f)
  seqfiles[has_user_seqs == FALSE & src %in% c("unrev") & type == "EC", use := TRUE]
  seqfiles[has_user_seqs == FALSE & has_unrevEC_seqs == FALSE & src %in% c("unrev") & type == "reaName", use := TRUE]
}

# (g)
seqfiles[fex == FALSE | file_size == 0, use := FALSE]

# for debugging
# seqfiles[use == TRUE, table(type, src)]

#-------------------------------------------------------------------------------
# (6) Combine all fasta files into one (keeping track of file source in headers)
#-------------------------------------------------------------------------------
allseqfiles <- unique(seqfiles[use == TRUE, paste0(sdir,"/",taxonomy,"/",file)])
allseqs <- mclapply(allseqfiles, mc.cores = n_threads, FUN = function(sf) {
  tmpseqs <- readAAStringSet(sf)
  files <- basename(sf)
  last_dir <- basename(dirname(sf))
  sf_short <- file.path(last_dir, files)
  names(tmpseqs) <- paste0(sf_short,"|",names(tmpseqs)) # TODO only e.g. rev/EC.fasta
  return(tmpseqs)
})
allseqs <- do.call("c",allseqs)
allseqs <- allseqs[width(allseqs) >= 80]

if(verbose >= 1)
  cat("Number of reference sequences used for alignments:",length(allseqs),"\n")

#-------------------------------------------------------------------------------
# (6) Data export
#-------------------------------------------------------------------------------

if(length(allseqs) == 0)
  allseqs <- AAStringSet()
writeXStringSet(allseqs, filepath = "query.faa")

# sequence headers
if(length(allseqs) > 0) {
  seq_headers <- data.table(header = names(allseqs))
  seq_headers[, file := sub("\\|.+$","",header)]
  seq_headers[, header := sub("^.+\\.fasta\\|","",header)]
} else {
  seq_headers <- data.table(header = character(0L),
                            file = character(0L))
}

# save central temporary data structures
save(pwyrea, reaec, seqfiles, seq_headers, file = "prealignment_data.RData")

# # debugging
# writeXStringSet(allseqs, filepath = "~/tmp/refs.faa")
# save(pwyrea, reaec, seqfiles, seq_headers, file = "~/tmp/prealignment_data.RData")
