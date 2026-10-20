suppressPackageStartupMessages(library(data.table))
library(stringr)

args = commandArgs(trailingOnly=TRUE)

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

# positional arguments:
# 1. allPwy
# 2. key (search term/pattern)
# 3. Search type / column
# 4. true/false (exclude Super-Pathways?)
# 5. taxonomy range as NCBI taxon ID

allpwy <- fread(args[1], sep = "\t", quote = FALSE)
keys <- unlist(strsplit(args[2],"\\|"))
stype <- args[3]
excludeSuper <- as.logical(args[4])
taxid <- ifelse(args[5] == "NA", NA_integer_, as.integer(args[5]))

if(stype == "id") {
  allpwy <- allpwy[id %in% keys]
} else if(stype == "hierarchy") {
  keys <- tolower(keys)
  keep <- allpwy[, any(tolower(unlist(strsplit(hierarchy,"\\|"))) %in% keys) | any(tolower(unlist(strsplit(hierarchy,"\\;"))) %in% keys), by = id][V1 == TRUE, id]
  allpwy <- allpwy[id %in% keep]
} else {
  keypattern <- paste0("\\b(",keys,")\\b")
  allpwy <- allpwy[grepl(keypattern, id, ignore.case = TRUE, perl = TRUE) |
                     grepl(keypattern, hierarchy, ignore.case = TRUE, perl = TRUE) |
                     grepl(keypattern, name, ignore.case = TRUE, perl = TRUE) |
                     grepl(keypattern, altname, ignore.case = TRUE, perl = TRUE)]
}

if(excludeSuper) {
  allpwy <- allpwy[!(as.logical(superpathway) == TRUE | grepl("Super-Pathways", hierarchy))]
}

if(any(allpwy$reaId == "" | is.na(allpwy$reaId))) {
  indrm <- which(allpwy$reaId == "" | is.na(allpwy$reaId))
  cat("Pathways ignored because no reactions found:",length(indrm),"\n")
  cat(paste0("Pathways without reactions:\n",paste(allpwy$id[indrm], collapse = ", ")))
  allpwy <- allpwy[-indrm]
}

# save.image(file = "~/test.Rdata") # debug line

# taxonomy range filter
if(is.na(taxid) || taxid != 0) {
  # (1) get all tax ids that are assigned to at least one pathway:
  pres_taxids <- as.integer(unique(gsub("\\|","",sub("TAX-","",unlist(str_split(allpwy$taxrange,","))))))
  pres_taxids <- pres_taxids[!is.na(pres_taxids)]

  # (2) get all ancestral taxon IDs
  taxdb <- fread(paste0(script.dir,"/../dat/taxonomy_full.tsv.gz"))
  domids <- taxdb[rank == "domain", taxid]
  anc_ids <- integer(0L)
  if(!is.na(taxid)) {
    tmpid <- taxid
    while(TRUE) {
      # print(tmpid)
      if(tmpid %in% domids || is.na(tmpid))
        break
      tmpid <- taxdb[taxid == tmpid, parent_taxid]
      anc_ids <- c(anc_ids, tmpid)
    }
  }

  # (3) get all child tax IDs
  child_ids <- integer(0L)
  if(is.na(taxid))
    taxid <- c(2, 2157) # Prokaryota case
  tmpids <- taxid
  while(TRUE) {
    tmpids <- tmpids[tmpids %notin% child_ids]
    child_ids <- c(child_ids, tmpids)
    if(length(tmpids) == 0)
      break
    tmpids <- taxdb[parent_taxid %in% tmpids, taxid]
  }

  # combine all IDs and retain only the ones that are also in "pres_taxids"
  comb_reltax <- c(taxid, anc_ids, child_ids, 131567)  # 131567 is "cellular organisms"
  rel_ids <- intersect(comb_reltax, pres_taxids)

  # perfom the actual taxonomy filter
  taxregex <- paste0("\\|TAX-",rel_ids,"\\|", collapse = "|")
  allpwy <- allpwy[grepl(taxregex, taxrange) | is.na(taxrange) | taxrange == ""]
}

fwrite(allpwy, file = args[1], sep = "\t", quote = FALSE, col.names = FALSE)


#
# tmp_meta <- fread("dat/meta_pwy.tbl")
# tmp_meta[taxrange == "" | is.na(taxrange)]
#
# tmp_cust <- fread("dat/custom_pwy.tbl")
# tmp_cust[taxrange == "" | is.na(taxrange), .(id, name, reaEc)]
