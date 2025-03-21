# Get gapseq reaction database hits for reaction queries by EC and/or name. The
# function also identifies alternative ECs (for later reference sequence
# collection).
#
# Arguments:
# `rea` - Reaction ID in query data source (in most cases: metacyc rection ID)
# `reaName` - Reaction name
# `ec` - Enzyme Commission (EC) number
# `database` - Target database. "seed" or "vmh" (Support for vmh might be stopped soon)
#
getDBhit <- function(rea, reaName, ec, database) {
  dbhit <- c()

  ec <- unique(str_split(ec, "/")[[1]])

  # Check if EC numbers are complete and valid
  EC_test <- grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$",ec)

  # EC with escaped points
  ec_esc <- gsub(".","\\.", ec, fixed = TRUE)

  # get KEGG ID(s) for metacyc reaction (comma-separated)
  kegg <- metaRea[id == paste0("|",rea,"|"), kegg]

  # 1) search in target reaction DB by query EC (works with multiple EC input)
  if(any(EC_test)) {
    if(database == "vmh") {
      dbhit <- reaDB1[grepl(paste0(ec_esc[EC_test],"($| |,|;)", collapse = "|"), ecnumber), abbreviation]
    } else if(database == "seed") {
      dbhit <- seedEC[grepl(paste0(ec_esc[EC_test],"$", collapse = "|"), `External ID`), `MS ID`]
      dbhit <- unlist(strsplit(dbhit, "\\|"))
    }
  }

  # 2) search in target reaction DB by query KEGG ID
  if(length(kegg) == 1 && kegg != "") {
    if(database == "vmh") {
      hittmp <- reaDB1[grepl(gsub(",","|",kegg), keggId), abbreviation]
      dbhit <- c(dbhit, hittmp)
    } else if(database == "seed") {
      hittmp <- reaDB5[grepl(gsub(",","|",kegg), other), seed]
      hittmp <- unlist(strsplit(hittmp, "\\|"))
      dbhit <- c(dbhit, hittmp)
    }
  }

  # 3) search in target reaction DB by query alternative EC
  altec <- character(0L)
  altec_src <- character(0L)
  if(sum(EC_test) > 0) {
    for(i in 1:length(ec)) {
      if(EC_test[i]) {
        altectmp <- altecdb[grepl(paste0(ec_esc[i],"($|,)"),altecdb)]
        altectmp <- unlist(strsplit(altectmp,","))
        altectmp <- altectmp[altectmp %notin% ec]
        altec_srctmp <- rep("altec.csv", length(altectmp))

        if(is.null(altectmp)) {
          brendaec <- brenda[grepl(paste0(ec_esc[i],"($|,| |\t|;)"),brenda)]
          brendaec <- paste(brendaec, collapse = ";")
          brendaec <- str_extract_all(brendaec, "[0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+(?=$|,|;|\\.)")[[1]]
          bralttmp <- unique(brendaec[brendaec %notin% ec])
          if(length(bralttmp) <= 3) {
            # Too many matches (>3) would suggest ambiguity, so they are discarded.
            altectmp <- bralttmp
            altec_srctmp <- rep("brenda", length(altectmp))
          } else {
            altectmp <- character(0L)
            altec_srctmp <- character(0L)
          }
        }

        altec <- c(altec, altectmp)
        altec_src <- c(altec_src, altec_srctmp)

        if(length(altectmp) > 0) {
          altec_esc <- gsub(".","\\.", altectmp, fixed = TRUE)
          altec_esc_comb <- paste0(altec_esc,"(,|$|;)", collapse = "|")
          if(database == "vmh") {
            hittmp <- reaDB1[grepl(altec_esc_comb, ecnumber), abbreviation]
            dbhit <- c(dbhit, hittmp)
          } else if(database == "seed") {
            hittmp <- seedEC[`External ID` %in% altectmp, `MS ID`]
            hittmp <- unlist(strsplit(hittmp, "\\|"))
            dbhit <- c(dbhit, hittmp)
          }
        }
      }
    }
  }

  # 4) search in target reaction DB by metacycID (only when target database is vmh)
  if(database == "vmh" && length(rea) == 1 && rea != "") {
    hittmp <- reaDB2[grepl(paste0("META:",rea,"($|;)"),database_links),bigg_id]
    dbhit <- c(dbhit, hittmp)
  }

  # 5) match reaction using mnxref namespace
  if(database == "vmh") {
    hittmp <- reaDB6[other == rea, bigg]
    dbhit <- c(dbhit, hittmp)
  } else if(database == "seed") {
    hittmp <- reaDB5[other == rea, seed]
    dbhit <- c(dbhit, hittmp)
  }

  # 6) match reaction using custom/alternative enzyme-name - seedID mapping only
  if(database == "seed" && length(reaName) == 1 && reaName != "") {
    hittmp <- seedEnzymesNames[enzyme.name == reaName, seed.rxn.id]
    dbhit <- c(dbhit, hittmp)
  }

  return(list(dbhit = paste(sort(unique(dbhit)), collapse = " "),
              altec = paste(altec, collapse = "/"),
              altec_src = paste(altec_src, collapse = ",")))
}


# make sure script path is defined
if(!exists("script.dir"))
  stop("Script path not defined for 'getDBhit()' function")

# load database files
metaRea <- fread(paste0(script.dir, "/../dat/meta_rea.tbl"))
reaDB1 <- fread(paste0(script.dir, "/../dat/vmh_reactions.tsv"))
reaDB2 <- fread(paste0(script.dir, "/../dat/bigg_reactions.tbl"))
reaDB3 <- fread(paste0(script.dir, "/../dat/seed_reactions_corrected.tsv"))
reaDB4 <- fread(paste0(script.dir, "/../dat/mnxref_seed.tsv"))
reaDB5 <- fread(paste0(script.dir, "/../dat/mnxref_seed-other.tsv"))
reaDB6 <- fread(paste0(script.dir, "/../dat/mnxref_bigg-other.tsv"))
seedEC <- fread(paste0(script.dir, "/../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"))
altecdb <- readLines(paste0(script.dir, "/../dat/altec.csv"))
brenda <- readLines(paste0(script.dir, "/../dat/brenda_ec_edited.csv"))
seedEnzymesNames <- fread(paste0(script.dir, "/../dat/seed_Enzyme_Name_Reactions_Aliases.tsv"))
