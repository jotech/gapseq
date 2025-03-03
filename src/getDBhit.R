# Get database hits for reaction queries (by EC and/or name)
#
# Arguments:
# `rea` - Reaction ID in query data source (in most cases: metacyc rection ID)
# `reaName` - Reaction name
# `ec` - Enzyme Commission (EC) number
# `database` - Target database. "seed" or "vmh" (Support for vmh might be stopped soon)
#
getDBhit(rea, reaName, ec, database) {
  # Check if EC number is complete and valid
  EC_test <- grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$",ec)

  # EC with escaped points
  ec_esc <- gsub(".","\\.", ec, fixed = TRUE)

  # get KEGG ID(s) for metacyc reaction (comma-separated)
  kegg <- metaRea[id == paste0("|",rea,"|"), kegg]

  # 1) search in target reaction DB by query EC
  if(EC_test) {
    if(database == "vmh") {
      dbhit <- reaDB1[grepl(paste0(ec_esc,"($| |,|;)"), ecnumber), abbreviation]
    } else if(database == "seed") {
      dbhit <- seedEC[grepl(paste0(ec_esc,"$"), `External ID`), `MS ID`]
      dbhit <- unlist(strsplit(dbhit, "\\|"))
    }
  }

  # 2) search in target reaction DB by query KEGG ID
  if(length(kegg) == 1 && kegg != "") {
    if(database == "vmh") {
      hittmp <- reaDB1[grepl(gsub(",","|",kegg), keggId), abbreviation]
      dbhit <- c(dbhit, hittmp)
    } else if(database == "seed") {
      hittmp <- seedEC[grepl(paste0(ec_esc,"$"), `External ID`), `MS ID`]
      hittmp <- unlist(strsplit(hittmp, "\\|"))
      dbhit <- c(dbhit, hittmp)
    }
  }

  # 3) search in target reaction DB by query alternative EC
  if(EC_test) {
    altec <- altecdb[grepl(paste0(ec_esc,"($|,)"),altecdb)]
    altec <- unlist(strsplit(altec,","))
    altec <- altec[altec != ec]

    if(is.null(altec)) {
      brendaec <- brenda[grepl(paste0(ec_esc,"($|,| |\t|;)"),brenda)]
      brendaec <- paste(brendaec, collapse = ";")
      brendaec <- str_extract_all(brendaec, "[0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+(?=$|,|;|\\.)")[[1]]
      altec <- unique(brendaec[brendaec != ec])
    }

    if(length(altec) > 0) {
      altec_esc <- gsub(".","\\.", altec, fixed = TRUE)
      altec_esc_comb <- paste0(altec_esc,"(,|$|;)", collapse = "|")
      if(database == "vmh") {
        hittmp <- reaDB1[grepl(altec_esc_comb, ecnumber), abbreviation]
        dbhit <- c(dbhit, hittmp)
      } else if(database == "seed") {
        hittmp <- seedEC[`External ID` %in% altec, `MS ID`]
        dbhit <- c(dbhit, hittmp)
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

  return(list(dbhit = order(unique(dbhit)),
              altec = altec))
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
