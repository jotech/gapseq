download_ecdb <- function() {
  tmpxml <- tempfile()
  message("Downloading latest EC database from enzyme-database.org")
  fdl <- download.file("https://www.enzyme-database.org/downloads/enzyme-data.xml.gz", tmpxml)
  if(fdl != 0)
    stop("Download of 'enzyme-data.xml.gz' failed.")

  data <- read_xml(tmpxml)

  message("Parsing database XML")
  xml <- xmlParse(data)
  qwe <- xmlToList(xml)

  ecDB_date <- qwe[[1]]$table_structure$options["Create_time"]
  message("EC database date: ",ecDB_date)

  qwe <- qwe$database

  # 2 - citations
  # 4 - classes, subclasses, subsubclasses
  # 6 - enzymes
  # 8 - actions
  # 12 - references

  message("Collecting data for all ECs")
  enz <- qwe[[6]]
  n <- length(enz) - 1

  enzl <- list();
  for(i in 1:n) {
    message("\r",i,"/",n, appendLF = FALSE)
    row <- do.call("rbind", enz[[i]])
    tmp <- data.table(entry = i,
                      text = unname(unlist(row[, 1])),
                      attr = unname(unlist(row[, 2])))
    tmp <- tmp[attr != "true"]

    enzl[[i]] <- dcast(tmp, entry ~ attr, value.var = "text")
  }
  message(" done")

  enzDT <- rbindlist(enzl, fill = TRUE)
  enzDT[diagram == "diagram", diagram := NA_character_]
  enzDT[glossary == "glossary", glossary := NA_character_]
  enzDT[accepted_name == "accepted_name", accepted_name := NA_character_]
  enzDT[other_names == "other_names", other_names := NA_character_]
  enzDT[sys_name == "sys_name", sys_name := NA_character_]
  enzDT[cas_num == "cas_num", cas_num := NA_character_]
  enzDT[reaction == "reaction", reaction := NA_character_]
  enzDT[comments == "comments", comments := NA_character_]

  # number #8 we also need!
  act <- qwe[[8]]
  n <- length(act)-1

  actl <- list();
  message("Collecting data for all created/deleted/tranferred actions")
  for(i in 1:n) {
    message("\r",i,"/",n, appendLF = FALSE)
    row <- do.call("rbind", act[[i]])
    tmp <- data.table(entry = i,
                      text = unname(unlist(row[, 1])),
                      attr = unname(unlist(row[, 2])))
    tmp <- tmp[attr != "true"]

    actl[[i]] <- dcast(tmp, entry ~ attr, value.var = "text")
  }
  message(" done")
  actDT <- rbindlist(actl, fill = TRUE)
  actDT[note == "note", note := NA_character_]
  #actDT <- actDT[status %in% c("3-0","4-0")]

  out <- merge(enzDT[,.(ec_num, accepted_name, other_names, sys_name, cas_num, class, subclass, subsubclass, last_change_EC = last_change,reaction, comments)],
               actDT[,.(ec_num, action, history, status, note, last_change_action = last_change)], by = "ec_num")

  out[, ec_num_new := NA_character_]

  out[action == "deleted", ec_num_new := str_match(note, "\\\"new\\\"\\>EC \\s*(.*?)\\s*\\</a\\>")[,2]]
  out[action == "deleted" & is.na(ec_num_new) & str_count(note, "EC [0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+") == 1, ec_num_new := str_match(note, "EC \\s*([0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+?)\\s*( |\\.|,|$)")[,2]]

  out[action == "transferred", ec_num_new := str_match(note, "\\\"new\\\"\\>EC \\s*(.*?)\\s*\\</a\\>")[,2]]
  out[action == "transferred" & is.na(ec_num_new), ec_num_new := str_match(note, "\\>EC \\s*(.*?)\\s*\\</a\\>")[,2]]
  out[action == "transferred" & is.na(ec_num_new) & str_count(note, "EC [0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+") == 1, ec_num_new := str_match(note, "EC \\s*([0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+?)\\s*( |\\.|,|$)")[,2]]
  out[action == "transferred" & is.na(ec_num_new) & str_count(note, " [0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+") == 1, ec_num_new := str_match(note, " \\s*([0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+?)\\s*( |\\.|,|$)")[,2]]

  out[ec_num == ec_num_new, ec_num_new := NA_character_]

  message("Summary: ", paste0(names(out[, table(action)]),": ",out[, table(action)], collapse = "; "))

  return(out)
}

correctEC_pathways <- function(ecdb) {
  message("Updating ECs in dat/meta_pwy.tbl and dat/custom_pwy.tbl ...")
  pwym <- fread("dat/meta_pwy.tbl"); pwym$source <- "metacyc"
  pwyc <- fread("dat/custom_pwy.tbl"); pwyc$source <- "custom"
  pwym[, id := gsub("^\\||\\|$","", id)]
  pwyc[, id := gsub("^\\||\\|$","", id)]
  pwy <- rbind(pwym, pwyc)

  # some hard-coded correction of syntax errors in reaction names that would break string parsing (remove when meta_pwy.tbl is corrected/updated)
  pwy[id %in% c("PWY18C3-10","PWY18C3-12"), reaName := sub("sucrose-3-isobutanoyl-4-isovaleryl-3;-isovaleroyltransferase",
                                                           "sucrose-3-isobutanoyl-4-isovaleryl-3'-isovaleroyltransferase",
                                                           reaName)]
  pwy[id == "PWY-6024", reaName := sub("vitexin 2;;-O-arabinoside 7-O-galactosyltransferase",
                                       "vitexin 2''-O-arabinoside 7-O-galactosyltransferase",
                                       reaName)]

  # some reaction names contain HTML character entities such as "&pi;". In those cases, escape the semicolon ("\\;")
  pwy[, reaName := gsub("&([a-zA-Z0-9#]+);", "&\\1\\\\;", reaName)]

  # In rare cases the reaction name separator ";" occurs within a reaction name, when inside brackets "(...)". These semicolons need to be escaped ("\\;") before string splitting.
  pwy[, reaName := gsub("\\(([^()]*?);([^()]*?)\\)", "\\(\\1\\\\;\\2\\)",reaName)]

  cnames <- colnames(pwy)

  pwy_l <- list()
  pwy_updates <- list(); k <- 1
  for(i in 1:nrow(pwy)) {

    # get pwy data table with n rows (n - number of reactions)
    dttmp <- pwy[i, .(id, name, altname, hierarchy, taxrange,
                      reaId = unlist(str_split(reaId, ",")),
                      reaEc = unlist(str_split(reaEc, ",")),
                      keyRea,
                      reaName = unlist(str_split(reaName, "(?<!\\\\);")),
                      superpathway, status, spont, source)]
    dttmp <- dttmp[, .(reaEc = unlist(str_split(reaEc,"/"))),
                   by = .(id, name, altname, hierarchy, taxrange,reaId,keyRea,reaName,superpathway, status, spont, source)]

    # here we could also replace reaction names with their recommended ones


    # update ECs
    dttmp$new_ec <- ecdb[dttmp$reaEc, ec_num_new]

    if(any(!is.na(dttmp$new_ec))) {
      pwy_updates[[k]] <- copy(dttmp)
      msg <-  paste0("  Changing EC ",dttmp[!is.na(new_ec), reaEc]," to EC ",
                     dttmp[!is.na(new_ec), new_ec]," (",pwy[i,id],", ",pwy[i,source],")")
      msg <- paste(msg, collapse = "\n")
      message(msg)
      dttmp[!is.na(new_ec), reaEc := new_ec]
      dttmp[, new_ec := NULL]
      k <- k + 1
    }

    # replace transferred/deleted ECs with the new ones

    # get it back in 1-row format
    dttmp <- dttmp[, .(reaEc = paste(reaEc,collapse = "/")),
                   by = .(id, name, altname, hierarchy, taxrange,reaId,keyRea,reaName,superpathway, status, spont, source)]

    dttmp <- dttmp[, .(reaId = paste(reaId, collapse = ","),
                       reaEc = paste(reaEc, collapse = ","),
                       reaName = paste(reaName, collapse = ";"),
                       reaNr = .N,
                       ecNr = sum(reaEc != "" & !is.na(reaEc))),
                   by = .(id, name, altname, hierarchy, taxrange, keyRea,
                          superpathway, status, spont, source)]
    dttmp <- dttmp[, ..cnames]
    pwy_l[[i]] <- copy(dttmp)
  }
  if(length(pwy_updates) == 0) {
    message("  No EC corrections necessary for metacyc/custom pathways.")
  }
  pwy_upd <- rbindlist(pwy_l)

  pwy_upd[, reaName := gsub("\\\\;",";",reaName)]

  fwrite(pwy_upd[source == "metacyc", -"source"], file = "dat/meta_pwy.tbl",
         sep = "\t", quote = FALSE)
  fwrite(pwy_upd[source == "custom", -"source"], file = "dat/custom_pwy.tbl",
         sep = "\t", quote = FALSE)

  return(TRUE)
}

correctEC_conflicts <- function(ecdb) {
  fileOI <- "dat/ec_conflicts.tsv"
  message("Updating ECs in ",fileOI," ...")
  eccdt <- fread(fileOI)

  # # adhoc test
  # eccdt[3, num_a := "2.1.1.86"]
  # eccdt[7, num_b := "1.9.3.1"]
  # eccdt[14, num_b := "1.3.5.4"]

  eccdt[, num_a_new := ecdb[num_a, ec_num_new]]
  eccdt[, num_b_new := ecdb[num_b, ec_num_new]]

  if(all(is.na(eccdt$num_a_new)) & all(is.na(eccdt$num_b_new))) {
    message("  No EC corrections necessary in ec conflicts table.")
  }

  if(any(!is.na(eccdt$num_a_new))) {
    msg <-  paste0("  Changing EC ",eccdt[!is.na(num_a_new), num_a]," to EC ",
                   eccdt[!is.na(num_a_new), num_a_new])
    msg <- paste(msg, collapse = "\n")
    message(msg)
    eccdt[!is.na(num_a_new), num_a := num_a_new]
  }
  if(any(!is.na(eccdt$num_b_new))) {
    msg <-  paste0("  Changing EC ",eccdt[!is.na(num_b_new), num_b]," to EC ",
                   eccdt[!is.na(num_b_new), num_b_new])
    msg <- paste(msg, collapse = "\n")
    message(msg)
    eccdt[!is.na(num_b_new), num_b := num_b_new]
  }

  eccdt[, num_b_new := NULL]
  eccdt[, num_a_new := NULL]

  fwrite(eccdt, file = "dat/ec_conflicts.tsv",
         sep = "\t", quote = FALSE)

  return(TRUE)
}

correctEC_seedrxn <- function(ecdb) {
  fileOI <- "dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv"
  message("Updating ECs in ",fileOI," ...")
  dtec <- fread(fileOI)
  dtec[, new_ec := ecdb[`External ID`, ec_num_new]]
  dtec[, new_ec_exists := FALSE]
  existing_ecs <- dtec$`External ID`
  dtec[!is.na(new_ec), new_ec_exists := new_ec %in% existing_ecs]

  update_ecs <- unique(dtec[!is.na(new_ec),new_ec])

  new_entries <- list(); k <- 1
  for(uec in update_ecs) {
    tmp <- dtec[new_ec == uec]
    rxn <- sort(unique(unlist(str_split(tmp$`MS ID`, "\\|"))))
    old_ecs <- unique(tmp$`External ID`)
    old_ecs <- old_ecs[old_ecs != uec]
    new_entries[[k]] <- data.table(`MS ID` = paste(rxn, collapse = "|"),
                                   `Old MS ID` = NA_character_,
                                   `External ID` = uec,
                                   Source = "Enzyme Class")
    message(" Changing EC(s) ",paste(old_ecs, collapse = ", ")," to EC ",uec)
    k <- k + 1
  }

  if(length(new_entries) > 0) {
    dtec <- rbind(dtec[is.na(new_ec)],
                  rbindlist(new_entries), fill = TRUE)
  } else {
    message("  No EC corrections necessary in ec to rxn mapping table.")
  }

  # Since we are here: Sort entries nicely
  dtec[, L1 := str_match(`External ID`,"^([0-9]?)")[,2]]
  dtec[, L2 := str_match(`External ID`,"^[0-9]\\.([0-9]+?)")[,2]]
  dtec[, L3 := str_match(`External ID`,"^[0-9]\\.[0-9]+\\.([0-9]+?)")[,2]]
  dtec[, L4 := str_match(`External ID`,"^[0-9]\\.[0-9]+\\.[0-9]+\\.([0-9]+)$")[,2]]

  dtec <- dtec[order(L1, L2, L3, L4, `External ID`)]
  dtec[, L1 := NULL];dtec[, L2 := NULL];dtec[, L3 := NULL];dtec[, L4 := NULL]

  fwrite(dtec, fileOI, sep = "\t", quote = FALSE)

  return(TRUE)
}

correctEC_exceptions <- function(ecdb) {
  fileOI <- "dat/exception.tbl"
  message("Updating ECs in ",fileOI," ...")
  dtex <- fread(fileOI, skip = 1)

  header <- readLines(fileOI, n = 1)

  dtex[, new_ec := ecdb[`enzyme/reaction`, ec_num_new]]

  if(any(!is.na(dtex$new_ec))) {
    msg <-  paste0("  Changing EC ",dtex[!is.na(new_ec), `enzyme/reaction`]," to EC ",
                   dtex[!is.na(new_ec), new_ec])
    msg <- paste(msg, collapse = "\n")
    message(msg)
    dtex[!is.na(new_ec), `enzyme/reaction` := new_ec]
  } else {
    message("  No EC corrections necessary in exception table.")
    return(TRUE)
  }
  dtex[, new_ec := NULL]

  writeLines(header, fileOI)
  fwrite(dtex, fileOI, sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE)

  return(TRUE)
}

correctEC_userseqs <- function(ecdb) {
  message("Checking if fasta files in dat/seq/<Bacteria|Archaea>/user/ need to be renamed according to the EC number.")
  userseqs <- c(dir("dat/seq/Bacteria/user", full.names = TRUE, pattern = "[0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+\\.fasta$"),
                dir("dat/seq/Archaea/user", full.names = TRUE, pattern ="[0-9]\\.[0-9]+\\.[0-9]+\\.[0-9]+\\.fasta$"))
  dtusr <- data.table(path = userseqs,
                      ec = sub("\\.fasta","",basename(userseqs)))
  dtusr[, new_ec := ecdb[ec, ec_num_new]]
  dtusr$dir <- NA_character_

  if(any(!is.na(dtusr$new_ec))) {
    dtusr[!is.na(new_ec), new_path := paste0(dirname(path),"/",new_ec,".fasta")]
    msg <-  paste0("  Renaming EC ",dtusr[!is.na(new_ec), path]," to ",
                   dtusr[!is.na(new_ec), new_path])
    msg <- paste(msg, collapse = "\n")
    message(msg)
    for(i in which(!is.na(dtusr$new_ec))) {
      if(file.exists(dtusr[i, new_path])) {
        warning(paste0("Target file ",dtusr[i, new_path]," already exists."))
      } else {
        res <- file.rename(dtusr[i, path],dtusr[i, new_path])
      }
    }
  } else {
    message("  No EC corrections necessary for user sequences.")
  }

  return(TRUE)
}
