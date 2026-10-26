library(data.table)
library(parallel)
archive <- "/mnt/nuuk/Resources/PMN/plantcyc.7z"
pwyfile <- "plantcyc/17.0.0/data/pathways.dat"
rxnfile <- "plantcyc/17.0.0/data/reactions.dat"
clsfile <- "plantcyc/17.0.0/data/classes.dat"
enzfile <- "plantcyc/17.0.0/data/enzrxns.dat"

# A central helper function:
get_info <- function(datlines, field, collapse_sep = NA, presuf = NA) {
  if(field == "SPONTANEOUS?")
    field <- "SPONTANEOUS\\?"
  idx <- grep(paste0("^",field, " - "), datlines)
  if(length(idx) == 0)
    return(NA_character_)

  content <- sub(paste0("^",field, " - "),"", datlines[idx])

  if(!is.na(presuf))
    content <- paste0(presuf, content, presuf)

  if(is.na(collapse_sep))
    return(content)

  return(paste(content, collapse = collapse_sep))
}

# Pathways
cmd <- sprintf(
  "7z x -so %s %s",
  shQuote(archive),
  shQuote(pwyfile)
)

con <- pipe(cmd, open = "r")
pwylines <- readLines(con, warn = FALSE)
close(con)
pwylist <- split(pwylines, cumsum(pwylines == "//"))
pwylist <- pwylist[-length(pwylist)]

# Reactions
cmd <- sprintf(
  "7z x -so %s %s",
  shQuote(archive),
  shQuote(rxnfile)
)

con <- pipe(cmd, open = "r")
rxnlines <- readLines(con, warn = FALSE)
close(con)
rxnlist <- split(rxnlines, cumsum(rxnlines == "//"))
rxnlist <- rxnlist[-length(rxnlist)]

# Classes
cmd <- sprintf(
  "7z x -so %s %s",
  shQuote(archive),
  shQuote(clsfile)
)

con <- pipe(cmd, open = "r")
clslines <- readLines(con, warn = FALSE)
close(con)
clslist <- split(clslines, cumsum(clslines == "//"))
clslist <- clslist[-length(clslist)]

# Enzymatic reactions
cmd <- sprintf(
  "7z x -so %s %s",
  shQuote(archive),
  shQuote(enzfile)
)

con <- pipe(cmd, open = "r")
enzlines <- readLines(con, warn = FALSE)
close(con)
enzlist <- split(enzlines, cumsum(enzlines == "//"))
enzlist <- enzlist[-length(enzlist)]

#-------------------------------------------------------------------------------
# Extract data to data.tables
#-------------------------------------------------------------------------------

# Classes
clsinfo <- mclapply(clslist, FUN = function(datlines) {
  out <- c()
  out["id"] <- get_info(datlines,
                        "UNIQUE-ID")
  out["name"] <- get_info(datlines,
                          "COMMON-NAME")
  out["parent_id"] <- get_info(datlines,
                               "TYPES")
  out <- data.table(id = out["id"],
                    name = out["name"],
                    parent_id = out["parent_id"])
  return(out)
})
clsinfo <- rbindlist(clsinfo)
clsinfo[is.na(name), name := id]
setkey(clsinfo, "id")

# Enzymatic reactions
enzinfo <- mclapply(enzlist, FUN = function(datlines) {
  out <- c()
  out["id"] <- get_info(datlines,
                        "UNIQUE-ID")
  out["name"] <- get_info(datlines,
                          "COMMON-NAME")
  out["reaction"] <- get_info(datlines,
                              "REACTION")
  out <- data.table(id = out["id"],
                    name = out["name"],
                    reaction = out["reaction"])
  return(out)
})
enzinfo <- rbindlist(enzinfo)
setkey(enzinfo, "reaction")

# Reactions
rxninfo <- mclapply(rxnlist, FUN = function(datlines) {
  out <- c()
  out["id"] <- get_info(datlines,
                        "UNIQUE-ID")
  out["name"] <- get_info(datlines,
                          "COMMON-NAME")
  # out["name"] <- clsinfo[out["id"], name]
  out["ec"] <- paste(unique(get_info(datlines,
                                     "EC-NUMBER")), collapse = "/")
  out["spont"] <- get_info(datlines,
                           "SPONTANEOUS?","/")
  out <- data.table(id = out["id"],
                    name = out["name"],
                    ec = out["ec"],
                    spont = out["spont"])

  out$enzname <- NA_character_

  tmpenzn <- enzinfo[out$id, name]

  if(length(tmpenzn) > 0)
    out$enzname <- tmpenzn[1]

  return(out)
})
rxninfo <- rbindlist(rxninfo)
rxninfo[!is.na(enzname), name := enzname]
rxninfo[is.na(name), name := id]
rxninfo[, enzname := NULL]
setkey(rxninfo, "id")

# Pathways
pwy_eoi <- data.table(varid = c("id",
                                "name",
                                "altname",
                                "keyRea",
                                "types",
                                "taxrange"),
                      field = c("UNIQUE-ID",
                                "COMMON-NAME",
                                "SYNONYMS",
                                "KEY-REACTIONS",
                                "TYPES",
                                "TAXONOMIC-RANGE"),
                      sep = c(NA_character_,
                              NA_character_,
                              ",",
                              ",",
                              ",",
                              ","),
                      presuf = c(NA_character_,
                                 NA_character_,
                                 NA_character_,
                                 NA_character_,
                                 NA_character_,
                                 "|"))

pwyinfo <- mclapply(pwylist, FUN = function(datlines) {
  out <- c()
  for(j in 1:nrow(pwy_eoi)) {
    out[pwy_eoi[j,varid]] <- get_info(datlines,
                                      pwy_eoi[j,field],
                                      pwy_eoi[j,sep],
                                      pwy_eoi[j,presuf])
  }
  #out <- data.table(t(out))

  # get reaction IDs
  reaids <- get_info(datlines, "REACTION-LIST")

  # get hierarchy
  tmptypes <- get_info(datlines, "TYPES")
  tmptypes <- tmptypes[tmptypes != "Super-Pathways"]

  itypes <- tmptypes
  otypes <- tmptypes

  while(length(itypes)>0) {
    itypes <- clsinfo[id %in% itypes, parent_id]
    itypes <- itypes[itypes %notin% otypes]
    otypes <- c(itypes, otypes)
  }
  out <- data.table(id = out["id"],
                    name = out["name"],
                    altname = out["altname"],
                    keyRea = out["keyRea"],
                    taxrange = out["taxrange"],
                    hierarchy = paste0("|",otypes,"|", collapse = ","),
                    superpathway = grepl("Super-Pathways", out["types"]),
                    reaId = reaids,
                    rea_is_pwy = reaids %notin% rxninfo$id)

  return(out)
})
names(pwyinfo) <- unlist(lapply(pwyinfo, FUN = function(x) x$id[1]))

# resolvinbg nested pathways
k <- 1
while(any(rbindlist(pwyinfo)$reaId %notin% rxninfo$id) & k <= 5) {
  print(k)
  pwyinfo <- mclapply(pwyinfo, FUN = function(pwy) {
    if(any(pwy$rea_is_pwy)) {
      for(i in 1:nrow(pwy)) {
        if(pwy[i, rea_is_pwy] == TRUE) {
          newrea <- paste(pwyinfo[[pwy[i, reaId]]]$reaId, collapse = ",")
          pwy[i, reaId := newrea]
        }
      }

      pwy[, rea_is_pwy := NULL]

      pwy <- pwy[
        , .(reaId = unlist(strsplit(reaId, ",", fixed = TRUE))),
        by = setdiff(names(pwy), "reaId")
      ]

      pwy[, rea_is_pwy := reaId %notin% rxninfo$id]

      return(pwy)
    } else {
      #pwy[, rea_is_pwy := reaId %notin% rxninfo$id]
      return(pwy)
    }
  })
  k <- k + 1
}

pwyinfo <- rbindlist(pwyinfo)


pwyinfo[, reaEc := rxninfo[reaId, ec]]
pwyinfo[reaEc == "NA", reaEc := ""]
pwyinfo[, reaEc := gsub("EC-","",reaEc)]
pwyinfo[, reaName := rxninfo[reaId, name]]

pwyinfo[, spontan := ifelse(rxninfo[reaId, spont] == "T", TRUE, FALSE)]
pwyinfo[is.na(spontan), spontan := FALSE]

# substitute html chars in raeaction names
pwyinfo[, reaName := gsub("&([a-zA-Z0-9#]+);", "\\1", reaName)]
pwyinfo[, reaName := gsub("<.*?>","", reaName)]

pwyinfo <- pwyinfo[, .(reaId = paste(reaId, collapse = ","),
                       reaEc = paste(reaEc, collapse = ","),
                       reaName = paste(reaName, collapse = ";"),
                       reaNr = .N,
                       ecNr = sum(reaEc != ""),
                       status = TRUE,
                       spont = paste(reaId[spontan], collapse = ",")),
                   by = .(id, name, altname, keyRea, taxrange, hierarchy, superpathway)]
pwyinfo <- pwyinfo[, .(id, name, altname, hierarchy, taxrange, reaId, reaEc, keyRea, reaName, reaNr, ecNr, superpathway, status, spont)]

fwrite(pwyinfo, "dat/plantcyc_pwy.tbl", sep = "\t", quote = FALSE)
