suppressMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

# positional arguments:
# 1. allPwy
# 2. key (search term/pattern)
# 3. Search type / column
# 4. true/false (include Super-Pathways?)

allpwy <- fread(args[1], sep = "\t", quote = FALSE)
keys <- unlist(strsplit(args[2],"\\|"))
stype <- args[3]
excludeSuper <- as.logical(args[4])

if(stype == "id") {
  allpwy <- allpwy[id %in% keys]
} else if(stype == "hierarchy") {
  keys <- tolower(keys)
  keep <- allpwy[, any(tolower(unlist(strsplit(hierarchy,"\\|"))) %in% keys), by = id][V1 == TRUE, id]
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
  cat("Pathways ignored because no reactions found:",length(indrm))
  cat(paste0("Pathways without reactions:\n",paste(allpwy$id[indrm], collapse = ", ")))
  allpwy <- allpwy[-indrm]
}

# save.image(file = "~/test.Rdata") # debug line

fwrite(allpwy, file = args[1], sep = "\t", quote = FALSE, col.names = FALSE)
