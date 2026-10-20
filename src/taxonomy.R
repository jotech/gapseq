args <- commandArgs(trailingOnly = TRUE)

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

if(length(args) != 2)
  stop("The script 'src/taxonomy.R' requires exactly two postional command line arguments")

# Check for the most common cases before loaded the large taxonomy database
if(args[1] == "get_entry" && args[2] == "auto") {
  cat("auto\t0")
  quit(save = "no", status = 0)
}
if(args[1] == "get_entry" && grepl("^(B|b)acteria$",args[2])) {
  cat("Bacteria\t2")
  quit(save = "no", status = 0)
}
if(args[1] == "get_entry" && grepl("^(A|a)rchaea$",args[2])) {
  cat("Archaea\t2157")
  quit(save = "no", status = 0)
}
if(args[1] == "get_entry" && grepl("^(P|p)rokaryota$",args[2])) {
  cat("Prokaryota\tNA")
  quit(save = "no", status = 0)
}

# not Bacteria/Archaea/Eukaryota|Prokaryota case
if(args[1] == "get_entry") {
  suppressPackageStartupMessages(library(data.table))
  library(stringr)

  taxdb <- fread(paste0(script.dir,"/../dat/taxonomy_full.tsv.gz"))
  is_taxid <- !is.na(suppressWarnings(as.integer(args[2])))

  if(is_taxid) {
    taxid   <- args[2]
    taxname <- taxdb[taxid == as.integer(args[2]), name]
    if(length(taxname) == 0) {
      message("Error: Taxonomy ID \"",taxid,"\" not found in taxonomy database.", sep = "")
      quit(save = "no", status = 1)
    }
  } else {
    taxid   <- taxdb[tolower(name) == tolower(args[2]), taxid]
    taxname <- taxdb[tolower(name) == tolower(args[2]), name]
    if(length(taxid) == 0) {
      tmp_strdist <- adist(tolower(args[2]), tolower(taxdb$name))
      closest_names <- taxdb$name[order(tmp_strdist)[1:3]]
      message("Error: Taxon name \"",args[2],"\" not found in taxonomy database. Closest taxon names:\n",
              paste0("   ", closest_names, collapse = "\n"),
              sep = "")
      quit(save = "no", status = 1)
    }
    if(length(taxid) > 1) {
      message("Error: Taxon name \"",args[2],"\" matches mor than 1 taxon ID in the database:\n",
              paste0("   ", taxid, collapse = "\n"),
              sep = "")
      quit(save = "no", status = 1)
    }
  }

  cat(taxname,taxid, sep = "\t")
}


