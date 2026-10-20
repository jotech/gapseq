library(data.table)
library(stringr)

# full NCBI taxonomy DB
taxdump <- tempfile()
taxdump_dir <- tempdir()
download.file("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", taxdump, method = "curl")

untar(taxdump, exdir = taxdump_dir)

names <- fread(paste0(taxdump_dir,"/names.dmp"), header = FALSE, sep = "|", quote = "\t")
names <- names[V4 == "scientific name"]
names[, V4 := NULL]
names[, V5 := NULL]
setnames(names, c("taxid", "name", "altname"))

nodes <- fread(paste0(taxdump_dir,"/nodes.dmp"), header = FALSE, sep = "|", quote = "\t")
nodes <- nodes[, .(taxid = V1, parent_taxid = V2, rank = V3)]

taxdt <- merge(names, nodes, by = "taxid")
taxdt[, altname := NULL]
taxdt <- taxdt[!grepl("(clinical|environmental) samples", name)]

# retain only taxa from the domains Bacteria, Archaea, Eukaryota
# Initialize domain assignments
domains <- as.character(taxdt[rank == "domain", taxid])
taxdt[, domain := NA_character_]
taxdt[, taxid := as.character(taxid)]
taxdt[, parent_taxid := as.character(parent_taxid)]
setkey(taxdt, taxid)
taxdt[, domain := fifelse(taxid %in% domains, taxid, NA_character_)]

# Propagate domain assignments from parents to children
repeat {
  missing <- which(is.na(taxdt$domain))

  parent_domain <- taxdt$domain[
    match(taxdt$parent_taxid[missing], taxdt$taxid)
  ]

  update <- which(!is.na(parent_domain))

  if (length(update) == 0L) break

  set(
    taxdt,
    i = missing[update],
    j = "domain",
    value = parent_domain[update]
  )
}

taxdt <- taxdt[!is.na(domain)]
setorder(taxdt, name, rank, taxid)
fwrite(taxdt, "dat/taxonomy_full.tsv.gz", sep = "\t", quote = FALSE, compressLevel = 9)
