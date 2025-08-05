suppressMessages(library(data.table))
library(stringr)

#-------------------------------------------------------------------------------
# Read alignment stats and sequence headers
#-------------------------------------------------------------------------------

alicols <- c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","sseq")
if(file.size("out.tsv") == 0) {
  message("No hits found after alignments.")
  outcols <- c("id","tc","sub","sub_gapseq","exid","rea",alicols[1:(length(alicols)-1)],"comment")
  out_dummy <- matrix(nrow = 0, ncol = length(outcols))
  colnames(out_dummy) <- outcols
  fwrite(data.table(out_dummy), "transporter.tbl", sep = "\t", quote = FALSE)
  quit(save = "no")
} else {
  alignments <- fread("out.tsv") # alignments <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/out.tsv")
}
setnames(alignments, alicols[1:ncol(alignments)])

fasta_header.small <- readLines("fasta_header.small") # fasta_header.small <- readLines("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/fasta_header.small")
tcdb_all <- fread("tcdb_all", header = FALSE) # tcdb_all <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/tcdb_all", header = FALSE)
setnames(tcdb_all, c("tcid","substrates"))
tcdb_all[, substrates := gsub("CHEBI\\:[0-9]+","", substrates)]
tcdb_all[, substrates := gsub("\\|","", substrates)]
tcdb_all[, substrates := gsub("\\([0-9](\\+|-)\\)","", substrates)] # remove charges
tcdb_all[, substrates := gsub("\\((R|S|-)\\)-","", substrates)] # remove charges
rmterms <- "\\bion\\b|\\bwater\\b|\\bmolecule\\b|\\bmetabolite\\b|inorganic cation|metal cation|organic cation|\\bcation\\b|\\banion\\b|hydron|electron|proton"
tcdb_all <- tcdb_all[,substrates := gsub(rmterms,"", substrates)]
tcdb_all[, substrates := gsub("(ribonucleotide)n+m","", substrates, fixed = TRUE)]
tcdb_all[, substrates := gsub("(+)-","", substrates, fixed = TRUE)]
tcdb_all[, substrates := gsub("(+)","", substrates, fixed = TRUE)]
tcdb_all[, substrates := gsub(";+",";", substrates)]
tcdb_all[, substrates := gsub("^;|;$","", substrates)]
tcdb_all[, substrates := str_squish(substrates)]
tcdb_all <- tcdb_all[!grepl("^( )*$", substrates)]
tcdb_all[, substrates := str_replace_all(substrates, "([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1")]

SUBkey <- str_trim(readLines("SUBkey")) # SUBkey <- str_trim(readLines("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/SUBkey"))
SUBid <- fread("SUBid", header = FALSE) # SUBid <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/SUBid", header = FALSE)
SUBid[, V1 := tolower(V1)]
allDB <- fread("allDB") # allDB <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/allDB")
allDB[, metid := sub("\\[e.\\]$","",exmet)]
allDB <- allDB[, .(rea = paste(unique(id), collapse = ",")), by = .(type, metid)]

#-------------------------------------------------------------------------------
# (1) Parse arguments and gapseq/script path; and load gapseq data files
#-------------------------------------------------------------------------------
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

metDB <- fread(paste0(script.dir,"/../dat/seed_metabolites_edited.tsv"))
metDB <- metDB[,.(metid = id, name)]
setkey(metDB, "metid")

bitcutoff <- as.numeric(args[1]) # bitcutoff <- 200
identcutoff <- as.numeric(args[2]) # identcutoff <- 0
nouse_alternatives <- args[3] == "true" # nouse_alternatives <- FALSE
verbose <- as.integer(args[4])

# Define TC types (types 5, 6, 8 and 9 are ignored)
TC_types <- c("1.Channels and pores",
              "2.Electrochemical potential-driven transporters",
              "3.Primary active transporters",
              "4.Group translocators")

blasthits <- alignments[pident >= identcutoff]
blasthits[, tc := str_match(blasthits$qseqid, "([1-4]\\.[A-z]\\.[0-9]+\\.[0-9]+\\.[0-9]+)")[,1]]
blasthits <- blasthits[!is.na(tc)]

# assign substrate names to TCIDs in Blast table
findsubs <- function(tc) {
  # (a) Find substrates based on TCID
  type.i <- as.numeric(substr(tc,1,1))
  type.str <- TC_types[type.i]
  substrates <- unique(unlist(str_split(tcdb_all[tcid == tc, substrates], ";")))
  substrates_clean <- tolower(substrates)
  sub_hit <- character(0L)
  if(length(substrates_clean) > 0) {
    sub_hit <- SUBkey[str_detect(tolower(SUBkey), regex(paste0("^\\b", substrates_clean, "\\b$", collapse = "|"), ignore_case = TRUE))]
  }

  if(length(sub_hit) == 0) {
    # (b) Find substrates in header sequence of TC
    fa_headers <- fasta_header.small[grepl(paste0(tc," "), fasta_header.small, fixed = TRUE)]
    # fa_wors <- str_extract_all(fa_headers, "\\b\\w+\\b")[[1]]
    sub_hit <- sapply(SUBkey, function(key) {
      # Create word-boundary-aware regex pattern (escape special characters)
      pattern <- paste0("\\b", gsub("([\\W])", "\\\\\\1", key), "\\b")
      any(grepl(pattern, fa_headers, ignore.case = TRUE))
    })
    sub_hit <- SUBkey[sub_hit]
  }
  return(paste(sub_hit, collapse = ";"))
}

blasthits[, sub := findsubs(tc), by = tc]

splitcol2rows_mget <- function(dtInput, col2split, sep){
  dtInput <- dtInput[, .(tmp.add.col = unlist(strsplit(get(col2split),sep,T))), by=names(dtInput)]

  dtInput[, c(col2split):=NULL];
  setnames(dtInput, 'tmp.add.col', col2split);
  return(dtInput);
}

blasthits <- splitcol2rows_mget(blasthits, "sub", ";")
blasthits[, subtmp := tolower(sub)]

blasthits <- merge(blasthits, SUBid, by.x = "subtmp", by.y = "V1", allow.cartesian = TRUE)
setnames(blasthits, "V2", "exid")

# remove query-target-exid duplicate entries (take only the one with the highest bitscore)
blasthits <- blasthits[order(qseqid, stitle,-bitscore)]
blasthits <- unique(blasthits, by = c("qseqid","stitle","exid"))

# find reaction reaction IDs (gapseq database)
blasthits <- blasthits[, metid := gsub("^EX_|_e.$","",exid)]
blasthits[, type.i := as.numeric(substr(tc,1,1))]
blasthits[, type := TC_types[type.i]]
blasthits <- merge(blasthits, allDB[, .(type, metid, rea)], by = c("type","metid"), all.x = TRUE)
blasthits[!is.na(rea), comment := "transporter"]

# add gapseq database metabolite names
blasthits <- merge(blasthits, metDB, by = "metid")

fwr <- blasthits[!is.na(rea) & bitscore >= bitcutoff][!duplicated(metid), sort(metid)]
fwr_subopt <- blasthits[!is.na(rea)][!duplicated(metid), sort(metid)]
fwor <- blasthits[is.na(rea) & bitscore >= bitcutoff & metid %notin% fwr][!duplicated(metid), sort(metid)]
fwor_alt <- blasthits[is.na(rea) & metid %notin% fwr_subopt][!duplicated(metid), sort(metid)]

if(!nouse_alternatives)
for(altsub in fwor_alt) {
  candidates <- sort(unlist(str_split(allDB[metid == altsub, rea], ",")))
  if(length(candidates) > 0) {
    blasthits[metid == altsub, rea := paste(candidates, collapse = ",")]
  }
  blasthits[metid == altsub, comment := "alt-transporter"]
}

blasthits$comment <- factor(blasthits$comment, levels = c("transporter","alt-transporter"))
blasthits <- blasthits[order(comment,tc, -bitscore)]

#-------------------------------------------------------------------------------
# (n-1) Report results in logs
#-------------------------------------------------------------------------------
if(verbose > 0) {
  # Found transporters
  if(length(fwr) > 0) {
    cat("\nFound transporter for substances:\n",
        paste(sort(metDB[fwr, name]), collapse = "\n"),"\n", sep = "")
  }

  # Found transporters without reactions in gapseq's DB
  if(length(fwor) > 0) {
    cat("\nFound transporter without reactions in gapseq's DB:\n",
        paste(sort(metDB[fwor, name]), collapse = "\n"),"\n", sep = "")
  }

  # Alternative transport reactions types assigned to found transporters
  altassigned <- blasthits[comment == "alt-transporter" & !is.na(rea), unique(name)]
  if(length(altassigned) > 0) {
    cat("\nAdded alternative transport reactions for:\n",
        paste(sort(altassigned), collapse = "\n"),"\n", sep = "")
  }
  cat("\n")
}

#-------------------------------------------------------------------------------
# (n) Export reaction and pathway tables
#-------------------------------------------------------------------------------

# clean up
blastcols <- c(alicols[1:ncol(alignments)],"comment")
outdt <- cbind(blasthits[,.(id = qseqid, tc, sub, sub_gapseq = name, exid, rea)],
               blasthits[, ..blastcols])

# export
fwrite(outdt, "transporter.tbl", sep = "\t", quote = FALSE)
