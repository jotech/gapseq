suppressMessages(library(data.table))
library(stringr)

#-------------------------------------------------------------------------------
# Read alignment stats and sequence headers
#-------------------------------------------------------------------------------

alicols <- c("qseqid","pident","evalue","bitscore","qcovs","stitle","sstart","send","sseq")
if(file.size("out.tsv") == 0) {
  alignments <- data.table(matrix(0,nrow = 1, ncol = length(alicols)-1))
  alignments[, V1 := as.character(V1)]
  alignments[, V6 := as.character(V6)]
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
tcdb_all[, substrates := gsub(";+",";", substrates)]
tcdb_all[, substrates := gsub("^;|;$","", substrates)]
tcdb_all[, substrates := str_squish(substrates)]
tcdb_all <- tcdb_all[!grepl("^( )*$", substrates)]

SUBkey <- str_trim(readLines("SUBkey")) # SUBkey <- str_trim(readLines("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/SUBkey"))
SUBid <- fread("SUBid", header = FALSE) # SUBid <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/SUBid", header = FALSE)
SUBid[, V1 := tolower(V1)]
allDB <- fread("allDB") # allDB <- fread("~/test/transl/tmpdirNew/HRGMv2_0240.faa_Rufax1/allDB")

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

bitcutoff <- as.numeric(args[1]) # bitcutoff <- 200
identcutoff <- as.numeric(args[2]) # identcutoff <- 0

#
TC_types <- c("1.Channels and pores",
              "2.Electrochemical potential-driven transporters",
              "3.Primary active transporters",
              "4.Group translocators")



blasthits <- alignments[pident >= identcutoff]
blasthits[, tc := str_match(blasthits$qseqid, "([1-4]\\.[A-z]\\.[0-9]+\\.[0-9]+\\.[0-9]+)")[,1]]

# assign substrate names to tcs in Blast table
findsubs <- function(tc) {
  # (a) Find substrates based on TCID
  type.i <- as.numeric(substr(tc,1,1))
  type.str <- TC_types[type.i]
  substrates <- unique(unlist(str_split(tcdb_all[tcid == tc, substrates], ";")))
  substrates_clean <- tolower(substrates)
  sub_hit <- character(0L)
  if(length(substrates_clean) > 0) {
    sub_hit <- SUBkey[str_detect(tolower(SUBkey), regex(paste0("\\b", substrates_clean, "\\b", collapse = "|"), ignore_case = TRUE))]
  }

  if(length(sub_hit) == 0) {
    # (b) Find substrates in header sequence of TC
    fa_headers <- fasta_header.small[grepl(tc, fasta_header.small, fixed = TRUE)]
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

blasthits <- merge(blasthits, SUBid, by.x = "subtmp", by.y = "V1")
setnames(blasthits, "V2", "exid")

# remove query-target duplicate entries (take only the one with the highest bitscore)
blasthits <- blasthits[order(qseqid, stitle,-bitscore)]
blasthits <- unique(blasthits, by = c("qseqid","stitle"))

# find reaction reaction IDs (gapseq database)
blasthits <- blasthits[, metid := gsub("^EX_|_e.$","",exid)]







TC_blasthits <- unique(sort(str_match(blasthits$qseqid, "([1-4]\\.[A-z]\\.[0-9]+\\.[0-9]+\\.[0-9]+)")[,1]))

DBmissing <- data.table(sub = character(0L),
                        exmet = character(0L),
                        )
DBmissing_bitscore <- list()

tc2db <- function(tc) {
  tc <- "2.A.1.10.1"

  # (a) Find substrates based on TCID
  type.i <- as.numeric(substr(tc,1,1))
  type.str <- TC_types[type.i]
  substrates <- unique(unlist(str_split(tcdb_all[tcid == tc, substrates], ";")))
  substrates_clean <- tolower(substrates)
  sub_hit <- character(0L)
  if(length(substrates_clean) > 0) {
    sub_hit <- SUBkey[str_detect(tolower(SUBkey), regex(paste0("\\b", str_split(substrates_clean, "\\s+")[[1]], "\\b"), ignore_case = TRUE))]
  }

  if(length(sub_hit) == 0) {
    # (b) Find substrates in header sequence of TC
    fa_headers <- fasta_header.small[grepl(tc, fasta_header.small, fixed = TRUE)]
    # fa_wors <- str_extract_all(fa_headers, "\\b\\w+\\b")[[1]]
    sub_hit <- sapply(SUBkey, function(key) {
      # Create word-boundary-aware regex pattern (escape special characters)
      pattern <- paste0("\\b", gsub("([\\W])", "\\\\\\1", key), "\\b")
      any(grepl(pattern, fa_headers, ignore.case = TRUE))
    })
    sub_hit <- SUBkey[sub_hit]
  }

  # get alignment hits for TC
  hits <- blasthits[grepl(paste0(gsub("\\.","\\\\.",tc),"$"), qseqid)][order(-bitscore)]
  hit_good <- hits[,any(bitscore >= bitcutoff)]

  # get metabolite IDs (gapseq database) for substances
  exid_all <- unique(SUBid[V1 %in% tolower(sub_hit)])

  # get reaction IDs for transporters
  for(exid in exid_all) {
    exmetID <- gsub("^EX_|_e.$","",exid)
    candidates <- allDB[type == type.str & grepl(exmetID, exmet), id]
    if(length(candidates) > 0) {

    } else {
      DBmissing <- c(DBmissing, data.t)
    }
  }

}

