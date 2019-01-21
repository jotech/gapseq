#!/usr/bin/Rscript

# first argument: input file
# second argument: output file

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("First argument should be an input file and the second one the output file.")
}
#print(args[1])

list.of.packages <- c("Biostrings", "stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(Biostrings))
suppressMessages(library(stringr))
suppressMessages(library(methods))

#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.6.5.3.fasta")
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.9.3.1.fasta")
seq <- readAAStringSet(args[1])
seq.id <- names(seq)

# duplicated hits in seq.id (247,249,300) with str_extract_all
#hits <- str_extract(seq.id, "chain [1-9]+([A-Z])?|subunit [1-9]+([A-Z])?\\b|subunit [A-Z]+\\b|\\b[A-Z]+ subunit|chain [A-Z]+\\b|(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) subunit")
hits <- str_extract(seq.id, "(subunit|chain|polypeptide) [1-9]+([A-Z])?\\b|(subunit|chain|polypeptide) [A-Z]+\\b|\\b[A-Z]+ (subunit|chain|polypeptide)|(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) (subunit|chain|polypeptide)")
hits <- gsub("subunit|chain|polypeptide", "Subunit", hits)

# change order (alpha subunit => subunit alpha)
hits <- str_replace(hits, "([A-Z]+) Subunit", "Subunit \\1")
hits <- str_replace(hits, "(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) Subunit", "Subunit \\1")
#hits[1949]

# latin numbers
latin <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV")
pattern0 <- paste0('(?i)\\b',rev(latin))
replace0 <- as.character(15:1)
hits <- str_replace_all(hits, setNames(replace0, pattern0))

# letter numbers
pattern1 <- paste0('(?i)\\b',LETTERS[1:26], '\\b')
replace1 <- as.character(1:26)
hits <- str_replace_all(hits, setNames(replace1, pattern1))

# greek numbers
greek <- c("alpha","beta","gamma","delta","epsilon","zeta","eta","theta","iota","kappa","lambda","my","ny","omikron","pi","rho","sigma")
pattern2 <- paste0('(?i)\\b',greek, '\\b')
replace2 <- as.character(1:17)
hits <- str_replace_all(hits, setNames(replace2, pattern2))
#hits[1949]

# other corrections
hits <- str_replace(hits, "(Subunit [0-9]+)[A-Z]", "\\1")


#length(unique(hits))
#table(hits)
#seq.id[which(hits=="Subunit 12")]

# remove subunit hits with low amount of sequences (~false hits)
hits.tab <- table(hits)
#print(hits.tab)
if( mean(hits.tab) >= 10 ){
  hits.low <- names(hits.tab)[which( hits.tab < 5 )]
  hits[which(hits %in% hits.low)] <- NA
}
#seq.id[which(is.na(hits))][1]

#print(table(hits))

names(seq) <- paste(seq.id, hits)
#writeXStringSet(seq, "~/uni/gapseq/dat/seq/Bacteria/unipac90/1.6.5.3.fasta")
writeXStringSet(seq, args[2])
