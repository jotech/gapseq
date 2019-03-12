#!/usr/bin/Rscript

# first argument: input file
# second argument: output file

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("First argument should be an input file and the second one the output file.")
}
#print(args[1])

list.of.packages <- c("stringr") # cran
list.of.packages.ext <- c("Biostrings") # bioconductor
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
new.packages.ext <- list.of.packages[!(list.of.packages.ext %in% installed.packages()[,"Package"])]
if(length(new.packages.ext)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}



suppressMessages(library(Biostrings))
suppressMessages(library(stringr))
suppressMessages(library(methods))

#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.6.5.3.fasta") # complex1
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.9.3.1.fasta") # complex4
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.2.4.1.fasta") # PDH
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/4.1.1.39.fasta") # rubisco
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/1.2.7.1.fasta") # PFOR
#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/6.6.1.2.fasta") # b12

#seq <- readAAStringSet("~/uni/gapseq/dat/seq/Bacteria/unipac90/2.2.1.6.fasta") # 

seq <- readAAStringSet(args[1])
seq.id <- names(seq)

# duplicated hits in seq.id (247,249,300) with str_extract_all
#hits <- str_extract(seq.id, "chain [1-9]+([A-Z])?|subunit [1-9]+([A-Z])?\\b|subunit [A-Z]+\\b|\\b[A-Z]+ subunit|chain [A-Z]+\\b|(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) subunit")
#hits <- str_extract(seq.id, "(subunit|chain|polypeptide) [1-9]+([A-Z])?\\b|(subunit|chain|polypeptide) [A-Z]+\\b|\\b[A-Z]+ (subunit|chain|polypeptide)|(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) (subunit|chain|polypeptide)")

# TODO:
# - wrongly matches: Respiratory-chain NADH dehydrogenase --> because of 'chain' ...)
# - gene names -> causing too many artefacts...
# - how to handle sub-sub-complexes?

# patterns
com.synonymes <- "(\\bsubunit\\b|\\bchain\\b|\\bpolypeptide\\b|\\bcomponent\\b)"
com.pat1  <- paste0(com.synonymes, " [1-9]+([A-Z])?\\b")
com.pat2  <- paste0(com.synonymes, " [A-Z]+\\b")
com.pat3  <- paste0("\\b[A-Z]+ ", com.synonymes)
com.pat4  <- paste0("(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) ", com.synonymes)
com.pat5  <- paste0(com.synonymes, " (alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)")
com.pat6  <- paste0(com.synonymes, " [A-Z][A-z]+\\b")
com.pat7  <- paste0("(large|small) ", com.synonymes)
com.pat8  <- paste0("(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)-", com.synonymes)
hits <- str_extract(seq.id, paste0(com.pat1,"|",com.pat2,"|",com.pat3,"|",com.pat4,"|",com.pat5,"|",com.pat6,"|",com.pat7,"|",com.pat8))

#str_extract(seq.id[863], paste0(com.pat1,"|",com.pat2,"|",com.pat3,"|",com.pat4,"|",com.pat5))

hits <- gsub("subunit|chain|polypeptide|component", "Subunit", hits)

# change order (alpha subunit => subunit alpha)
hits <- str_replace(hits, "([A-Z]+) Subunit", "Subunit \\1")
hits <- str_replace(hits, "(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma) Subunit", "Subunit \\1")
hits <- str_replace(hits, "(alpha|beta|gamma|delta|epsilon|zeta|eta|theta|iota|kappa|lambda|my|ny|omikron|pi|rho|sigma)-Subunit", "Subunit \\1")
hits <- str_replace(hits, "(large|small) Subunit", "Subunit \\1")
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

#
pattern3 <- c("large", "small")
replace3 <- as.character(1:2)
hits <- str_replace_all(hits, setNames(replace3, pattern3))

# other corrections
hits <- str_replace(hits, "(Subunit [0-9]+)[A-z]", "\\1") # how to handle sub-sub-complexes?


#length(unique(hits))
#table(hits)
#seq.id[which(hits=="Subunit 8")][1:5]
#hits[which(hits=="Subunit 8")][1:5]

# remove subunit hits with low amount of sequences (~false hits)
hits.tab <- table(hits)
#print(hits.tab)
if( dim(hits.tab)>0 & mean(hits.tab) >= 10 ){
  hits.low <- names(hits.tab)[which( hits.tab < 5 )]
  hits[which(hits %in% hits.low)] <- NA
}
#seq.id[which(is.na(hits))][1]

# try to get high quality (i.e. ordered subunits) and remove everything else if coverage is high enough
hits.tab <- table(hits)
hits.sel <- grep("Subunit [0-9]+", names(hits.tab), value=T)
hits.sel.cov <- sum(hits.tab[which(names(hits.tab) %in% hits.sel)])
hits.sel.quali <- hits.sel.cov / length(seq.id)
if( hits.sel.quali >= 0.66 ){
  hits[which(!hits %in% hits.sel)] <- NA
}

#print(table(hits))

names(seq) <- paste(seq.id, hits)
#writeXStringSet(seq, "~/uni/gapseq/dat/seq/Bacteria/unipac90/1.6.5.3.fasta")
writeXStringSet(seq, args[2])
