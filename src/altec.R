library(data.table)
library(stringr)

html <- readLines("~/uni/gapseq/dat/EC-Supplements-24-2018.html")

ec_hit <- str_extract_all(html, "[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+")
ec_hit.flat <- sapply(ec_hit, paste)
keyword_hit <- str_extract(html, "transferred")

altec <- data.table(ec.numbers=ec_hit.flat, ec.count=sapply(ec_hit, length), keyword=keyword_hit)

altec[ !is.na(keyword) & ec.count>1 ]

