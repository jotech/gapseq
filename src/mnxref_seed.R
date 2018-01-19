library(data.table)
library(stringr)

mnxref <- read.table("~/uni/dat/db/mnxref/reac_prop.tsv", header=F, sep="\t", stringsAsFactors=F)
idx <- grep("^seed:", mnxref$V6)
mnxSeed <- mnxref[idx,]
#idx2 <- which(mnxSeed$V4=="true")
#filtered <- mnxSeed[idx2,]; colnames(filtered)=c("MNX_ID", "Equation","Description", "Balance", "EC", "Source")
filtered <- mnxSeed[,c(1,4,5,6)]; colnames(filtered)=c("MNX_ID", "Balance", "EC", "Source")
filtered$Source <- gsub("^seed:","", filtered$Source)

seedDB <- fread("~/uni/dat/db/seed/reactions.tsv", header=T, sep="\t")
idx <- match(filtered$Source, seedDB$id)
hit <- str_match(seedDB$abbreviation[idx], "R[0-9]{5}")
filtered$kegg <- ifelse(is.na(hit), "", as.character(hit))

write.table(filtered, file="~/uni/gapseq/dat/mnxref_seed.tsv", quote=F, row.names=F, sep="\t")
dim(filtered)
head(filtered)
