library(data.table)
library(stringr)

#mnxref <- read.table("~/uni/dat/db/mnxref/reac_prop.tsv", header=F, sep="\t", stringsAsFactors=F)
mnxref <- fread("~/uni/gapseq/dat/mnxref_reac_xref.tsv", skip=365)


mnxref.seed <- mnxref[`#XREF` %like% "seed" ]; setnames(mnxref.seed, old="#XREF", new="seed"); mnxref.seed[,seed:=str_remove(seed, "seed:")]
mnxref.meta <- mnxref[`#XREF` %like% "metacyc" ]; setnames(mnxref.meta, old="#XREF", new="other"); mnxref.meta[,other:=str_remove(other, "metacyc:")]
mnxref.kegg <- mnxref[`#XREF` %like% "kegg" ]; setnames(mnxref.kegg, old="#XREF", new="other"); mnxref.kegg[,other:=str_remove(other, "kegg:")]

mnxref.dic  <- merge(mnxref.seed, mnxref.meta, by="MNX_ID", )
mnxref.dic  <- rbind(mnxref.dic, merge(mnxref.seed, mnxref.kegg, by="MNX_ID", ))

fwrite(mnxref.dic, "~/uni/gapseq/dat/mnxref_seed-other.tsv", sep = "\t")
