library(data.table)
library(stringr)

#mnxref <- read.table("~/uni/dat/db/mnxref/reac_prop.tsv", header=F, sep="\t", stringsAsFactors=F)
mnxref <- fread(cmd='grep -v "^#" ~/uni/gapseq/dat/mnxref_reac_xref.tsv', sep="\t", header=F, col.names = c("XREF", "MNXID", "rxn"))

mnxref.seed <- mnxref[XREF %like% "seed" ]; setnames(mnxref.seed, old="XREF", new="seed"); mnxref.seed[,seed:=str_remove(seed, "seed:")]
mnxref.meta <- mnxref[XREF %like% "metacyc" ]; setnames(mnxref.meta, old="XREF", new="other"); mnxref.meta[,other:=str_remove(other, "metacyc:")]
mnxref.kegg <- mnxref[XREF %like% "kegg" ]; setnames(mnxref.kegg, old="XREF", new="other"); mnxref.kegg[,other:=str_remove(other, "kegg:")]
mnxref.bigg <- mnxref[XREF %like% "bigg" ]; setnames(mnxref.bigg, old="XREF", new="bigg"); mnxref.bigg[,bigg:=str_remove(bigg, "bigg:")]

mnxref.dic  <- merge(mnxref.seed, mnxref.meta, by="MNXID")
mnxref.dic  <- rbind(mnxref.dic, merge(mnxref.seed, mnxref.kegg, by="MNXID"))
#fwrite(mnxref.dic[,.(MNXID, seed, other)], "~/uni/gapseq/dat/mnxref_seed-other.tsv", sep = "\t")

mnxref.dic2  <- merge(mnxref.bigg, mnxref.meta, by="MNXID")
mnxref.dic2 <- rbind(mnxref.dic2, merge(mnxref.bigg, mnxref.kegg, by="MNXID"))
#fwrite(mnxref.dic2[,.(MNXID, bigg, other)], "~/uni/gapseq/dat/mnxref_bigg-other.tsv", sep = "\t")
