library(data.table)
library(stringr)

kegg1 <- fread("~/uni/dat/db/kegg/ec_pathway.map", header=F, col.names = c("ec", "pwy"))
kegg2 <- fread("~/uni/dat/db/kegg/pathway", header=F, col.names = c("pwy", "name"))
kegg3 <- fread("~/uni/dat/db/kegg/reaction_pathway.map", header=F, col.names = c("rxn", "pwy"))
kegg4 <- fread("~/uni/dat/db/kegg/reaction_ec.map", header=F, col.names = c("ec", "rxn"))
kegg5 <- fread("~/uni/dat/db/kegg/reaction", header=F, col.names = c("rxn", "rxn.name"))

kegg.pwy <- data.table()
for(p in kegg2$pwy){
  pid <- gsub("path:","",p)
  pname <- kegg2[pwy==p,name]
  #prxn <- gsub("ec:","rxn-",paste0(kegg1[pwy==p,ec],collapse = ","))
  prxn <- gsub("rn:","",paste0(kegg3[pwy==p,rxn],collapse = ","))
  prxn.id <- paste0("rn:",unlist(str_split(prxn,",")))
  #pec <- gsub("ec:","",paste0(kegg1[pwy==p,ec],collapse = ","))
  pec  <- gsub("ec:","",kegg4[match(prxn.id, rxn), ec])
  pec  <- ifelse(is.na(pec),"", pec)
  prxn.name <- str_extract(kegg5[match(prxn.id, rxn), rxn.name], ".*?(?=;)")
  
  kegg.pwy <- rbind(kegg.pwy, data.table(id=pid, name=pname, altname="", hierarchy="kegg", taxrange="", reaId=prxn, reaEc=paste(pec,collapse=","), keyRea="", reaName=paste0(prxn.name,collapse=";"), reaNr=length(pec), ecNr=length(na.omit(pec)), superpathway=FALSE, status=TRUE))
}

fwrite(kegg.pwy, file="~/uni/gapseq/dat/kegg_pwy.tbl", sep="\t", quote=F)
