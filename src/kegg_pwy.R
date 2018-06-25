library(data.table)

kegg1 <- fread("~/uni/dat/db/kegg/ec_pathway.map", header=F, col.names = c("ec", "pwy"))
kegg2 <- fread("~/uni/dat/db/kegg/pathway", header=F, col.names = c("pwy", "name"))

kegg.pwy <- data.table()
for(p in kegg2$pwy){
  pec <- gsub("ec:","",paste0(kegg1[pwy==p,ec],collapse = ","))
  pid <- gsub("path:","",p)
  pname <- kegg2[pwy==p,name]
  prxn <- gsub("ec:","rxn-",paste0(kegg1[pwy==p,ec],collapse = ","))
  kegg.pwy <- rbind(kegg.pwy, data.table(id=pid, name=pname, altname="", hierarchy="kegg", taxrange="", reaId=prxn, reaEc=pec, keyRea=""))
}

fwrite(kegg.pwy, file="~/uni/gapseq/dat/kegg_pwy.tbl", sep="\t", quote=F)
