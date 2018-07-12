library(data.table)
library(stringr)

get_hierarchy <- function(kegg6){
  hierarchy <- data.table()
  for(i in seq_along(kegg6)){
    lin <- kegg6[i]
    if ( startsWith(lin, "A") ) sys <- lin
    if ( startsWith(lin, "B") ) sub <- lin
    if ( startsWith(lin, "C") ) {
      pwy <- lin
      hierarchy <- rbind(hierarchy, data.table(sys,sub,pwy))
    }
  }
  hierarchy$pwy <- paste0("map",str_extract(hierarchy$pwy, "[0-9]{5}"))
  hierarchy$sys <- str_remove_all(hierarchy$sys, "(A|<b>|<\\/b>)")
  hierarchy$sub <- str_remove_all(hierarchy$sub, "^B[:space:]+")
  return(hierarchy)
}

kegg1 <- fread("~/uni/dat/db/kegg/ec_pathway.map", header=F, col.names = c("ec", "pwy"))
kegg2 <- fread("~/uni/dat/db/kegg/pathway", header=F, col.names = c("pwy", "name"))
kegg3 <- fread("~/uni/dat/db/kegg/reaction_pathway.map", header=F, col.names = c("rxn", "pwy"))
kegg4 <- fread("~/uni/dat/db/kegg/reaction_ec.map", header=F, col.names = c("ec", "rxn"))
kegg5 <- fread("~/uni/dat/db/kegg/reaction", header=F, col.names = c("rxn", "rxn.name"))
kegg6 <- get_hierarchy(scan("~/uni/dat/db/kegg/br08901.keg", what = "character", sep = "\n", quiet = T))


kegg.pwy <- data.table()
#for(p in head(kegg2$pwy,10)){
for(p in kegg2$pwy){
  pid <- gsub("path:","",p)
  pname <- kegg2[pwy==p,name]
  
  # get reaction/ec using ec-pathway relationship
  pec1 <- gsub("ec:","", kegg1[pwy==p,ec])
  if( length(pec1) > 0 ) {
    prxn.id1 <- kegg4[match(paste0("ec:",pec1), ec), rxn]
  }else{
    prxn.id1 <- NULL; pec1 <- NULL
  }
    
  # get reaction/ec using reaction-pathway relationship
  prxn2 <- gsub("rn:","",paste0(kegg3[pwy==p,rxn],collapse = ","))
  if( prxn2 != ""){
    prxn.id2 <- paste0("rn:",unlist(str_split(prxn2,",")))
    pec2 <- gsub("ec:","",kegg4[match(prxn.id2, rxn), ec])
  }else{
    prxn.id2 <- NULL; pec2 <- NULL
  }
  
  # join both
  idx <- !duplicated(c(pec1,pec2))
  prxn.id <- c(prxn.id1,prxn.id2)[idx]
  pec     <- c(pec1,pec2)[idx]
  pec  <- ifelse(is.na(pec),"", pec)
  prxn    <- paste0(gsub("rn:","",prxn.id), collapse = ",")
  
  prxn.name <- str_extract(kegg5[match(prxn.id, rxn), rxn.name], ".*?(?=;)")
  prxn.name <- ifelse(is.na(prxn.name),gsub("rn:","",prxn.id), prxn.name)
  idx  <- match(pid, kegg6$pwy)
  hier <- paste("kegg",kegg6$sys[idx], kegg6$sub[idx], sep = ";") 
  p
  cat(p, pname, length(pec),length(prxn.name),"\n")
  kegg.pwy <- rbind(kegg.pwy, data.table(id=pid, name=pname, altname="", hierarchy=hier, taxrange="", reaId=prxn, reaEc=paste(pec,collapse=","), keyRea="", reaName=paste0(prxn.name,collapse=";"), reaNr=length(pec), ecNr=length(na.omit(pec)), superpathway=FALSE, status=TRUE))
}

fwrite(kegg.pwy, file="~/uni/gapseq/dat/kegg_pwy.tbl", sep="\t", quote=F)

