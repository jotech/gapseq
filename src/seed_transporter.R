#
# Search in seed reaction database for transport reactions and identify transporter type and transported substances
#

library(data.table)
library(stringr)


setwd("~/uni/gapseq/")
reaDB <- fread("dat/seed_reactions_corrected.tsv")
reaDB.new <- fread("dat/seed_reactions_new.tsv")
reaEC <- fread("dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv", sep="\t")
metDB <- fread("dat/seed_metabolites.tsv")

#test <- reaDB[is_transport==1]
test <- merge(reaDB[is_transport==1], reaDB.new[,.(id,aliases)], by="id", all.x=T)

#
# Find type of transporter in reaction database
#

patterns <- list("1.Channels and pores" = c("channel", "pore"),
     "2.Electrochemical potential-driven transporters" = c("uniport", "symport", "antiport", "permease", "gradient"),
     "3.Primary active transporters" = c("ABC", "ATPase", "ATP"),
     "4.Group translocators" = c("PTS"),
     "Diffusion" = c("diffusion"),
     "skip" = c("extracellular", "degradation", "reductase"))

test[grep(paste(patterns[[1]], collapse = "|"), name),type1:=names(patterns)[1]]
test[grep(paste(patterns[[2]], collapse = "|"), name),type1:=paste0(type1,";",names(patterns)[2])]
test[grep(paste(patterns[[3]], collapse = "|"), name),type1:=paste0(type1,";",names(patterns)[3])]
test[grep(paste(patterns[[4]], collapse = "|"), name),type1:=paste0(type1,";",names(patterns)[4])]
test[grep(paste(patterns[[5]], collapse = "|"), name),type1:=paste0(type1,";",names(patterns)[5])]

test[grep(paste(patterns[[1]], collapse = "|"), aliases.y),type2:=names(patterns)[1]]
test[grep(paste(patterns[[2]], collapse = "|"), aliases.y),type2:=paste0(type2,";",names(patterns)[2])]
test[grep(paste(patterns[[3]], collapse = "|"), aliases.y),type2:=paste0(type2,";",names(patterns)[3])]
test[grep(paste(patterns[[4]], collapse = "|"), aliases.y),type2:=paste0(type2,";",names(patterns)[4])]
test[grep(paste(patterns[[5]], collapse = "|"), aliases.y),type2:=paste0(type2,";",names(patterns)[5])]

test[,type:=ifelse(is.na(type1), type2, type1)]
test$type <- gsub("NA;|;NA","", test$type)

# duplicate rows for reactions with multiple types
type.split <- str_split(test$type, ";")
type.split.ln <- sapply(type.split, length)
test.new <- data.table()
for( i in which(type.split.ln>1) ){
  tmp <- test[i, cbind(unlist(type.split[i]), .SD)] # duplicate rows
  tmp[,type:=V1]
  test.new <- rbind(test.new, tmp[,-"V1"])
}
test <- rbind(test[-which(type.split.ln>1),], test.new)
test[id=="rxn05581", .(id, type)]


# entry stoichiometry field: stoichiometry_coeff:seed_id:compartment:x:full_name
test[stoichiometry %like% "ATP" & stoichiometry %like% "ADP" & is.na(type), type:=names(patterns)[3]] # atpases
test[stoichiometry %like% "cpd00067:0" & stoichiometry %like% "cpd00067:1" & is.na(type), type:=names(patterns)[2]] # proton antiport
test[stoichiometry %like% "cpd00971:0" & stoichiometry %like% "cpd00971:1"& is.na(type), type:=names(patterns)[2]] # sodium antiport
test[stoichiometry %like% "Phosphoenolpyruvate" & stoichiometry %like% "Pyruvate" & is.na(type), type:=names(patterns)[4]] # pts system

Nmets <- str_count(test$stoichiometry, "cpd")
test[Nmets == 2 & is.na(test$type), type:=names(patterns)[5]] # diffusion reactions

# debug: 208/5062 transporter without transporter type
# 20190801: 234/5280 transporter without type
test[is.na(type)]
sum(is.na(test$type))
dim(test)
table(test$type)





#
# Find substance which is transported
#

mets  <- str_extract_all(test[,stoichiometry], "cpd[0-9]+")
exmet <- sapply(mets, function(m){paste(m[duplicated(m)],collapse=";")})
test$exmet <-  exmet # substances that change compartments

conmet <- gsub(":1","",str_extract_all(test[exmet=="" & type==names(patterns)[4]]$stoichiometry, "cpd[0-9]+:1"))
test[exmet=="" & type==names(patterns)[4], exmet:=conmet] # pts substances

# if more then one transported metabolite is detected remove obviouse ones (h,na,water,...)
lex <- sapply(str_split(test[,exmet],";"),length)
test[lex>=2 & exmet %like% "cpd00067", exmet:=gsub("cpd00067","",exmet)] # h+
test[lex>=2 & exmet %like% "cpd00971", exmet:=gsub("cpd00971","",exmet)] # na+
test[lex>=2 & exmet %like% "cpd00002", exmet:=gsub("cpd00002","",exmet)] # atp
test[lex>=2 & exmet %like% "cpd00001", exmet:=gsub("cpd00001","",exmet)] # h2o
test[lex>=2 & exmet %like% "cpd00010", exmet:=gsub("cpd00010","",exmet)] # coa
test$exmet <- gsub("^;","",gsub(";$","",gsub(";;",";",test$exmet)))

test[lex>=2]

names <- sapply(test$exmet, function(m){
  idx <- match(unlist(str_split(m,";")), metDB$id)
  paste(metDB$name[idx], collapse = ";")
})
test$exmetnames <- names

# debug 114/5062 transporter without transported substance detected
# 20190801: debug 264/5280
test[exmet==""]
nrow(test[exmet==""])
test[exmet=="",name]
# 
lex <- sapply(str_split(test[,exmet],";"),length)
table(lex)
length(which(lex>=2))
test[lex>=2]

#
# export
#
#export <- test[!is.na(type) & exmet!="",.(id,name,type,exmet,exmetnames,stoichiometry,direction)]
export <- test[!is.na(type) & exmet!="",.(id,name,type,exmet,exmetnames)]
export$exmet <- gsub("(cpd[0-9]+)","\\1\\[e0\\]",export$exmet) # add compartment tag
#fwrite(export, file="~/uni/gapseq/dat/seed_transporter.tbl", sep="\t", quote = F)



#
# OLD? tcdb 
#

library(Biostrings)
tcdb <- readAAStringSet("/home/jo/uni/gapseq/dat/tcdb.fasta")
names(tcdb)
allMets <- unique(sort(unlist(str_split(export$exmetnames,";"))))
allMets <- allMets[-grep("\\(", allMets)] # TODO: how to handle substance names with brackets?? (causing problems in grep search)
allMets <- allMets[which(sapply(allMets, str_length) > 2)]
ignore <- allMets
allMets <- 
key <- paste0(allMets, collapse = "|")
key2 <- paste0("\\b",paste0(allMets, collapse = "\\b|\\b"), "\\b")
t <- grep(key, names(tcdb),perl=T,value=F)
t2 <- str_extract(names(tcdb), key)


grep("chloride", names(tcdb),value=T)




