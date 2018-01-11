# Transporter sequences from tcdb
# match hits to agora db => import

# 1. read and process agora database
#   - get_AgoraDB()
#   - get_TransorterType()
# 2. get blast results
#   - transporter_blast()
#     * if parameter only_rel is set, then only transporter of substances 
#       defined in ~/uni/gapseq/dat/substances.csv will be use => speed up
# 3. add transporter
#   - transporter_add() uses agoraDB and blast results as input and adds mising transport reactions to model file

library(Biostrings)
library(foreach)
library(doParallel)
#install_github("mhahsler/rBLAST")
#install_github("jotech/sybil") # modified sybil with some additional functions
library(rBLAST)
library(data.table)


# wrapper function
# TODO: usage of orgID?
find_transporter <- function(model, fasta){
  
  # get databases ( can be updated with update_db() )
  rel_tr <- read.table("~/uni/gapseq/dat/reltrdb.csv", sep="\t", header=T)
  refdb  <- read.table("~/uni/gapseq/dat/reftrdb.csv", sep="\t", header=T)
  
  # blast
  blastdb <- transporter_blast(fasta, orgID="myorg", rel_transporter=rel_tr)
  
  # add transporter
  model   <- transporter_add(model=model, orgID="myorg", blastdb=blastdb, refdb=refdb)
  
  return(model)
}


# Transporter_blast should have a a column model with model id
# TODO: harmonize default layout for exchange reactions (EX_x(e))
#
# model sybil model
# orgID organism ID needed if blast results are for many organism
# blastdb blast results obtained by transporter_blast()
# refdb agora database used as template for transporter insertion. 
transporter_add <- function(model, orgID, blastdb, refdb, cutoff_bits=50, cutoff_ident=70, verbose=TRUE){
  # TODO: check for organism
  
  refdb <- data.table(refdb)
  class <- get_TransporterClasses()
 
  # add exchange reaction if not available
  # TODO: only add exchange if transporter is found ?!?
  exex     <- findExchReact(model)@react_id
  foundex  <- unlist(str_split(unique(blastdb$exid), ";"))
  missedex <- setdiff(foundex, exex)
  if(length(missedex)>0){
    met <- as.character(refdb$exmet[match(missedex, refdb$exid)])
    na <- which(!missedex %in% refdb$exid) # if diffusion reaction is not found
    if(length(na)>0){
      if(verbose) cat("No reference found for:", paste(missedex[na],collapse=","),"\n")
      met <- met[-na]
    }
    if(verbose) cat("Add exchange reactions for:", paste(met,collapse=","),"\n")
    lb <- rep(-1000, length(met))
    ids <- paste0("EX_", gsub("\\[","(",gsub("\\]",")",met))) # set default exchange layout
    model <- addMultiReact(model, ids=ids, mets=met, Scoefs=rep(-1, length(met)), lb=lb)
  }
  
  # add diffusion reactions 
  # TODO: should diffusion reactions be hard-coded?
  exex     <- findExchReact(model)@react_id
  diffex     <- c("EX_o2(e)"="O2t", "EX_h2(e)"="H2td", "EX_h2o(e)"="H2Ot", "EX_etoh(e)"="ETOHt") # echange reaction + corresponding diffusion reaction
  misseddiff <- setdiff(names(diffex), exex)
  if(length(misseddiff)>0){
    met <- as.character(refdb$exmet[match(misseddiff, refdb$exid)])
    na <- which(!misseddiff %in% refdb$exid) # if diffusion reaction is not found
    if(length(na)>0){
      if(verbose) cat("No reference found for:", paste(misseddiff[na],collapse=","),"\n")
      met <- met[-na]
    }
    if(verbose) cat("Add diffusion reactions for:", paste(met,collapse=","),"\n")
    lb    <- rep(-1000, length(met))
    ids   <- paste0("EX_", gsub("\\[","(",gsub("\\]",")",met))) # set default exchange layout
    model <- addMultiReact(model, ids=ids, mets=met, Scoefs=rep(-1, length(met)), lb=lb)
  }
  # add pseudo transporter for compounds with diffusion 
  misseddiff <- setdiff(unname(diffex), model@react_id)
  hit   <- match(misseddiff, refdb$tr)
  idx   <- hit[which(!is.na(hit))]
  if(length(idx)>0) {
      ids   <- as.character(refdb$tr[idx])
      rstr  <- as.character(refdb$rea[idx])
      rname <- as.character(refdb$name[idx])
      if(verbose) cat("Add pseudo transporter for compounds with diffusion:", paste(ids,collapse=","),"\n")
      model <- addReactByString(model, ids=ids, react_str=rstr, reactName=rname)
  }

  
  blastdb <- data.table(blastdb)[org==orgID & Bits>=cutoff_bits & Perc.Ident>=cutoff_ident]
  if(nrow(blastdb)==0){
    warning(paste("No significant transporter hit found for", model@mod_desc))
    return(model)
  }
  

  # summary blast hits per external metabolite
  uniMets <- unique(blastdb$exname)
  dfMets <- data.frame()
  for(i in seq_along(uniMets)){
    met <- uniMets[i]
    trtypes <- paste(unique(blastdb[exname==met]$name), sep=",")
    tc <- unique(blastdb[exname==met]$tc)
    trtype <-  names(class)[as.integer(str_extract(tc, "."))] # get transporter type by first digit of TC number
    dfMets <- rbind(dfMets, data.frame(met=met, table(trtype)))
  }
  if(verbose) print(dfMets)
  
  
  # try to associate hits with agora transporter
  uniHit <- unique(blastdb$tc)
  idx <- match(uniHit, refdb$tc)
  class <- get_TransporterClasses()
  dfhits <- data.frame()
  for(i in seq_along(uniHit)){
    hit <- uniHit[i]
    ref <- which(blastdb$tc==hit)
    ex <- as.character(unique(blastdb$exid[ref]))
    trname <- unique(blastdb$name[ref])
    met <- unique(blastdb$exname[ref])
    nrtype <- as.integer(str_extract(hit, "."))
    if(!nrtype %in% seq_along(class)) next # TC number is beyond scope of this transporter adding script (e.g. 9.B.142.2.9 ~ modification of LPS)
    trtype <- names(class)[nrtype] # get transporter type by first digit of TC number
    agorahit=NA
    if(!is.na(idx[i])){ # found direct link
      agorahit <- refdb$tr[idx[i]]
    }else{
      if(trtype == names(class)[1]){ # if transporter is channel or pore add a diffusion-like transport reaction
        hit3 <- refdb[exid==ex & transporter=="diffusion-like"]
        if(nrow(hit3)>0) agorahit <- hit3$tr else agorahit=NA
      }
      if(is.na(agorahit)){ # try to get agora transporter by tranporter type matching
        hit2 <- refdb[exid==ex & type==trtype]
        if(nrow(hit2)==1){
          agorahit <- hit2$tr
        }else if(nrow(hit2)>1){
          agorahit <- paste(hit2$tr, sep=",")
        }
      }
    }
    dfhits <- rbind(dfhits, data.frame(met=met, tc=hit, type=trtype, rid=agorahit))
  }
  # all unique hits incl. no reference hits
  #dfhits[-duplicated(dfhits$rid),]
  
  # missed metabolites (transporter sequence found for metabolite but no hit in agora db)
  missed <- setdiff(unlist(str_split(dfMets$met, ";")), unlist(str_split(dfhits$met[!is.na(dfhits$rid)], ";")))
  if(verbose & length(missed)>0) cat("Missed metabolites:", paste(missed,collapse=","), "\n")
  
  # reactions to be added
  r2add <- setdiff(unique(dfhits$rid[!is.na(dfhits$rid)]), react_id(model)) # found reactions
  if(verbose & length(r2add)>0) cat("Reactions which will be added:\n", paste(r2add,collapse=","),"\n")
  
  # add reactions
  hit <- match(r2add, refdb$tr)
  if(length(is.na(hit))>0) cat("Not found in reference database:", paste(r2add[which(is.na(hit))]), "\n")
  idx <- which(!is.na(hit))
  model <- addReactByString(model, ids=r2add[idx], react_str=as.character(refdb$rea[hit[idx]]))

  return(model)  
}

update_db <- function(){
  rel_sub  <- rel_substances()
  rel_tr   <- rel_transporter(rel_sub)
  write.table(rel_tr, "~/uni/gapseq/dat/reltrdb.csv", row.names = FALSE, sep="\t")
  
  ref_models   <- readRDS("~/uni/dat/mod/r/agora_contrained_western_anaerob_SWcorrected.RDS")
  refdb    <- get_RefDB(ref_models, rel_sub)
  refdb <- get_TransorterType(refdb)
  write.table(refdb, "~/uni/gapseq/dat/reftrdb.csv", row.names = FALSE, sep="\t")
}


#
# TODO: if only_rel=F then no exchange ids will be set
#
transporter_blast <- function(fasta, orgID, rel_transporter){
  query <- Biostrings::readAAStringSet("~/uni/gapseq/dat/tcdb.fasta")
  
  # limit sequences to ones which have known export reactions
  idx <- match(rel_transporter$fullname, names(query))
  query <- query[idx]
  
  sink("/dev/null") # win: sink("NUL")
  setwd(tempdir())
  makeblastdb(fasta, dbtype = "nucl", args="-out mydb")
  blastdb <- blast(db="./mydb", type="tblastn")
  sink() 
  
  cores=detectCores(); cl <- makeCluster(cores); registerDoParallel(cl)
  distJ <- distJob(length(query), cores)
  
  print(system.time(df_transport <- foreach(i=1:cores,.combine=rbind) %dopar% {
    library(rBLAST) # hast to be reloaded to make 'predict' work in parallel mode
    jobs <- distJ[[i]]
    predict(blastdb, query[jobs])
  }))
  stopCluster(cl)
  
  df_transport <- merge(df_transport, rel_transporter, by.x="QueryID", by.y="id")
  df_transport$org <- orgID
  
  return(df_transport)
}


#string mining for transporter classes in agora db
get_TransorterType <- function(ref_tr, rm_skip=T, rm_unknown=F, rm_notfound=T){
  
  preDat <- read.csv("~/uni/gapseq/dat/vmh2tc.csv", sep="\t", header=T, stringsAsFactors=F)
  newC <- data.frame("transporter"=vector("character", length=nrow(ref_tr)), "type"=vector("character", length=nrow(ref_tr)), "tc"=vector("character", length=nrow(ref_tr)), stringsAsFactors=F)
  hit <- match(ref_tr$tr, preDat$rea)
  idx <- !is.na(hit)
  if(sum(idx)>0){
    newC$transporter[idx] <- preDat$type[hit[idx]]
    newC$tc[idx]          <- preDat$tc[hit[idx]]
  }
  
  class <- get_TransporterClasses()
  for(i in seq_along(class)){
    hit   <- str_extract(tolower(ref_tr$name), paste(tolower(class[[i]]), collapse = "|"))
    if(all(is.na(hit))) next
    Cname <- ifelse(!is.na(hit), names(class)[i], NA)
    newC$transporter <- ifelse(newC$transporter=="", ifelse(!is.na(hit), hit, ""),    newC$transporter)
    newC$type        <- ifelse(newC$type=="",        ifelse(!is.na(Cname),Cname, ""), newC$type       )
  }
  ref_tr <- cbind(ref_tr, newC)
  
  # try to find transporter from reaction name string
  # TODO: assumes compartent to be labeled with [x]
  for(i in which(newC$transporter=="")){
    rea_str <- ref_tr$rea[i]
    met     <- ref_tr$exmet[i]
    hit     <- gregexpr(paste0("\\b",gsub("\\[.\\]","",met),"\\b"), rea_str)[[1]] # search for metabolite without compartment id
    if(length(hit) == 2 & ref_tr$mets[i] == 2) ref_tr$transporter[i] <- "diffusion-like"
    else if(length(hit) == 2) ref_tr$transporter[i] <- "unknown"
  }
  
  skip    <- which(ref_tr$type=="skip")
  notf    <- which(ref_tr$transporter=="")
  unknown <- which(ref_tr$transporter=="unknown")
  #ref_tr[notf,]
  #ref_tr[unknown,]
  
  rmRows <- c()
  if(rm_skip    & length(skip)    > 0) rmRows <- c(rmRows, skip)    # remove entries which seems to be unrelevant
  if(rm_unknown & length(unknown) > 0) rmRows <- c(rmRows, unknown) # remove entries which have unknown transporter type
  if(rm_notfound& length(notf)    > 0) rmRows <- c(rmRows, notf)    # remove entries which have not hit at all (i.e. unknown status)
  if(length(rmRows)>0) ref_tr <- ref_tr[-rmRows,]
  return(ref_tr)
}


# read reference database from model files
# TODO: exmet wrong if transport reactions has more then one relevant external metabolite (e.g. NOr1mq: no[e], n2o[e] not just water!)
get_RefDB <- function(ref_models, rel_sub){
  #ref_models   <- readRDS("~/uni/dat/mod/r/agora_contrained_western_anaerob_SWcorrected.RDS")
  rel_ex <- unique(rel_sub$exid)
  
  ref_tr <- data.frame("tr"=character(),"exid"=character(), "exmet"=character(), "name"=character(), "rea"=character(), "mets"=character())
  for(mod in ref_models){
    ex   <- findExchReact(mod)
    sel  <- intersect(ex@react_id, rel_ex)
    met  <- ex[which(ex@react_id %in% sel),]@met_pos
    hit  <- which(Matrix::colSums(S(mod)[met,])!=0)
    tr  <- react_id(mod)[hit]
    sel2<- which(!tr %in% sel & !tr %in% ref_tr$tr)
    if(length(sel2)==0) next
    exmet<-sapply(hit[sel2], function(i){
      met_id(mod)[met][which(S(mod)[met,i] != 0)][1] # TODO could be more then one match!!
    })
    reaN<- react_name(mod)[hit][sel2]
    reaS<- printReaction(mod, react=tr[sel2], printOut=F)
    reaS<- gsub("^.*\\\t","",reaS)
    if(length(hit[sel2])>1) reaM<- Matrix::colSums(abs(S(mod)[,hit[sel2]]))
    else reaM<- sum(abs(S(mod)[,hit[sel2]]))
    exid <- ex@react_id[match(exmet, ex@met_id)]
    ref_tr <- rbind(ref_tr, data.frame("tr"=tr[sel2],"exid"=exid,"exmet"=exmet, "name"=reaN, "rea"=reaS, "mets"=reaM))
  }
  
  ref_tr <- rbind(ref_tr, data.frame(tr="FORt",exid="EX_for(e)",exmet="for[e]", name="formate transporter via permease", rea="for[e] <==> for[c]", mets="2"))
  
  return(ref_tr)
}


#
# Helper functions
#


# get relevant substances
rel_substances <- function(){
  subs <- read.csv("~/uni/gapseq/dat/sub2pwy.csv", header=T)
  idx1 <- which(subs$exid != "")
  idx2 <- which(subs$altname != "")
  subsdf <- rbind(data.frame(name=subs$name[idx1], exid=subs$exid[idx1]), 
                  data.frame(name=subs$altname[idx2], exid=subs$exid[idx2]))
  return(subsdf)
}


# extract relevant transporter from tcdb given a list of substances to focus on
rel_transporter <- function(rel_sub){
  seq    <- Biostrings::readAAStringSet("~/uni/dat/seq/transporter/tcdb.fasta")
  seqlong <- seq@ranges@NAMES
  seqshort <- gsub( " .*$", "", seqlong)
  seqname  <- sub(".+? ", "", seqlong); seqname  <- gsub(" OS=.*$","",gsub(" \\(EC .*$","", gsub(" \\[.*\\]","",seqname)))
  tcnr <- str_extract(seqshort, "[0-9]+\\.[A-Z]+\\.[0-9]+\\.[0-9]+\\.[0-9]+")
  
  # find transporters for relevant substances
  
  seqExId   <- vector(mode="character", length=length(seq))
  seqExName <- vector(mode="character", length=length(seq))
  for(i in 1:nrow(rel_sub)){
    sub  <- as.character(rel_sub$name[i])
    exid <- as.character(rel_sub$exid[i])
    idx <- grep(paste0("\\b",sub,"\\b"), seqname, ignore.case = T)
    idx <- idx[which(regexpr(exid, seqExId[idx], fixed=T) == -1)] # remove hits which have already entries
    if(length(idx)==0){
      next # if nothing found or already found
    } 
    seqExId[idx]   <- ifelse(seqExId[idx]==""  , exid, paste(seqExId[idx]  , exid, sep=";"))
    seqExName[idx] <- ifelse(seqExName[idx]=="", sub , paste(seqExName[idx], sub , sep=";"))
  }
  idx <- which(seqExId != "")
  
  seqdf <- data.table(id=seqshort[idx], name=seqname[idx], tc=tcnr[idx], exid=seqExId[idx], exname=seqExName[idx], fullname=seqlong[idx])
  return(seqdf)
}


# Transporter classes according to TC classification system
get_TransporterClasses <- function(){
  return(list("1.Channels and pores" = c("channel", "pore"),
              "2.Electrochemical potential-driven transporters" = c("uniport", "symport", "antiport", "permease", "gradient"),
              "3.Primary active transporters" = c("ABC", "ATPase", "ATP"),
              "4.Group translocators" = c("PTS"),
              "skip" = c("extracellular", "degradation")))
}

# group jobs for cpu cores
distJob <- function(Njobs, Nworker){
  tot <- 1:Njobs
  bins <- seq(from = 1,to = Njobs, by=Njobs/Nworker)
  idx <- findInterval(tot, bins)
  distJ <- lapply(1:Nworker, function(i){
    tot[which(idx==i)]})
  return(distJ)
}

