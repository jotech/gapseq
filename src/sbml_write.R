#!/usr/bin/Rscript

library(methods)
library(sybil, quietly=T)
library(data.table, quietly=T)
library(stringr, quietly=T)
library(sybilSBML, quietly=T)
source("~/uni/sybil/R/addReactByString.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
    stop("Needs file with new reactions to be added...")
}

mod <- readRDS("./model.RDS")
reactions <- scan(args[1], what="character", quiet=T)
path <- args[2]

check_react <- function(react_str, compart=c("[c]","[e]")){
  Cmet <- stringr::str_extract_all(react_str, "\\[.\\]")
  valid <- sapply(Cmet, function(c){
    all(c %in% compart)
  })
  return(valid)
}

addReactions <- function(model, candidates){
  vmh <- read.csv(paste0(path,"/dat/vmh_reactions.csv"), sep=",", stringsAsFactors=F)
  bigg <- data.table::fread(paste0(path,"/dat/bigg_reactions.tbl"), sep="\t", stringsAsFactors=F)
  
  modRea <- gsub("_.$","",react_id(model))
  newRea <- setdiff(candidates, modRea)
  
  idx <- match(newRea, vmh$abbreviation)
  if(!all(is.na(idx))){
    react_ids <- vmh$abbreviation[idx]
    react_str <- vmh$formula[idx]
    react_names <- vmh$description[idx]
  }else {react_ids<-c(); react_str<-c(); react_names<-c()}
  idx2 <- match(newRea[which(is.na(idx))], bigg$bigg_id) # get reaction from other database if not found
  if(!all(is.na(idx2))){
    react_ids <- c(react_ids, bigg$bigg_id[idx2])
    bigg_str <- gsub("_(.)\\b","\\[\\1\\]",bigg$reaction_string[idx2])
    react_str <- c(react_str, bigg_str)
    react_names <- c(react_names, bigg$name[idx2])
  }
  ok <- check_react(react_str) # remove reactions from other compartments
  if(length(ok>0)){
    test = addReactByString(model, ids=react_ids[ok], react_str=react_str[ok], reactName=react_names[ok])
    cat("Reactions added: ", setdiff(test@react_id, model@react_id),"\n")
    return(test)
  }
  else{
    warning("No reaction added to model!")
    return(model)
  }
}

cat("Already in the model: ", intersect(mod@react_id,reactions),"\n")
newmod <- addReactions(mod, reactions)
react_attr(newmod)[,1] <- factor(react_attr(newmod)[,1]) # needed to get factor in data.table otherwise export will fail
saveRDS(newmod, "./newmod.RDS")
cat(sybilSBML::writeSBML(newmod, filename="modelnew.xml"), "\n\n")
