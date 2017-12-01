library(sybil)


findrBiomass <- function(model, keys=c("biom")){
  ex_pos <- sybil::findExchReact(model)@react_pos
  rbio <- vector()
  for(k in keys){
    idx <- grep(k, sybil::react_id(model), ignore.case = TRUE)
    if(length(idx)==0) idx <- grep(k, sybil::react_name(model), ignore.case = TRUE)
    if(length(idx)>0) rbio <- c(rbio, idx)
  }
  if(length(rbio)==0) return(NULL)
  rbio <- setdiff(rbio, ex_pos) # exclude exchange reactions
  #return(sybil::react_id(model)[rbio])
  return(rbio)
}

auxotrophy <- function(model, checkNoEx=FALSE){
  idx <- findrBiomass(model)
  met <- met_id(model)[which(rowSums(S(model)[,idx]) < 0)] # get metabolites from biomass reactions
  ex  <- findExchReact(model)
  idx <- match(gsub("\\[.\\]","",met), gsub("\\[.\\]","",ex@met_id), nomatch=0)
  
  mex <- ex@react_id[idx] # components from biomass which have exchange reaction => test with fva
  fva <- fluxVar(model, percentage=10, react = mex) 
  max <- maxSol(fva, "lp_obj")
  aux <- mex[which(max < 0)]
  
  if(checkNoEx){
    noex<- met[is.na(idx)] # components from biomass which have no exchange reaction => test if could be produced by metabolic model
    aux2 <- sapply(noex, function(m){
      midx <- match(m, met_id(model))
      mid  <- met_id(model)[midx]
      rid <- paste0("Test_",mid)
      testmod <- addReact(model, id=rid, met=mid, Scoef=c(-1)) # add dummy reaction and try to optimize flux
      testmod <- changeObjFunc(testmod, react=rid)
      sol     <- optimizeProb(testmod, retOptSol=FALSE)
      if(sol$ok==0 & sol$obj>1e-3) TRUE else FALSE 
      aux <- c(aux, aux2[which(aux2==FALSE)])
    })
  }
  return(aux)
}

cs <- c("EX_glc(e)", "EX_fru(e)", "EX_malt(e)", "EX_malthx(e)", "EX_man(e)", "EX_acgam(e)", "EX_gam(e)", "EX_rib_D(e)", "EX_sucr(e)", "EX_gal(e)", 
"EX_lcts(e)", "EX_cellb(e)", "EX_tre(e)", "EX_rib_D(e)", "EX_xyl_D(e)", "EX_arab_L(e)", "EX_mnl(e)", "EX_melib(e)", "EX_raffin(e)", 
"EX_glcn(e)", "EX_salcn(e)", "EX_glcur(e)", "EX_rmn(e)", "EX_galur(e)", "EX_arab_D(e)", "EX_fuc_L(e)", "EX_acnam(e)", "EX_acgal(e)", 
"EX_lyx_L(e)", "EX_glyc(e)", "EX_cit(e)", "EX_tur(e)", "EX_ch4(e)", "EX_etoh(e)", "EX_chsterol(e)", "EX_xylt(e)")


csource <- function(model, cs=cs){
  cs    <- as.character(cs)
  ex    <- findExchReact(model)
  idx   <- match(cs, ex@react_id)
  
  lowbnd(model)[grep("^EX_", react_id(model))] <- -1000 # set all exchanges to -INF
  ex_id <- cs[which(!is.na(idx))]
  fva   <- sybil::fluxVar(model, percentage=10, react = ex_id, verboseMode=0) 
  ex_min<- minSol(fva, "lp_obj")
  ex_cs <- ifelse(ex_min < 0, TRUE, FALSE)
  usage <- sapply(seq_along(idx), function(i){
    if(is.na(idx[i])) NA
    else ex_cs[match(cs[i], ex_id)]})
  df_cs <- data.frame("csource"=cs, "usage"=usage)
  return(df_cs)
}





exchange <- data.frame()
rast2id <- models_descr[,c("id", "Code")]
for( mod in mods ){
  lowbnd(mod)[grep("^EX_", react_id(mod))] <- -1000 # set all exchanges to -INF
  #lowbnd(mod)[grep("EX_o2", react_id(mod), fixed = T)] <- 0
  ex <- findExchReact(mod)
  cs       <- as.character(csource$Exchange)
  fp       <- as.character(fermprod$Exchange)
  ex_id       <- intersect(ex@react_id, c(cs, fp))
  ex_cs_ind   <- which(ex_id %in% cs)
  ex_fp_ind   <- which(ex_id %in% fp)
  fva         <- fluxVar(mod, percentage=10, react = ex_id) 
  ex_max      <- maxSol(fva, "lp_obj")
  ex_min      <- minSol(fva, "lp_obj")
  ex_cs       <- ifelse(ex_min[ex_cs_ind] < 0, TRUE, FALSE)  
  ex_fp       <- ifelse(ex_max[ex_fp_ind] > 0, TRUE, FALSE)
  all_cs      <- ex_cs[match(cs, ex_id[ex_cs_ind])]
  all_fp      <- ex_fp[match(fp, ex_id[ex_fp_ind])]
  all_cs[is.na(all_cs)] <- FALSE
  all_fp[is.na(all_fp)] <- FALSE
  typ         <- c(rep("carbon source", length(all_cs)), rep("fermentation product", length(all_fp)))
  name        <- c(as.character(csource$Name), as.character(fermprod$Name))
  #rast_id     <- gsub(".sbml","",mod@mod_desc)
  #mod_id      <- as.character(rast2id$Code[which(rast2id$id %in% rast_id)])
  mod_id      <- mod_id(mod)
  exchange    <- rbind(exchange, data.frame("mod"=mod_id, "sub"=c(cs,fp), "name"=name, "ex"=c(all_cs, all_fp), "typ"=typ))
}
