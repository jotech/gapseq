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
  met <- met_id(model)[which(Matrix::rowSums(S(model)[,idx, drop=F]) < 0)] # get metabolites from biomass reactions
  ex  <- findExchReact(model)
  idx <- match(gsub("\\[.*\\]","",met), gsub("\\[.*\\]","",ex@met_id), nomatch=0) # remove comparment tag
  
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
  print(paste(length(aux), "auxotrophies found"))
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
  
  # TODO: supress annoying warnings
  fva   <- sybil::fluxVar(model, percentage=10, react = ex_id, verboseMode=-1)
  
  ex_min<- minSol(fva, "lp_obj")
  ex_cs <- ifelse(ex_min < 0, TRUE, FALSE)
  usage <- sapply(seq_along(idx), function(i){
    if(is.na(idx[i])) NA
    else ex_cs[match(cs[i], ex_id)]})
  df_cs <- data.frame("csource"=cs, "usage"=usage)
  return(df_cs)
}


# check for essential substances needed to growth
getEssentials <- function(model, limit=10){
  ex <- sybil::findExchReact(model)  
  var_r <- sybil::fluxVar(model, react=ex@react_id, percentage=limit)
  ex_max <- sybil::maxSol(var_r, "lp_obj")
  min_id  <- ex@react_id[which(ex_max<0)]

  return(min_id)
}


#model <- models$myb71
model <- test
ub <- model@uppbnd
ex <- findExchReact(model)
#ess <- min_id
lpobject <- sybil::sysBiolAlg(model, algorithm="fba")


ess <- c(min_id, "EX_glc(e)", "EX_glu_L(e)", "EX_xan(e)", "EX_btn(e)", "EX_urea(e)", "EX_gam(e)", "EX_q8(e)", "EX_o2(e)", "EX_cobalt2(e)",
         "EX_q8(e)", "EX_leu_L(e)", "EX_bhb(e)", "EX_ser_L(e)", "EX_pro_L(e)", "EX_alahis(e)", "EX_no2(e)", "EX_tre(e)", "EX_hxan(e)",
         "EX_metala(e)", "EX_alahis(e)", "EX_tsul(e)", "EX_fum(e)", "EX_fru(e)", "EX_glyc3p(e)", "EX_bhb(e)", "EX_metala(e)")

lb <- model@lowbnd
lb[ex@react_pos] <- 0
ess <- c("EX_cu2(e)","EX_zn2(e)","EX_fe3(e)","EX_mg2(e)","EX_k(e)", "EX_mn2(e)","EX_cobalt2(e)", "EX_cl(e)","EX_so4(e)","EX_ca2(e)",
         "EX_thm(e)", "EX_mqn7(e)","EX_glyc3p(e)","EX_gam(e)", "EX_xan(e)",
         #"EX_ribflv(e)","EX_ile_L(e)","EX_val_L(e)","EX_lys_L(e)","EX_spmd(e)","EX_glytyr(e)","EX_glyphe(e)","EX_btn(e)",
         "EX_ocdca(e)", "EX_lanost(e)"#,
         #"EX_glc(e)", "EX_o2(e)"
         )

#ess <- c()

idx <- match(ess, react_id(model))
lb[idx] <- -100
#model@lowbnd <- lb
#sol <- optimizeProb(model, retOptSol=F)
sol <- optimizeProb(lpobject, react=1:length(lb), ub=ub, lb=lb)
message("state:",sol$stat, "\tr=", sol$obj)

names(sol$fluxes) <- react_id(model)
sol$fluxes[ex@react_pos][which(round(sol$fluxes[ex@react_pos],3)<0)]
sol$fluxes[idx]
#setdiff(names(sol$fluxes[ex@react_pos][which(round(sol$fluxes[ex@react_pos],3)<0)]), ess)

shadow <- getRowsDualGLPK(lpobject@problem@oobj); names(shadow) <- react_id(model)
head(sort(shadow[ex@react_pos]),20)
