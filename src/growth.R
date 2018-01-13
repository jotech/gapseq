library(sybil)

findCarbonSource <- function(mod){
  library(igraph)
  #data("Ec_core")
  #mod <- Ec_core
  myb71_med <- c("EX_cu2(e)","EX_zn2(e)","EX_fe3(e)","EX_pi(e)","EX_mg2(e)","EX_k(e)", "EX_mn2(e)","EX_cobalt2(e)", "EX_cl(e)","EX_so4(e)","EX_ca2(e)",
                 "EX_glc(e)", "EX_o2(e)",
                 "EX_thm(e)", "EX_mqn7(e)",
                 #"EX_pyr(e)", 
                 "EX_gam(e)", 
                 "EX_xan(e)", "EX_lanost(e)")
  mod <- set_diet(myb71_new, myb71_med)
  check_growth(mod)
  
  sol <- optimizeProb(mod, retOptSol=F)
  tol <- 1e-03
  nzf <- which(abs(sol$fluxes)>tol) # non zero fluxes
  nM  <- length(mod@met_id)
  ignore <- mod@met_id[tail(sort(rowSums(as.matrix(abs(S(mod)))), index.return=T)$ix,10)]
  mBio <- which(S(mod)[,findrBiomass(mod)]!=0)
  ex   <- findExchReact(mod)
  mEx  <- ex@met_pos; fEx  <- round(sol$fluxes[ex@react_pos], -log10(tol))
  
  g <- make_empty_graph(n = nM, directed = TRUE)
  g <- set_vertex_attr(g,name="id",value=mod@met_id)
  g <- set_vertex_attr(g,name="type",value="internal")
  g <- set_vertex_attr(g,name="type",index=mBio,value="biomass")
  g <- set_vertex_attr(g,name="type",index=mEx[which(fEx<0)],value="uptake")
  g <- set_vertex_attr(g,name="type",index=mEx[which(fEx>0)],value="production")
  for(i in nzf){
    idx <- which(S(mod)[,i]!=0)
    sig <- ifelse(sol$fluxes[i]>0,1,-1) # negative fluxes swap educt and product
    idxE <- which(sig*S(mod)[,i] < 0) # educts
    idxP <- which(sig*S(mod)[,i] > 0) # products
    if(length(idxE)==0 | length(idxP)==0) next
    comb <- expand.grid(idxE,idxP) # combinatorics
    for(j in 1:nrow(comb)){
      edu <- comb$Var1[j]
      pro <- comb$Var2[j]
      g <- add_edges(g, edges=c(edu,pro), weight=abs(sol$fluxes[i]), reaction=mod@react_id[i])
    }
    
  }
  
  g <- delete.vertices(g, which(vertex_attr(g, "id") %in% ignore)) # delete some substances (h2o, h, ...)
  g <- delete.vertices(g, which(degree(g) == 0)) # delete isolates
  
  idxL <- as.numeric(V(g)); names(idxL) <- vertex_attr(g,"id")
  labV <- ifelse(vertex_attr(g, "type")=="internal", NA, vertex_attr(g, "id"))
  sizeV<- ifelse(vertex_attr(g, "type")=="internal", 0, 15)
  colV <- ifelse(!vertex_attr(g, "type") %in% c("production","uptake"), "yellow", ifelse(vertex_attr(g, "type")=="production", "green", "red"))
  plot(g, vertex.label=labV, vertex.size=sizeV, vertex.label.cex=0.6, vertex.color=colV, edge.arrow.size=0.5)
  
  # influence of uptake substances
  vEx  <- V(g)[vertex_attr(g, "type") == "uptake"]
  vBio <- V(g)[vertex_attr(g, "type") == "biomass"]
  dist <- distances(g, v=vEx, to = vBio, mode="out")
  rownames(dist) <- vertex_attr(g, "id")[vEx]; colnames(dist) <- vertex_attr(g, "id")[vBio]
  dist[dist==Inf] <- 0; dist <- round(dist,1); dist
  sort(rowSums(dist, na.rm = T))
  sort(rowMeans(dist, na.rm = T))
  dist[which(rownames(dist) %in% c("glc_D[e]", "gam[e]")),]
  
  # path from uptake substance to biomass
  idxB <- idxL["gam6p[c]"]
  idxS <- idxL["glc_D[e]"]
  path <- all_shortest_paths(g, from=idxS, to=idxB, mode="out")
  vertex_attr(g, "id")[path$res[[1]]]
  printReaction(mod,react=edge_attr(g, "reaction")[path$res[[1]]])
  
  # neighborhood of a substance
  idxL[neighbors(g, idxL["glc_D[e]"], mode="out")]
  incident(g, idxL["glc_D[e]"])
  edge_attr(g,"reaction", index = incident(g, idxL["glc_D[e]"]))
  edge_attr(g,"weight", index = incident(g, idxL["glc_D[e]"]))
  
  # plot neighborhood of a substance
  gsub <- subgraph.edges(g, eids=incident(g, idxL["glc_D[e]"]))
  plot(gsub, 
       edge.label=paste(edge_attr(gsub,"reaction"),round(edge_attr(gsub,"weight"),1)), edge.label.cex=0.7,
       vertex.label=vertex_attr(gsub,"id"))
}


set_diet <- function(model, diet, lb=-1000){
  ex <- findExchReact(model)
  lowbnd(model)[ex@react_pos] <- 0
  idx <- match(diet, ex@react_id)
  lowbnd(model)[ex@react_pos[idx]] <- lb
  return(model)
}


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

auxotrophy <- function(model, checkNoEx=FALSE, rmMetals=FALSE, useNames=FALSE){
  metals <- c("EX_so4(e)", "EX_ca2(e)", "EX_cu2(e)", "EX_fe2(e)", "EX_zn2(e)", "EX_fe3(e)", "EX_mg2(e)", "EX_k(e)", "EX_cl(e)", "EX_mn2(e)", "EX_cobalt2(e)",
              "EX_cpd00048_e0", "EX_cpd00063_e0", "EX_cpd10516_e0", "EX_cpd00058_e0", "EX_cpd00034_e0", "EX_cpd00254_e0", "EX_cpd00205_e0", "EX_cpd00099_e0","EX_cpd00030_e0")
  idx <- findrBiomass(model)
  met <- met_id(model)[which(Matrix::rowSums(S(model)[,idx, drop=F]) < 0)] # get metabolites from biomass reactions
  ex  <- findExchReact(model)
  idx <- match(gsub("\\[.*\\]","",met), gsub("\\[.*\\]","",ex@met_id), nomatch=0) # remove comparment tag
  
  mex   <- ex@react_id[idx] # components from biomass which have exchange reaction => test with fva
  mname <- react_name(model)[ex@react_pos[idx]]
  fva   <- fluxVar(model, percentage=10, react = mex) 
  max   <- maxSol(fva, "lp_obj")
  if (useNames) aux   <- mname[which(max < 0)] else aux <- mex[which(max < 0)]
  
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
  if(rmMetals) aux <- setdiff(aux, metals)
  print(paste(length(aux), "auxotrophies found"))
  return(sort(aux))
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
getEssentials <- function(model, ignore=NULL, limit=10, useNames=F){
  ex <- sybil::findExchReact(model)  
  var_r <- sybil::fluxVar(model, react=ex@react_id, percentage=limit)
  ex_max <- sybil::maxSol(var_r, "lp_obj")
  
  if(useNames) min_id <- sort(model@react_name[ex@react_pos[which(ex_max<0)]])
    else min_id  <- sort(ex@react_id[which(ex_max<0)])
  
  print(paste(length(min_id), "essential substances found"))
  if(length(ignore)==0) return(min_id)
  else return(setdiff(min_id, ignore))
}


trimMedium <- function(model, medium=NULL, limit=10, useNames=F){
  ex <- sybil::findExchReact(model)
  if(is.null(medium)) medium <- ex@react_id
  med <- intersect(ex@react_id, medium)
  lb <- model@lowbnd
  lb[ex@react_pos] <- 0
  idx <- match(med, react_id(model))
  lb[idx] <- -100
  model@lowbnd <- lb
  
  var_r <- sybil::fluxVar(model, react=med, percentage=limit)
  plot(var_r)
  if( all(round(var_r@lp_obj,3)==0) ) warning("all fva solutions are almost zero..")
  ex_max <- sybil::maxSol(var_r, "lp_obj")
  
  not_essential  <- med[which(!round(ex_max,6)<0)]
  if(useNames) not_essential <- model@react_name[match(not_essential, model@react_id)]
  
  return(c(not_essential, setdiff(medium,ex@react_id)))
}

trimMediumRand <- function(model, medium=NULL){
  ex <- sybil::findExchReact(model)  
  if(is.null(medium)) medium <- ex@react_id
  med <- intersect(ex@react_id, medium)
  lb <- model@lowbnd
  ub <- model@uppbnd
  lpobject <- sybil::sysBiolAlg(model, algorithm="fba")
  
  medium_new <- med
  for(i in 1:100){
    lb[ex@react_pos] <- 0
    idx   <- match(medium_new, react_id(model))
    rmIdx <- sample(1:length(idx),1)
    rmMet <- medium_new[rmIdx]
    idx   <- idx[-rmIdx]
    lb[idx] <- -100
    sol <- optimizeProb(lpobject, react=1:length(lb), ub=ub, lb=lb)
    if(round(sol$obj,6)>0){
      medium_new <- setdiff(medium_new, rmMet)
    }
  }
  return(c(setdiff(med, medium_new), setdiff(medium,ex@react_id)))
}


check_growth <- function(model, medium=NULL, uptake=-100, printShadow=F, printExchanges=F, returnSol=F){
  ex <- findExchReact(model)
  if(is.null(medium)){
    medium <- ex@react_id[which(ex@lowbnd < 0)]
    uptake <- ex@lowbnd[which(ex@lowbnd < 0)]
  } 
  lpobject <- sybil::sysBiolAlg(model, algorithm="fba")
  ub <- model@uppbnd
  lb <- model@lowbnd
  
  #ess <- c(min_id, "EX_glc(e)", "EX_glu_L(e)", "EX_xan(e)", "EX_btn(e)", "EX_urea(e)", "EX_gam(e)", "EX_q8(e)", "EX_o2(e)", "EX_cobalt2(e)",
  #         "EX_q8(e)", "EX_leu_L(e)", "EX_bhb(e)", "EX_ser_L(e)", "EX_pro_L(e)", "EX_alahis(e)", "EX_no2(e)", "EX_tre(e)", "EX_hxan(e)",
  #         "EX_metala(e)", "EX_alahis(e)", "EX_tsul(e)", "EX_fum(e)", "EX_fru(e)", "EX_glyc3p(e)", "EX_bhb(e)", "EX_metala(e)")

  lb[ex@react_pos] <- 0
  idx <- match(medium, react_id(model))
  lb[idx] <- uptake
  sol <- optimizeProb(lpobject, react=1:length(lb), ub=ub, lb=lb)
  message("state:",sol$stat, "\tr=", round(sol$obj,6))
  
  if(printExchanges){
    names(sol$fluxes) <- react_id(model)
    sol$fluxes[ex@react_pos][which(round(sol$fluxes[ex@react_pos],3)<0)]
    print(sort(round(sol$fluxes[idx],1)))
  }
  
  if(printShadow){ # TODO: not working length(shadow) != length(react_id) wtf??
    shadow <- getRowsDualGLPK(lpobject@problem@oobj); names(shadow) <- react_id(model)
    print(head(sort(shadow[ex@react_pos]),20))
  }
  if(returnSol) return(sol)
}

id2names <- function(mod, ids){
  idx <- match(ids, mod@react_id)
  return(mod@react_name[idx])
}
