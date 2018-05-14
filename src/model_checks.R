check_biosynthesis <- function(model, medium=NA, verbose=TRUE, useNames=FALSE, db="seed", cutoff=1e-6){
  if( sybil::SYBIL_SETTINGS("SOLVER") == "cplexAPI" ) solver_ok=1 else if( sybil::SYBIL_SETTINGS("SOLVER") == "glpkAPI" ) solver_ok=5
  
  if( is.na(medium) )
    medium <- paste0("EX_",c('cpd00001','cpd00007','cpd00009','cpd00027','cpd00076','cpd00030','cpd00034','cpd00048','cpd00058','cpd00063','cpd00067','cpd00099','cpd00149','cpd00205','cpd00254','cpd00531','cpd00971','cpd01012','cpd01048','cpd10515','cpd10516','cpd11595','cpd00013'), "_e0")
  model <- set_diet(model, diet = medium)
  if( db=="seed" ){
    aa <- c("cpd00035[c0]", "cpd00051[c0]", "cpd00132[c0]", "cpd00041[c0]", "cpd00084[c0]", "cpd00023[c0]", "cpd00033[c0]", "cpd00119[c0]", "cpd00322[c0]", "cpd00107[c0]", "cpd00039[c0]", "cpd00060[c0]", "cpd00064[c0]", "cpd00066[c0]", "cpd00129[c0]", "cpd00054[c0]", "cpd00161[c0]", "cpd00065[c0]", "cpd00069[c0]", "cpd00156[c0]", "cpd00266[c0]", "cpd00053[c0]")
    nb <- c("cpd00128[c0]", "cpd00207[c0]", "cpd00307[c0]", "cpd00151[c0]", "cpd00092[c0]", "cpd00246[c0]", "cpd00309[c0]")
    nt <- c("cpd00002[c0]", "cpd00052[c0]", "cpd00038[c0]", "cpd00056[c0]", "cpd00062[c0]")
    co <- c("cpd00305[c0]", "cpd00220[c0]", "cpd00016[c0]", "cpd00003[c0]", "cpd00644[c0]", "cpd00104[c0]", "cpd00393[c0]", "cpd00166[c0]", "cpd00028[c0]")
  } else stop("Only seed db support yet.")
  df_met <- rbind(data.frame(met=aa, typ="Amino acids"), data.frame(met=nb, typ="Nucleic bases"), data.frame(met=nt, typ="Nucleoside triphosphate"), data.frame(met=co, typ="Cofactors"))
  
  df <- data.frame()
  if( verbose ) cat("Checking for production of compounds:", nrow(df_met), "\n")
  for(i in 1:nrow(df_met)){
    if( verbose ) cat("\r",i,"/",nrow(df_met))
    metId  <- df_met$met[i]
    metIdx <- match(metId, model@met_id)
    if( is.na(metIdx) )
      next
    metName<- met_name(model)[metIdx]
    metTyp <- df_met$typ[i]
    tmpMod <- addReact(model, id="tmpObj", met=metId, Scoef=-1, lb=0, ub=1000)
    tmpMod <- changeObjFunc(tmpMod, react="tmpObj")
    sol    <- optimizeProb(tmpMod, retOptSol=FALSE)
    status <- sol$stat==solver_ok & round(sol$obj,cutoff) > 0
    df <- rbind(df, data.frame(id=metId, name=metName, typ=metTyp, producable=status))  
  }
  return(df)
}

