generate_rxn_stoich_hash <- function(stoich,rev) {
  out <- c()
  for(i in 1:length(stoich)) {
    rxn.info  <- str_split(unlist(str_split(string = stoich[i],pattern = ";")), pattern = ":", simplify = T)
    if(rxn.info[1,1]!="") {
    
    met.ids   <- paste0(rxn.info[,2],"[",rxn.info[,3],"]")
    met.scoef <- as.integer(rxn.info[,1])
    s.order   <- order(met.ids)
    met.ids   <- met.ids[s.order]
    met.scoef <- met.scoef[s.order]
    if(rev[i]=="=" & met.scoef[1]<0)
      met.scoef <- met.scoef * -1
    if(rev[i]=="<") {
      met.scoef <- met.scoef * -1
      rev[i]==">"
    }
    
    st.out <- paste(paste(met.scoef,met.ids,sep=":"), collapse = ";")
    st.out <- paste(st.out, rev[i], sep="#")
    out <- c(out,st.out)
    } else {
      out <- c(out, "ERR")
    }
  }
  return(out)
}



