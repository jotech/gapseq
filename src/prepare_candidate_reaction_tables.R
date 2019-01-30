prepare_candidate_reaction_tables <- function(blast.res, transporter.res, high.evi.rxn.BS, for.GS.draft = FALSE) {
  # Read reaction blast results
  require(data.table)
  
  dt <- fread(blast.res, header=T, stringsAsFactors = F)
  dt <- dt[,.(rxn, name, ec, tc = NA_character_, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send, pathway, status, pathway.status, seed = dbhit)]
  
  # Read transporter blast results
  dt.trans <- fread(transporter.res, header=T, stringsAsFactors = F)
  dt.trans <- dt.trans[,.(rxn = id, name = paste("transport",tc,sub,sep="-"), ec = NA_character_, tc, qseqid, pident, evalue, bitscore, qcovs, stitle, 
                          sstart, send, pathway = NA_character_, status = NA_character_, pathway.status = NA_character_, seed = rea)]
  dt.trans[bitscore >= high.evi.rxn.BS, status := "good_blast"]
  dt.trans[bitscore <  high.evi.rxn.BS, status := "bad_blast"]
  
  # A specific fix for the issue with reactions 1.3.8.1 and 1.3.8.13
  if("1.3.8.1" %in% dt$ec & "1.3.8.13" %in% dt$ec) {
    rm.ids <- c()
    one.hits.id <- which(dt$ec == "1.3.8.1" & !is.na(dt$bitscore))
    two.hits.id <- which(dt$ec == "1.3.8.13" & !is.na(dt$bitscore))
    if(length(one.hits.id)>0 & length(two.hits.id)>0) {
      for(i in one.hits.id) {
        for(j in two.hits.id) {
          ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
          if(ol > 0.2) {
            if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
              rm.ids <- c(rm.ids, j)
            if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
              rm.ids <- c(rm.ids, i)
          }
        }
      }
    }
    rm.ids <- unique(rm.ids)
    if(length(rm.ids > 0))
      dt <- dt[-rm.ids]
  }
  
  # A specific fix for the issue with reactions 4.1.2.9 and 4.1.2.22
  if("4.1.2.9" %in% dt$ec & "4.1.2.22" %in% dt$ec) {
    rm.ids <- c()
    one.hits.id <- which(dt$ec == "4.1.2.9" & !is.na(dt$bitscore))
    two.hits.id <- which(dt$ec == "4.1.2.22" & !is.na(dt$bitscore))
    if(length(one.hits.id)>0 & length(two.hits.id)>0) {
      for(i in one.hits.id) {
        for(j in two.hits.id) {
          ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
          if(ol > 0.2) {
            if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
              rm.ids <- c(rm.ids, j)
            if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
              rm.ids <- c(rm.ids, i)
          }
        }
      }
    }
    rm.ids <- unique(rm.ids)
    if(length(rm.ids > 0))
      dt <- dt[-rm.ids]
  }
  
  # prepare reaction/transporter blast results table
  dt <- splitcol2rows_mget(dt, "seed", " ")
  dt.trans <- splitcol2rows_mget(dt.trans, "seed", ",")
  dt <- rbind(dt, dt.trans)
  
  #
  # construct gapfill candidate reactions and their weight
  #
  if(!for.GS.draft) {
    dt.cand <- copy(dt)
  } else {
    dt.cand <- copy(dt[bitscore < high.evi.rxn.BS | (bitscore >= high.evi.rxn.BS & status == "bad_blast")])
    dt.cand[bitscore > high.evi.rxn.BS, bitscore := bitscore / 2] # apply penalty to the bitscore of exception reactions
  }
  dt.cand[, max.bs := max(bitscore), by = "seed"]
  dt.cand <- dt.cand[max.bs == bitscore]
  dt.cand <- dt.cand[!duplicated(seed)]
  dt.cand[, max.bs := NULL]
  dt.cand[bitscore >= high.evi.rxn.BS, bitscore := high.evi.rxn.BS - 1]
  dt.cand[, weight := 1 - (bitscore / high.evi.rxn.BS)]
  
  return(list(dt = dt, dt.cand = dt.cand))
}

# a little helper function to calculate the overlap of two genetic regions
calc_seq_overlap <- function(astart, aend, bstart, bend) {
  out <- numeric(length = length(bstart))
  for(i in seq_along(bstart)){
    out[i] <- length(intersect(astart:aend,bstart[i]:bend[i]))/((length(astart:aend)+length(bstart[i]:bend[i]))/2)
  }
  return(out)
}

splitcol2rows_mget <- function(dtInput, col2split, sep){
  dtInput <- dtInput[, .(tmp.add.col = unlist(strsplit(get(col2split),sep,T))), by=names(dtInput)]
  
  dtInput[, c(col2split):=NULL];
  setnames(dtInput, 'tmp.add.col', col2split); 
  return(dtInput);
}