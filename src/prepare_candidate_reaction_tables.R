prepare_candidate_reaction_tables <- function(blast.res, transporter.res, high.evi.rxn.BS, min.bs.for.core, curve.alpha = 1) {
  # Read reaction blast results
  require(data.table)
  
  dt <- fread(blast.res, header=T, stringsAsFactors = F)
  dt <- dt[,.(rxn, name, ec, tc = NA_character_, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send, pathway, status, 
              pathway.status, seed = dbhit, complex, exception, complex.status)]
  
  # Read transporter blast results
  dt.trans <- fread(transporter.res, header=T, stringsAsFactors = F)
  dt.trans <- dt.trans[,.(rxn = id, name = paste("transport",tc,sub,sep="-"), ec = NA_character_, tc, qseqid, pident, evalue, bitscore, qcovs, stitle, 
                          sstart, send, pathway = NA_character_, status = NA_character_, pathway.status = NA_character_, seed = rea, 
                          complex = NA_character_, exception = 0, complex.status = NA_integer_)]
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
  dt.cand <- copy(dt)
  
  #
  # 1. Handling reactions associated with complexes (i.e. multi-mer enzymes)
  #
  dt.cand.clpx <- copy(dt.cand[!is.na(complex)]) 
  dt.cand.clpx[is.na(bitscore), bitscore := 0] # no blast hit -> NA -> 0 BS
  dt.cand.clpx[, max.bs := max(bitscore, na.rm=TRUE), by = c("seed","complex")]
  dt.cand.clpx <- dt.cand.clpx[max.bs == bitscore] # get only the best hits
  dt.cand.clpx <- dt.cand.clpx[!duplicated(paste(seed,complex, sep = "$"))]
  dt.cand.clpx[, max.bs := NULL]
  # filter out "Subunit undefined" rows if they map to a range of a defined subunit of the enzyme
  dt.cand.clpx <- dt.cand.clpx[order(stitle,seed,complex,-bitscore)] # TODO: Make sure the subunits with notation "subunit undefined" is last row per seed-reaction id "seed"
  dt.cand.clpx[, rm := F]
  dt.cand.clpx$itmp <- 1:nrow(dt.cand.clpx)
  # Calculating overlaps
  for(i in 1:nrow(dt.cand.clpx)) {
    astart <- dt.cand.clpx[i, sstart]
    aend <- dt.cand.clpx[i, send]
    astitle <- dt.cand.clpx[i, stitle]
    aseed <- dt.cand.clpx[i, seed]
    dt_g_tmp <- dt.cand.clpx[i < itmp & rm==F & stitle == astitle & seed == aseed & !is.na(sstart)]
    if(nrow(dt_g_tmp)>0) {
      dt_g_tmp[calc_seq_overlap(astart, aend, sstart, send) > 0.5, rm := T]
      ind.rm <- dt_g_tmp[rm==T, itmp]
      dt.cand.clpx[itmp %in% ind.rm, rm := T]
    }
  }
  dt.cand.clpx <- dt.cand.clpx[rm == F]
  dt.cand.clpx[, rm := NULL]
  dt.cand.clpx[, itmp := NULL]
  dt.cand.clpx[complex != "Subunit undefined" & is.na(complex.status), min.bs := min(bitscore, na.rm = T), by = "seed"] # to be conservative, choose the lowest bitscore for complexes which were not predicted to be present
  dt.cand.clpx <- dt.cand.clpx[bitscore == min.bs | is.na(min.bs)]
  dt.cand.clpx[, min.bs := NULL]
  # for complexes whose presence was predicted choose the subunit with the hightest bitscore as reference.
  dt.cand.clpx[complex != "Subunit undefined" & complex.status == 1, max.bs := max(bitscore, na.rm = T), by = "seed"]
  dt.cand.clpx <- dt.cand.clpx[bitscore == max.bs | is.na(max.bs)]
  dt.cand.clpx[, max.bs := NULL]
  # for remaining undefined subunits choose the best hit and handle it as monomer
  dt.cand.clpx[complex == "Subunit undefined", max.bs := max(bitscore, na.rm = T), by = "seed"]
  dt.cand.clpx <- dt.cand.clpx[bitscore == max.bs | is.na(max.bs)]
  dt.cand.clpx[, max.bs := NULL]
  
  # 2. Handling of reactions which have pathway topology evidences
  dt.cand.topo <- copy(dt.cand[status %in% c("bad_blast","no_blast") & pathway.status %in% c("full","treshold","keyenzyme")])
  dt.cand.topo[, bitscore := high.evi.rxn.BS]
  dt.cand.topo <- dt.cand.topo[!duplicated(seed)]
  
  # 3. Handling of single reactions associated to monomer enzymes
  dt.cand.mono <- copy(dt.cand[is.na(complex)])
  dt.cand.mono[, all.na := all(is.na(bitscore)), by = c("seed")]
  dt.cand.mono <- dt.cand.mono[all.na == FALSE] # remove reactions & complex subunits which have no bitscore at all
  dt.cand.mono[, all.na := NULL]
  dt.cand.mono[, max.bs := max(bitscore, na.rm=TRUE), by = "seed"]
  dt.cand.mono <- dt.cand.mono[max.bs == bitscore]
  dt.cand.mono[, max.bs := NULL]
  dt.cand.mono <- dt.cand.mono[!duplicated(seed)]
  
  # 3.5. Handling of spontaneous reactions
  dt.cand.spon <- copy(dt[status=="spontaneous" & is.na(bitscore)])
  dt.cand.spon[, bitscore := high.evi.rxn.BS]
  dt.cand.spon <- dt.cand.spon[!duplicated(seed)]
  
  # 4. Combine all three parts and choose highest support
  dt.cand <- rbindlist(list(dt.cand.clpx, dt.cand.mono, dt.cand.topo, dt.cand.spon))
  dt.cand[, max.bs := max(bitscore, na.rm=TRUE), by = "seed"]
  dt.cand <- dt.cand[max.bs == bitscore]
  dt.cand[, max.bs := NULL]
  dt.cand <- dt.cand[!duplicated(seed)]
  
  # 5. transfortm best bitscores to weights
  dummy.weight = 100
  dt.cand[, bs.tmp := bitscore]
  dt.cand[bs.tmp > high.evi.rxn.BS, bs.tmp := high.evi.rxn.BS] # if bitscore is higher than the bitscore threshold then assign a weight, that woulb be close to 0. 
  dt.cand[, weight := 1 - ((bs.tmp - min.bs.for.core) / (high.evi.rxn.BS - min.bs.for.core))] # assign normalised weight
  dt.cand[weight > 1, weight := 1]
  dt.cand[, weight := weight^curve.alpha]
  dt.cand[, weight := weight * dummy.weight]
  dt.cand[, weight := (weight + 0.005) / (1 + 0.005)]
  dt.cand[, bs.tmp := NULL]

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
