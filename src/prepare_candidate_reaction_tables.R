prepare_candidate_reaction_tables <- function(blast.res, transporter.res, high.evi.rxn.BS, min.bs.for.core, curve.alpha = 1) {
  # Read reaction blast results
  require(data.table)
  
  dt <- fread(blast.res, header=T, stringsAsFactors = F, blank.lines.skip = T)
  
  dt <- dt[pathway != "|PWY-6168|"] # pathway exists only in fungi
  dt <- dt[pathway != "|PWY-1042|"] # pathway exists only in plants
  dt <- dt[pathway != "|PWY-6970|"] # pathway exists only in plants
  dt <- dt[pathway != "|PWY-4983|"] # pathway exists only in mammals
  dt <- dt[pathway != "|PWY-6845|"] # pathway exists only in mammals
  
  dt <- dt[,.(rxn, name, ec, tc = NA_character_, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send, pathway, status, 
              pathway.status, seed = dbhit, complex, exception, complex.status)]
  
  # Read transporter blast results
  dt.trans <- fread(transporter.res, header=T, stringsAsFactors = F, blank.lines.skip = T)
  dt.trans <- dt.trans[,.(rxn = id, name = paste("transport",tc,sub,sep="-"), ec = NA_character_, tc, qseqid, pident, evalue, bitscore, qcovs, stitle, 
                          sstart, send, pathway = NA_character_, status = NA_character_, pathway.status = NA_character_, seed = rea, 
                          complex = NA_character_, exception = 0, complex.status = NA_integer_)]
  dt.trans[bitscore >= high.evi.rxn.BS, status := "good_blast"]
  dt.trans[bitscore <  high.evi.rxn.BS, status := "bad_blast"]
  
  # # A specific fix for the issue with reactions 1.3.8.1 and 1.3.8.13
  # if("1.3.8.1" %in% dt$ec & "1.3.8.13" %in% dt$ec) {
  #   rm.ids <- c()
  #   one.hits.id <- which(dt$ec == "1.3.8.1" & !is.na(dt$bitscore))
  #   two.hits.id <- which(dt$ec == "1.3.8.13" & !is.na(dt$bitscore))
  #   if(length(one.hits.id)>0 & length(two.hits.id)>0) {
  #     for(i in one.hits.id) {
  #       for(j in two.hits.id) {
  #         ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
  #         if(ol > 0.2) {
  #           if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, j)
  #           if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, i)
  #         }
  #       }
  #     }
  #   }
  #   rm.ids <- unique(rm.ids)
  #   if(length(rm.ids > 0))
  #     dt <- dt[-rm.ids]
  # }
  # 
  # # A specific fix for the issue with reactions 4.1.2.9 and 4.1.2.22
  # if("4.1.2.9" %in% dt$ec & "4.1.2.22" %in% dt$ec) {
  #   rm.ids <- c()
  #   one.hits.id <- which(dt$ec == "4.1.2.9" & !is.na(dt$bitscore))
  #   two.hits.id <- which(dt$ec == "4.1.2.22" & !is.na(dt$bitscore))
  #   if(length(one.hits.id)>0 & length(two.hits.id)>0) {
  #     for(i in one.hits.id) {
  #       for(j in two.hits.id) {
  #         ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
  #         if(ol > 0.2) {
  #           if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, j)
  #           if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, i)
  #         }
  #       }
  #     }
  #   }
  #   rm.ids <- unique(rm.ids)
  #   if(length(rm.ids > 0))
  #     dt <- dt[-rm.ids]
  # }
  # 
  # # A specific fix for the issue with reactions 2.6.1.13 and 2.6.1.11
  # if("2.6.1.13" %in% dt$ec & "2.6.1.11" %in% dt$ec) {
  #   rm.ids <- c()
  #   one.hits.id <- which(dt$ec == "2.6.1.13" & !is.na(dt$bitscore))
  #   two.hits.id <- which(dt$ec == "2.6.1.11" & !is.na(dt$bitscore))
  #   if(length(one.hits.id)>0 & length(two.hits.id)>0) {
  #     for(i in one.hits.id) {
  #       for(j in two.hits.id) {
  #         ol <- calc_seq_overlap(dt[i,sstart],dt[i,send],dt[j,sstart],dt[j,send])
  #         if(ol > 0.2) {
  #           if(dt[i, bitscore] > max(dt[two.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, j)
  #           if(dt[j, bitscore] > max(dt[one.hits.id, bitscore]))
  #             rm.ids <- c(rm.ids, i)
  #         }
  #       }
  #     }
  #   }
  #   rm.ids <- unique(rm.ids)
  #   if(length(rm.ids > 0))
  #     dt <- dt[-rm.ids]
  # }
  # 
  
  # This function checks if an gene was assigned to two different ec numbers, which however catalyse different reactions that are
  # usually catalysed by individual enzymes. Thus, the EC asigment with the lower bitscore is dismissed.
  resolve_common_EC_conflicts <- function(ec1, ec2, dt) {
    if(ec1 %in% dt$ec & ec2 %in% dt$ec) {
      rm.ids <- c()
      one.hits.id <- which(dt$ec == ec1 & !is.na(dt$bitscore))
      two.hits.id <- which(dt$ec == ec2 & !is.na(dt$bitscore))
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
    return(dt)
  }
  
  resolve_common_TC_conflicts <- function(tc1, tc2, dt) {
    if(tc1 %in% dt$tc & tc2 %in% dt$tc) {
      rm.ids <- c()
      one.hits.id <- which(dt$tc == tc1 & !is.na(dt$bitscore))
      two.hits.id <- which(dt$tc == tc2 & !is.na(dt$bitscore))
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
    return(dt)
  }
  
  # specific reaction conflict fixes
  dt <- resolve_common_EC_conflicts("1.3.8.1","1.3.8.13", dt)
  dt <- resolve_common_EC_conflicts("4.1.2.9","4.1.2.22", dt)
  dt <- resolve_common_EC_conflicts("2.6.1.13","2.6.1.11", dt)
  dt <- resolve_common_EC_conflicts("2.3.1.29","2.3.1.37", dt)
  dt <- resolve_common_EC_conflicts("2.3.1.29","2.3.1.50", dt)
  dt <- resolve_common_EC_conflicts("2.3.1.37","2.3.1.50", dt)
  dt <- resolve_common_EC_conflicts("1.17.3.2","1.17.1.4", dt) # xanthine oxidase vs xanthine dehydrogenase
  dt <- resolve_common_EC_conflicts("1.3.3.6","1.3.1.8", dt) # acyl-CoA oxidase vs acyl-CoA dehydrogenase
  dt <- resolve_common_EC_conflicts("1.2.7.1","1.2.1.51", dt) # NADP-dependent Pyruvate dehydrogenase vs FMN-dependent PDH
  dt <- resolve_common_EC_conflicts("2.6.1.11","2.6.1.19", dt) # acetylornithine transaminase VS 4-aminobutyrate—2-oxoglutarate transaminase
  dt <- resolve_common_EC_conflicts("1.1.1.100","1.1.1.36", dt) # 3-oxoacyl-ACP reductase VS acetoacetyl-CoA reductase
  dt <- resolve_common_EC_conflicts("4.1.3.36","4.2.1.150", dt) # 1,4-dihydroxy-2-naphthoyl-CoA synthase VS (S)-3-hydroxybutanoyl-CoA dehydrogenase
  dt <- resolve_common_EC_conflicts("1.2.1.76","1.2.1.10", dt) # succinate semialdehyde dehydrogenase VS acetaldehyde dehydrogenase
  
  # specific transporter conflict fixes
  dt.trans <- resolve_common_TC_conflicts("1.a.8.2.1","2.a.14.1.3", dt.trans)
  dt.trans <- resolve_common_TC_conflicts("1.a.8.2.7","2.a.14.1.3", dt.trans)
  
  # Due to BRENDA's alternative ECs theres a mismatch of metacyc reactions to seed reaction for EC 2.6.1.36 and EC 2.6.1.13 ... remove the mismatches.
  dt <- dt[!(rxn == "L-LYSINE-AMINOTRANSFERASE-RXN" & grepl("rxn00467|rxn20496|rxn33315", seed))]
  
  if(dt[!duplicated(rxn) & rxn == "NAD-SYNTH-GLN-RXN",status] %in% c("no_blast","bad_blast") &
      dt[!duplicated(rxn) & rxn == "NAD-SYNTH-NH3-RXN",status] == "good_blast") {
    dt[rxn == "NAD-SYNTH-GLN-RXN", pathway.status := NA_character_]
  }
  if(dt[!duplicated(rxn) & rxn == "NAD-SYNTH-NH3-RXN",status] %in% c("no_blast","bad_blast") &
     dt[!duplicated(rxn) & rxn == "NAD-SYNTH-GLN-RXN",status] == "good_blast") {
    dt[rxn == "NAD-SYNTH-NH3-RXN", pathway.status := NA_character_]
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
  dummy.weight <- 100
  dt.cand[, bs.tmp := bitscore]
  dt.cand[bs.tmp > high.evi.rxn.BS, bs.tmp := high.evi.rxn.BS] # if bitscore is higher than the bitscore threshold then assign a weight, that woulb be close to 0. 
  dt.cand[, weight := (bs.tmp - high.evi.rxn.BS) * ((0.005 - dummy.weight)/(high.evi.rxn.BS - min.bs.for.core)) + 0.005]
  dt.cand[weight < 0.005, weight := 0.005]
  dt.cand[weight > dummy.weight, weight := 100]
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
