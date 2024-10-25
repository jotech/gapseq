parse_BMjson <- function(json.file, gs_mets) {
  # gs_mets should be a data.table of the database metabolites
  # json.file should be path to the json-biomass file
  
  # chekcing if required packages are there
  if(any(!(c("cobrar","jsonlite") %in% rownames(installed.packages())))) {
    missing_packages <- c("cobrar","jsonlite")[!(c("cobrar","jsonlite") %in% rownames(installed.packages()))]
    stop(paste0("Following required R-package(s) is/are missing: ",paste(missing_packages, collapse = ", ")))
  }
  suppressMessages(require(cobrar))
  suppressMessages(require(jsonlite))
  suppressMessages(require(data.table))

  # Read json bm file
  bmjs <- read_json(json.file)

  # get name of metabolite groups
  g_names <- unlist(lapply(bmjs$met_groups, function(x) x$group_name))
  names(bmjs$met_groups) <- g_names
  
  # check summed groups mass
  g_smass <- sum(unlist(lapply(bmjs$met_groups, function(x) x$mass)))
  if(abs(1-g_smass) > 0.0010)
    warning(paste0("Summed mass fractions do not add up to 100% (i.e. 1 gDW):\n", 
                   paste(g_names,
                         round(unlist(lapply(bmjs$met_groups, function(x) x$mass)), 
                               digits = 5), 
                         collapse = "\n")))
  
  # check if DNA, RNA and Protein is in Biomass
  if(!all(c("DNA","RNA","Protein") %in% g_names))
    stop("One of DNA, RNA, or Protein is missing in biomass formulation.")
  
  get_groupBM_table <- function(x) {
    bm_i <- x
    #print(bm_i$group_name)
    bm_i_fractions <- sum(unlist(lapply(bm_i$components, function(x) x$coef)))
    if(abs(1-bm_i_fractions) > 0.0001 & bm_i$unit_components == "MOLFRACTION")
      warning(paste0("Fractions of group \"",bm_i$group_name,"\" do not add up to 1 (i.e. 100%). Gapseq will rescale, but please consider cerrecing the biomass definition."))
    
    c_list <- list()
    
    for(i in 1:length(bm_i$components)) {
      bm_i_c <- bm_i$components[[i]]
      tmp_dt <- data.table(id = bm_i_c$id,
                           comp = bm_i_c$comp,
                           coef = bm_i_c$coef / bm_i_fractions)
      if("link" %in% names(bm_i_c)) {
        c_links <- unlist(strsplit(bm_i_c$link, "\\|"))
        c_links_mets <- gsub("\\:.*$", "", c_links)
        c_links_coef <- as.numeric(gsub("^.*\\:", "", c_links))
        tmp_dt2 <- data.table(id   = c_links_mets,
                              comp = bm_i_c$comp,
                              coef = c_links_coef * bm_i_c$coef / bm_i_fractions)
        tmp_dt <- rbind(tmp_dt, tmp_dt2)
      }
      c_list[[i]] <- tmp_dt
    }
    c_list <- rbindlist(c_list)
    c_list <- merge(c_list, gs_mets [, .(id, name, formula)], sort = F)
    c_list[, formula := gsub("R|R[0-9]+","",formula)]
    c_list <- c_list[, .(coef  = sum(coef)), by = c("id","name","comp","formula")]
    c_list[, mass := mass(formula)] 
    
    c_list[, g.m := mass * coef]
    denom <- sum(c_list$g.m)
    c_list[, Scoef := bm_i$mass * coef / denom * 1000]
    
    return(c_list)
  }
  
  grp_tab <- lapply(bmjs$met_groups, get_groupBM_table)
  grp_tab <- rbindlist(grp_tab, idcol = "met_group")
  #grp_tab[, sum(Scoef * mass / 1000)]
  
  # get scoefs for GAM
  gam <- bmjs$energy_GAM
  gam_dt <- data.table(id = c("cpd00001","cpd00002","cpd00008","cpd00009","cpd00067"),
                       Scoef = c(gam, gam, -gam, -gam, -gam),
                       comp = rep("c",5),
                       met_group = "GAM")
  gam_dt <- merge(gam_dt, gs_mets[,.(id, name, formula)])
  grp_tab <- rbind(grp_tab, gam_dt, fill = T)
  
  # sum up scoefs for the same meabolites
  grp_tab <- grp_tab[, .(Scoef = sum(Scoef), met_group = paste(unique(met_group), collapse = ", ")), 
                     by = c("id","name","comp","formula")]
  grp_tab[, Scoef := -Scoef]
  
  # include Biomass-Pseudo-Metabolite
  bmm_dt <- data.table(id    = "cpd11416",
                       name  = "Biomass",
                       comp  = "c",
                       Scoef = 1) 
  grp_tab <- rbind(grp_tab, bmm_dt, fill = T)
  
  # prepare names
  grp_tab[, id := paste0(id,"[",comp,"0]")]
  
  res <- list()
  res$id     <- bmjs$id
  res$name   <- bmjs$name
  res$domain <- bmjs$domain
  res$bmS    <- grp_tab
  
  return(res)
}
