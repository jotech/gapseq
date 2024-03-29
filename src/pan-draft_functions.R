### parse the option and return a list
parse_input_options <- function(opt, suffix) {
  # Check if opt is a directory
  if (dir.exists(opt)) {
    # List files in the directory
    files <- list.files(opt, full.names = TRUE, pattern = paste0(suffix, "$"))
    cat("Reading input from directory:\t", opt, "\n") # print directory

  } else {
    # Treat opt as a comma-separated list of file names or a pattern with wildcard
    if (grepl(",", opt) | length(opt) == 1) {
      cat("Reading input from specified filename pathway\n")
      # Treat opt as a comma-separated list of file names
      files <- unlist(strsplit(opt, ","))
    } else {
      stop("The input has not been interpreted properly\n")
    }

    # Get only files that have the specified suffix
    files <- files[grepl(paste0(suffix, "$"), files)]

    # If some of the files have a different suffix, return an error message
    if (length(files) == 0) {
      stop("No files with the specified suffix. Please check the file names.")
    }
    # Check if files exist
    if (!all(file.exists(files))) {
      stop("Some files do not exist. Please check the file names.")
    }
  }
  cat(files, "\n")
  return(files)
}

### load files (.RDS) from paths  
# specify suffix to save the prefix as id in a list
load_files_from_paths_for_RDS <- function(file_paths, suf) {
  file_list <- list() # Create an empty list to store the loaded files
  # Iterate over the file paths and read the corresponding files
  for (path in file_paths) {
    file_data <- readRDS(path)
    filename <- basename(path)
    filename_without_suffix <- sub(suf, "", filename)
    file_list[[filename_without_suffix]] <- file_data
  }
  return(file_list)
}

### load files (.tbl) from paths  
load_files_from_paths_for_tbl <- function(file_paths, suf) {
  file_list <- list() # Create an empty list to store the loaded files
  # Iterate over the file paths and read the corresponding files
  for (path in file_paths) {
    file_data <- fread(path)
    filename <- basename(path)
    filename_without_suffix <- sub(suf, "", filename)
    file_list[[filename_without_suffix]] <- file_data
  }
  return(file_list)
}

### Add annotation column to model attributes if not already there
add_annotation_column_to_attributes <- function(mod.orig) {
    if(!("annotation" %in% colnames(mod.orig@mod_attr))) {
    bm_ind <- which(mod.orig@react_id == "bio1")
    annostr <- ""
    if(grepl("Bacteria",mod.orig@react_name[bm_ind]))
        annostr <- "tax_domain:Bacteria"
    if(grepl("Archaea",mod.orig@react_name[bm_ind]))
        annostr <- "tax_domain:Archaea"
    
    mod.orig@mod_attr <- cbind(mod.orig@mod_attr,
                                data.frame(annotation = annostr))
    }
    return(mod.orig)
}

# TRANSFORM DICT into DATA.FRAME
pad_dict_to_dataframe <- function(dict) {
  # Find the maximum length among the lists in the dictionary
  max_length <- max(lengths(dict))
  # Pad the shorter lists with NA to match the maximum length
  padded_dict <- lapply(dict, function(x) {
    if (length(x) < max_length) {
      c(x, rep(NA, max_length - length(x)))
    } else {x}})
  # Convert the padded dictionary to a dataframe
  df <- as.data.frame(do.call(rbind, padded_dict))
  return(df)
}

# STANDARDIZE the duplicated MET_NAME in the RXN
# for every duplicated compound check if any element contains "-c0" (they will be favoured)

standardize_duplicated_met_name <- function(met_id2duplicated_met_name_df, met_id2rxn_id_df, info_all_rxns_mods) {
  for (cpd in rownames(met_id2duplicated_met_name_df)){ 

    # determine final name for the cpd
    matching_elements <- sapply(met_id2duplicated_met_name_df[cpd,], grepl, pattern = "-c0") # Use 'grepl()' to check if any element contains "-c0"
    if (any(matching_elements)) {
      cpd_name <- met_id2duplicated_met_name_df[cpd,which(matching_elements)[1]] # take only the first cmp that contains "-c0"  
    } else {
      cpd_name <- met_id2duplicated_met_name_df[cpd, 1]}
    cat(paste("modifying", cpd_name, "\n"))

    # identify the rxn that include that compound 
    rxnid_with_met_dt <- met_id2rxn_id_df[cpd,!(is.na(met_id2rxn_id_df[cpd,]))]
    
    for (rxn_id in as.list(rxnid_with_met_dt)) {
      for (idx in seq_along(info_all_rxns_mods)) {
        rxn <- info_all_rxns_mods[[idx]]
        if (rxn@react_id == rxn_id) {
          matching_indexes <- which(rxn@met_id == cpd)
          rxn@met_name[matching_indexes] <- as.character(cpd_name)
          info_all_rxns_mods[idx] <- list(rxn)
        }
      }
    }
  }
  return(info_all_rxns_mods)
}

# Custom MEDIAN function for weight normalization
# description: if the number of time a reaction is annotated in the pool of draft models in less than the total number of draft models,
# then modify the median calculation by adding n-times (equal to the number of draft models in the pool minus the number of time that reaction was detected by blast) 
# the maximum weight associated to a rxn (i.e. 100).
custom_median <- function(weight, num.pan, num.mod) {
  values_to_add <- rep(100., num.mod-num.pan[[1]]) # Add n values equal to 100
  weight <- c(weight, values_to_add)
  cm <- median(weight)
  return(cm) 
}

# Custom QUARTILE function for weight normalization
# description: if the number of time a reaction is annotated in the pool of draft models in less than the total number of draft models,
# then modify the weight calculation by adding n-times (equal to the number of draft models in the pool minus the number of time that reaction was detected by blast) 
# the maximum weight associated to a rxn (i.e. 100).
custom_quartile_weight <- function(weight, num.pan, num.mod, th) {
  values_to_add <- rep(100., num.mod-num.pan[[1]]) # Add n values equal to 100
  weight <- c(weight, values_to_add)
  weight <- weight[order(weight)]
  cm <- weight[round(num.mod*th)]
  return(cm) 
}

# Build the pan-Draft model
# input: 
#       1) subSet_rxn_df: list of present reaction ID in the selectes subset for the genomes to use to build the model
#       2) info_all_rxns_mods: list with the information regarding all the reactions
#       3) model_list: list of loaded models  
# return: pan-Draft
build_panDraft <- function(subSet_rxn_df, info_all_rxns_mods, mod_desc) {
  # SUBSET the RXN list
  rxn_subset <- subSet_rxn_df$rxn
  subset_info_rxns_mods <- list()
  # extract reaction information only of the subset
  for (rxn in info_all_rxns_mods) {
    if (rxn@react_id %in% rxn_subset){
      subset_info_rxns_mods <- c(subset_info_rxns_mods, rxn)
    }
  }
  ############################
  # Reconstruct the Pan-Draft model  
  cat("\nConstructing draft model... \n")
  pan.mod <- sybil::modelorg(name = "panDraf_model", id = "panDraf_model")

  pan.mod@mod_compart <- c("c0","e0","p0")
  mod_compart(pan.mod)
  # add metabolite compartment list
  # n.comp <- max(pan.mod@met_comp, na.rm = T)
  # if(n.comp == 2)
  #  pan.mod@mod_compart <- c("c0","e0")
  # if(n.comp == 3)
  #  pan.mod@mod_compart <- c("c0","e0","p0")

  pan.mod@react_attr <- data.frame(seed = character(0), rxn = character(0), name = character(0), ec = character(0), tc = character(0), 
                                qseqid = character(0), pident = numeric(0), evalue = numeric(0), bitscore = numeric(0), qcovs = numeric(0),
                                stitle = character(0), sstart = numeric(0), send = numeric(0), pathway = character(0), status = character(0), 
                                pathway.status = character(0), complex = character(0), exception = numeric(0), complex.status = numeric(0), gs.origin = logical(0), 
                                annotation = character(0), MNX_ID = character(0), seedID = character(0), keggID = character(0), biggID = character(0), 
                                biocycID = character(0), stringsAsFactors = F)

  # generate the subSystems matrix and the react_attr dataframe
  subsys_unique <- list()
  react_attr_list <- list()

  for (rxn in subset_info_rxns_mods){
    subsys_unique <- c(subsys_unique, rxn@subSys)
    subsys_unique <- as.character(unique(subsys_unique))
    react_attr_list <- c(react_attr_list, list(rxn@react_attr))
  }
  pan.mod@subSys <- Matrix::Matrix(F,nrow = 0, ncol = length(subsys_unique),sparse = T)
  colnames(pan.mod@subSys) <- subsys_unique
  react_attr_df <- do.call(rbind, lapply(react_attr_list, data.frame, stringsAsFactors = FALSE)) # Convert list of lists to dataframe

  # add one reaction at the time
  for (i in 1:length(subset_info_rxns_mods)) {
    rxn <- subset_info_rxns_mods[[i]]
    rxn.id <- rxn@react_id
    met.ids <- rxn@met_id
    met.scoef <- rxn@s
    is.rev <- rxn@react_rev
    met.comp <- rxn@met_comp
    rxn.uppbnd <- rxn@uppbnd
    rxn.lowbnd <- rxn@lowbnd
    rxn.name <- rxn@react_name
    met.name <- rxn@met_name
    gpr.tmp <- rxn@gpr
    dts.tmp <- rxn@subSys

    pan.mod <- sybil::addReact(model = pan.mod, 
                            id = rxn.id, 
                            met = met.ids,
                            Scoef = met.scoef,
                            reversible = is.rev, 
                            metComp = met.comp,
                            ub = rxn.uppbnd,
                            lb = rxn.lowbnd,
                            reactName = rxn.name, 
                            metName = met.name,
                            gprAssoc = gpr.tmp,
                            subSys = dts.tmp)

    # define the reactions attribute
    if (grepl("DM", rxn.id)) {
      pan.mod@react_attr[which(pan.mod@react_id == rxn.id),] <- subset(react_attr_df, seed == rxn.id) 
    } else {
      pan.mod@react_attr[which(pan.mod@react_id == rxn.id),] <- subset(react_attr_df, seed == strsplit(rxn.id, "_c0")) 
    }

    # In case no genes are associated to the first reaction added to the model,
    # no 'genes' 'grp' and 'gprRules' object are initiated by sybil, causing a
    # mismatch between length of the number of reactions and those three objects
    # in the final modelorg model. This mismatch causes libSBML to stop on a
    # segmentation error. Solution: "Manually" initiate the three objects in those
    # cases
    if(i == 1 && gpr.tmp == "") {
      pan.mod@genes[[1]] <- ""
      pan.mod@gpr <- ""
      pan.mod@gprRules <- ""
    }
  }

  pan.mod@mod_desc <- mod_desc # mod description of first loaded model
  # add annotation column to model attributes if not already there
  if(!("annotation" %in% colnames(pan.mod@mod_attr))) {
    bm_ind <- which(pan.mod@react_id == "bio1")
    annostr <- ""
    if(grepl("Bacteria",pan.mod@react_name[bm_ind]))
      annostr <- "tax_domain:Bacteria"
    if(grepl("Archaea",pan.mod@react_name[bm_ind]))
      annostr <- "tax_domain:Archaea"
    pan.mod@mod_attr <- data.frame(annotation = annostr)
  }

  # add gapseq version info to model object
  gapseq_version <- system(paste0(script.dir,"/.././gapseq -v"), intern = T)
  pan.mod@mod_desc <- gapseq_version
  
  cat("\tcompleted\n")
  return(pan.mod)
}

# Build the data.table of presence/absence reaction in a list of models
# input: 
#       1) model_list: list of models
# return: rxn2mod_dt
build_rxn2mod_dt <- function(model_list) {
  # Iterate over the rxn in each model to obtain a pre-dataset: list of lists of rxn
  rxn_List <- list()
  mod_List <- list()
  mod_id2mod_dict <- list() # dictionary to save association between mod ID and mod
  
  for (mod in model_list) {
    setattr(mod_id2mod_dict, mod@mod_id, mod) # set dict key-value
    
    # RXN list
    rxn_List <- c(rxn_List, mod@react_id)
    mod_idXreact_num <- c(mod@mod_id, rep(mod@mod_id,length(mod@react_id)-1)) 
    mod_List <- c(mod_List, mod_idXreact_num) 
  }
  
  rxn2mod_List <- list(rxn = rxn_List,
                       mod_id = mod_List)
  
  # RXN: Convert the list to a data.table
  dt <- data.table::copy(rxn2mod_List) 
  setDT(dt)
  
  dt[, presence := 1] # Add a column indicating the presence
  dt[, c("rxn", "mod_id") := lapply(.SD, as.character), .SDcols = c("rxn", "mod_id")] # Convert the columns to character
  rxn2mod_dt <- dcast(dt, rxn ~ mod_id, fun.aggregate = length, fill = 0) # Reshape the data.table
  
  return(list(rxn2mod_dt, mod_id2mod_dict))
}
