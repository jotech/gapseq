write_gapseq_sbml <- function(mod, out.id) {
  
  if( "sybilSBML" %in% rownames(installed.packages()) ){
    
    # handling of missing values/attributes
    if( any(is.na(mod@met_attr$charge)) ) mod@met_attr$charge[which(is.na(mod@met_attr$charge))] <- ""
    if( any(is.na(mod@met_attr$chemicalFormula)) ) mod@met_attr$chemicalFormula[which(is.na(mod@met_attr$chemicalFormula))] <- ""
    if( any( mod@met_attr$chemicalFormula=="null"))mod@met_attr$chemicalFormula[which(mod@met_attr$chemicalFormula=="null")]<- ""
    
   
    # handling of subunit names (and remove empty subunits)
    colnames(mod@subSys) <- gsub("^\\||\\|$","",colnames(mod@subSys))
    mod@subSys <- mod@subSys[, apply(mod@subSys,2,any)]
    
    # gpr terms
    mod@gpr <- gsub("\\&", "and", mod@gpr)
    mod@gpr <- gsub("\\|", "or",  mod@gpr)

    # Model notes to include gapseq version
    notes_str <- paste0("<notes>\n  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>",mod@mod_desc,"</p>\n  </html>\n</notes>")
    mod@mod_attr <- data.frame(notes = notes_str)
    
    # writing sbml
    sbml.o <- sybilSBML::writeSBML(mod, filename = paste0(out.id, ".xml"), 
                                   level = 3, version = 1, fbcLevel = 2, 
                                   printNotes = T, printAnnos = T)
    if(sbml.o==F)
      warning("Writing SBML-file failed.")
    
    
    # - - - - - - - - #
    # Patch xml file  #
    # - - - - - - - - #
    cat("Patching file ",out.id,".xml ...", sep = "")
    xml_lines <- readLines(paste0(out.id, ".xml"))
    n_lines <- length(xml_lines)
    # (1) following patch adds an attribute to the sbml-file syntax, that is required by SBML-file validator (http://sbml.org/Facilities/Validator)
    indtmp <- grep("fbc:required=\"false\">", xml_lines, fixed = TRUE)
    xml_lines[indtmp] <- gsub("fbc:required=\"false\">","fbc:required=\"false\" groups:required=\"false\">",
                              xml_lines[indtmp], fixed = TRUE)
    
    # (2) Following patch adds the "groups:id" attribute to each subsystem
    grp_block_start <- grep("<groups:listOfGroups>", xml_lines, fixed = TRUE)
    for(i in 1:ncol(mod@subSys)) {
      grpID <- colnames(mod@subSys)[i]
      grpID_deform <- paste0("subsys_",gsub("-","_", grpID, fixed = TRUE))
      
      grpSearch  <- paste0("groups:name=\"",grpID,"\"")
      grpReplace <- paste0("groups:id=\"",grpID_deform,"\" groups:name=\"",grpID,"\"")
      
      indtmp <- grep(grpSearch, xml_lines[(grp_block_start+1):n_lines])
      xml_lines[grp_block_start + indtmp] <- gsub(grpSearch, grpReplace,
                                                  xml_lines[grp_block_start + indtmp],
                                                  fixed = TRUE) 
    }
    
    # (3) add model name and id
    indtmp <- grep("<model fbc:strict=\"true\">", xml_lines, fixed = TRUE)
    xml_lines[indtmp] <- gsub("<model fbc:strict=\"true\">",
                              paste0("<model fbc:strict=\"true\" id=\"",
                                     gsub("\\.","_",mod@mod_id),
                                     "\" name=\"",mod@mod_name,"\">"),
                              xml_lines[indtmp], fixed = TRUE)
    
    # (4) Gene name correction
    ind_gp <- grep("fbc:geneProduct", xml_lines, fixed = T)
    
    xml_lines[ind_gp] <- gsub("__MINUS__","_", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__COLON__","_", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__DOT__","_", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__ONE__","1", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__TWO__","2", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__THREE__","3", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__FOUR__","4", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__FIVE__","5", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__SIX__","6", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__SEVEN__","7", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__EIGHT__","8", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__NINE__","9", xml_lines[ind_gp])
    xml_lines[ind_gp] <- gsub("__ZERO__","0", xml_lines[ind_gp])
    
    # rewrite patched sbml file
    writeLines(xml_lines, paste0(out.id, ".xml"))
    
    cat(" done\n")
    
  }else{
    print("R-package sybilSBML not found. Please install sybilSBML for sbml output.")
  }
}
