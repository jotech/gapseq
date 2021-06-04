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
    xml_lines <- readLines(paste0(out.id, ".xml"))
    # (1) following patch adds an attribute to the sbml-file syntax, that is required by SBML-file validator (http://sbml.org/Facilities/Validator)
    system(paste0("perl -pi -w -e 's/fbc:required=\"false\"/fbc:required=\"false\" groups:required=\"false\"/g;' ", out.id, ".xml"))
    
    # (2) Following patch adds the "groups:id" attribute to each subsystem
    grp_block_start <- grep("<groups:listOfGroups>", xml_lines, fixed = T)
    for(i in 1:ncol(mod@subSys)) {
      grpID <- colnames(mod@subSys)[i]
      grpID_deform <- paste0("subsys_",gsub("-","_", grpID, fixed = T))
      
      grpSearch  <- paste0("groups:name=\"",grpID,"\"")
      grpReplace <- paste0("groups:id=\"",grpID_deform,"\" groups:name=\"",grpID,"\"")
      
      system(paste0("perl -pi -w -e 's/",grpSearch,"/",grpReplace,"/g if $. > ",grp_block_start,"' ", out.id, ".xml"))
    }
    
    # (3) add model name and id
    system(paste0("perl -pi -w -e 's/<model fbc:strict=\"true\">/<model fbc:strict=\"true\" id=\"", mod@mod_id,
                  "\" name=\"",mod@mod_name,"\">/g;' ", out.id, ".xml"))
    
    
  }else{
    print("sybilSBML not found, please install sybilSBML for sbml output")
  }
}