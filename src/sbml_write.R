write_gapseq_sbml <- function(mod, out.id) {
  
  if( "sybilSBML" %in% rownames(installed.packages()) ){
    
    # handling of missing values/attributes
    if( any(is.na(mod@met_attr$charge)) ) mod@met_attr$charge[which(is.na(mod@met_attr$charge))] <- ""
    if( any(is.na(mod@met_attr$chemicalFormula)) ) mod@met_attr$chemicalFormula[which(is.na(mod@met_attr$chemicalFormula))] <- ""
    if( any( mod@met_attr$chemicalFormula=="null"))mod@met_attr$chemicalFormula[which(mod@met_attr$chemicalFormula=="null")]<- ""
    
    # handling of subunit names (and remove empty subunits)
    colnames(mod@subSys) <- gsub("^\\||\\|$","",colnames(mod@subSys))
    mod@subSys <- mod@subSys[, apply(mod@subSys,2,any)]
    
    # writing sbml
    sbml.o <- sybilSBML::writeSBML(mod, filename = paste0(out.id, ".xml"), 
                                   level = 3, version = 1, fbcLevel = 2, 
                                   printNotes = T, printAnnos = T)
    if(sbml.o==F)
      warning("Writing SBML-file failed.")
    
    # - - - - - - - - #
    # Patch xml file  #
    # - - - - - - - - #
    
    # following patch adds an attribute to the sbml-file syntax, that is required by SBML-file validator (http://sbml.org/Facilities/Validator)
    system(paste0("perl -pi -w -e 's/fbc:required=\"false\"/fbc:required=\"false\" groups:required=\"false\"/g;' ", out.id, ".xml"))
    
    # Following patch adds the "groups:id" attribute to each subsystem
    for(i in 1:ncol(mod@subSys)) {
      grpID <- colnames(mod@subSys)[i]
      grpID_deform <- paste0("subsys_",gsub("-","_", grpID, fixed = T))
      
      grpSearch  <- paste0("groups:name=\"",grpID,"\"")
      grpReplace <- paste0("groups:id=\"",grpID_deform,"\" groups:name=\"",grpID,"\"")
      
      system(paste0("perl -pi -w -e 's/",grpSearch,"/",grpReplace,"/g;' ", out.id, ".xml"))
    }
    
  }else{
    print("sybilSBML not found, please install sybilSBML for sbml output")
  }
}