get_gene_logic_string <- function(complex.str, gene) {
  require(data.table)
  require(stringr)
  if(length(complex.str) != length(gene))
    stop("ERROR: length of complex vector unequal to length of gene vector")
  if(any(duplicated(gene)))
    stop("ERROR: Duplicates in gene vector.")
  
  dt.gl <- data.table(complex = complex.str, gene = gene)
  
  su.n <- length(unique(dt.gl$complex))
  
  #dt.gl.cplx <- dt.gl[!is.na(complex) & complex != "Subunit undefined"]
  dt.gl.cplx <- dt.gl[!is.na(complex)]
  
  if(nrow(dt.gl.cplx) > 0) {
    strings_cx <- dt.gl.cplx[, paste(gene, collapse = " | "), by = complex]
    strings_cx[, V1 := paste0("(",V1,")")]
    
    str_cx <- paste(strings_cx$V1, collapse = " & ")
    str_cx <- paste0("(",str_cx,")")
  } else {
    str_cx <- NA
  }
  
  #dt.gl.mono <- dt.gl[is.na(complex) | complex == "Subunit undefined"]
  dt.gl.mono <- dt.gl[is.na(complex)]
  
  if(nrow(dt.gl.mono) > 0) {
    str_mo <- paste(dt.gl.mono$gene, collapse = " | ")
  } else {
    str_mo <- NA
  }
  
  if(!is.na(str_cx) & !is.na(str_mo))
    return(paste(str_cx,"|",str_mo))
  if(is.na(str_cx) & !is.na(str_mo))
    return(str_mo)
  if(!is.na(str_cx) & is.na(str_mo))
    return(str_cx)
  if(is.na(str_cx) & is.na(str_mo))
    return("")
}