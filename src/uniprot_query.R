library(httr)
library(stringr)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#
# Programmatic download of protein cluster sequences from EC / Protein name /  #
# Gene ID queries via Uniprot REST requires 3 steps:                           #
#  (1) Find UniprotKB ID/AC of query-matching proteins                         #
#  (2) Match IDs/Accessions from (1) to identity-based protein clusters        #
#      (UniRef_90/UniRef_50)                                                   #
#  (3) Retrieve representative Sequences of protein clusters                   #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~#

# query_type <- "ec" # (either "ec", "protein_name", "gene_exact" or "id").
# query_term <- "2.7.1.90"
batch_size_query_results <- 500
batch_size_clusters <- 200 # 200 should work: https://github.com/ebi-uniprot/uniprot-rest-api/issues/275#issuecomment-1173888616
max_attempts <- 10 # Maximum number of attempts per query to receive a status==200 response from uniprot website.

GET_retries <- function(url) {
  attempts <- 1
  get_success <- FALSE
  while(get_success == FALSE && attempts <= max_attempts) {
    res <- GET(url)
    attempts <- attempts + 1
    get_success <- res$status_code
  }
  if(get_success == FALSE) {
    quit(save = "no", status = 1)
  }
  
  return(res)
}
# rev_status <- "false"
# uniref_id <- "0.5"
# output_fasta_file <- "test.fasta"
# taxonomy <- "Bacteria"

# positional arguments:
# 1: query type (EC, protein_name, gene_exact)
# 2: query term
# 3: output file name
# 4: taxonomy
# 5: reviewed status (true or false)
# 6: UniRef identity clustering cutoff (0.9 or 0.5)

args = commandArgs(trailingOnly=TRUE)

# parse positional arguments
query_type <- args[1]
query_term <- utils::URLencode(args[2])
output_fasta_file <- args[3]
taxonomy <- args[4]
if(length(args)==6) {
  rev_status <- args[5]
  uniref_id <- args[6]
}

# check if all required arguments as provided
if (query_type!="id" & length(args)!=6) {
  stop("")
}
if (query_type=="id" & length(args)!=4) {
  stop("")
}

# parse taxonomy term
taxonomy <- unlist(strsplit(taxonomy, ","))
taxonomy <- paste0("(taxonomy_name:",taxonomy,")", collapse = "%20OR%20")
taxonomy <- paste0("(",taxonomy,")")

if(query_type != "id") {
  # (0) - Build uniprot entry URL 
  urli <- paste0("https://rest.uniprot.org/uniprotkb/search?format=list&query=(",
                 query_type,":",query_term,
                 ")%20AND%20",taxonomy,"%20AND%20(reviewed:",
                 rev_status,")&size=",
                 batch_size_query_results)
  
  # (1) - Find uniprot accession numbers of matching entries
  accessions <- c()
  total <- NA_integer_
  # cat("Retrieving sequence accession IDs ...\n")
  while(length(urli)==1) {
    ri <- GET_retries(urli)
    ri_header <- headers(ri)
    
    # get metrics and accessions
    total <- as.numeric(ri_header$`x-total-results`)
    accessions <- c(accessions, unlist(strsplit(httr::content(ri, as = "text", encoding = "UTF-8"),"\n")))
    cat("\r\t\t\t",length(accessions), "/", total)
    
    ind_next <- grep("\"next\"$",ri_header$link)
    
    next_url <-  str_match(ri_header$link[ind_next], "<\\s*(.*?)\\s*>")[,2]
    
    urli <- next_url
  }
  cat("\n")
  if(length(accessions) == 0) {
    cat(NULL, file = output_fasta_file)
    quit(save = "no", status = 0)
  }
  
  
  # (2) retrieve corresponding uniref90 cluster IDs (not sequences)
  batch_ids <- floor(1:length(accessions)/batch_size_clusters)
  cluster_ids <- c()
  k <- 0
  
  for(batchi in unique(batch_ids)) {
    k <- k + sum(batch_ids == batchi)
    cat("\r\t\t\t",k,"/",total)
    acc_concat <- paste0("%28uniprot_id%3A",accessions[batch_ids == batchi],"%29",
                         collapse = "%20OR%20")
    urlc <- paste0("https://rest.uniprot.org/uniref/search?format=list&query=%28%28",
                   acc_concat,"%29%20AND%20%28identity%3A",uniref_id,"%29%29",
                   "&size=", batch_size_query_results)
    
    ri <- GET_retries(urlc)
    cluster_ids <- c(cluster_ids,
                     unlist(strsplit(httr::content(ri, as = "text",
                                                   encoding = "UTF-8"),"\n")))
  }
  cluster_ids <- unique(cluster_ids)
  if(length(cluster_ids)==0)
    quit(save = "no", status = 1)
  cat(" (",length(cluster_ids),"cluster sequences found)\n")
  total_uniref <- length(cluster_ids)
  
  # (3) retrieve Uniref sequences
  batch_ids <- floor(1:length(cluster_ids)/batch_size_clusters)
  seqs <- c()
  k <- 0
  
  for(batchi in unique(batch_ids)) {
    k <- k + sum(batch_ids == batchi)
    cat("\r\t\t\t",k,"/",total_uniref)
    uniref_acc_concat <- paste0("id%3A",cluster_ids[batch_ids == batchi],
                                collapse = "%20OR%20")
    
    urls <- paste0("https://rest.uniprot.org/uniref/search?compressed=false&format=fasta&query=",
                   uniref_acc_concat,"&size=", batch_size_query_results)
    
    ri <- GET_retries(urls)
    
    seqs <- c(seqs, httr::content(ri, as = "text", encoding = "UTF-8"))
  }
  cat("\n")
  if(length(seqs)==0)
    quit(save = "no", status = 1)
  
  cat(seqs, sep = "", file = output_fasta_file)
  quit(save = "no", status = 0)
}


if(query_type == "id") {
  # (0) - Build uniprot entry URL 
  urli <- paste0("https://rest.uniprot.org/uniprotkb/search?compressed=false&format=fasta&query=%28%28accession%3A",
                 query_term,"%29%20AND%20",taxonomy,"%29")
  
  ri <- GET_retries(urli)
  
  seqi <- httr::content(ri, as = "text", encoding = "UTF-8")
  cat(seqi, sep = "", file = output_fasta_file)
  quit(save = "no", status = 0)
}

