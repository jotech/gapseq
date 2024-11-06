library(httr)
library(stringr)
library(parallel)

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
batch_size_clusters <- 100 # 200 should work: https://github.com/ebi-uniprot/uniprot-rest-api/issues/275#issuecomment-1173888616
max_attempts <- 10 # Maximum number of attempts per query to receive a status==2xx response from uniprot website.

GET_retries <- function(url) {
  require(httr)
  attempts <- 1
  get_success <- FALSE
  res <- NULL

  while(!get_success && attempts <= max_attempts) {
    tryCatch({
      res <- GET(url)
      get_success <- res$status_code >= 200 && res$status_code < 300
    }, error = function(e) {
      message(paste("Attempt", attempts, "failed with error:", conditionMessage(e)))
      # Error handling: reset res and continue to retry
      res <- NULL
    })

    if (!get_success) {
      attempts <- attempts + 1
      # cat(paste0(url,"\n   ATTEMPT: ",attempts-1,"; status: ", res$status_code,"\n\n"), file = "/home/silvio/Software/gapseq/test.log", append = TRUE)
      Sys.sleep(1)  # Pause a sec between attempts.
    }

  }
  if(!get_success) {
    cat(NULL, file = output_fasta_file)  # Create empty file
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
query_term <- utils::URLencode(args[2], reserved = TRUE)

if(query_type == "protein_name") {
  query_term <- paste0("%22",query_term,"%22")
}

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
taxonomy <- paste0("%28taxonomy_name:",taxonomy,"%29", collapse = "%20OR%20")
taxonomy <- paste0("%28",taxonomy,"%29")

if(query_type != "id") {
  # (0) - Build uniprot entry URL
  # urli <- "https://rest.uniprot.org/uniprotkb/search?fields=accession&format=list&query=%28%28ec%3A1.4.4.2%29+AND+%28%28taxonomy_name:Bacteria%29%29+AND+%28reviewed%3Afalse%29%29&size=500"
  urli <- paste0("https://rest.uniprot.org/uniprotkb/search?fields=accession&format=list&query=%28%28",
                 query_type,"%3A",query_term,"%29+AND+",taxonomy,"+AND+%28reviewed%3A",rev_status,"%29%29&size=500")

  # (1) - Find uniprot accession numbers of matching entries
  # Note: This uses pagination
  accessions <- character(0L)
  while(!is.null(urli)) {
    ri <- GET(urli)
    acctmp <- content(ri, as = "text", encoding = "UTF-8")
    acctmp <- unlist(str_split(acctmp, "\n"))
    acctmp <- acctmp[acctmp != ""]
    accessions <- c(accessions, acctmp)

    linkurl <- ri$headers$link
    if(is.null(linkurl)) {
      urli <- NULL
    } else {
      urli <- gsub("^<|>.*$","",linkurl)
    }
    cat("\r",length(accessions))
  }
  total <- length(accessions)
  if(length(accessions) == 0) {
    cat(NULL, file = output_fasta_file)
    quit(save = "no", status = 0)
  }
#   ri <- GET_retries(urli)
#   accessions <- content(ri, as = "text", encoding = "UTF-8")
#   accessions <- unlist(str_split(accessions, "\n"))[-1]
#   accessions <- accessions[accessions != ""]


  # (2) retrieve corresponding uniref90 cluster IDs (not sequences)
  acc_batches <- split(accessions, ceiling(seq_along(accessions)/batch_size_clusters))
  cl <- makeCluster(min(length(acc_batches),15))
  clusterExport(cl, c("GET_retries","uniref_id","max_attempts","output_fasta_file"))
  map_worker <- function(accvec) {
    require(stringr)
    acc_concat <- paste0("%28uniprot_id%3A",accvec,"%29",
                         collapse = "%20OR%20")
    urlc <- paste0("https://rest.uniprot.org/uniref/search?compressed=false&fields=id&format=tsv&query=%28%28",
                   acc_concat,"%29%20AND%20%28identity%3A",uniref_id,"%29%29&size=500")
    ri <- GET_retries(urlc)
    uniref_acc <- content(ri, as = "text", encoding = "UTF-8")
    uniref_acc <- unlist(str_split(uniref_acc, "\n"))[-1]
    uniref_acc <- uniref_acc[uniref_acc != ""]
  }
  cluster_ids <- parLapply(cl, acc_batches, map_worker)


  cluster_ids <- unique(unlist(cluster_ids))
  if(length(cluster_ids)==0) {
    cat(NULL, file = output_fasta_file)
    quit(save = "no", status = 1)
  }

  cat(" (",length(cluster_ids),"cluster sequences found )\n")
  total_uniref <- length(cluster_ids)

  # (3) retrieve Uniref sequences
  uniref_batches <- split(cluster_ids, ceiling(seq_along(cluster_ids)/batch_size_clusters))

  clust_worker <- function(clustvec) {
    require(stringr)
    acc_concat <- paste0("%28id%3A",clustvec,"%29",
                         collapse = "%20OR%20")
    urlc <- paste0("https://rest.uniprot.org/uniref/search?format=fasta&query=%28",
                   acc_concat,"%29&size=500")
    ri <- GET_retries(urlc)
    seqs <- content(ri, as = "text", encoding = "UTF-8")
  }
  cluster_seqs <- parLapply(cl, uniref_batches, clust_worker)
  stopCluster(cl)

  cluster_seqs <- unlist(cluster_seqs)
  if(length(cluster_seqs)==0) {
    cat(NULL, file = output_fasta_file)
    quit(save = "no", status = 1)
  }
  cat(cluster_seqs, sep = "", file = output_fasta_file)
  quit(save = "no", status = 0)
}

if(query_type == "id") {
  query_term <- unique(unlist(str_split(query_term,";")))

  seqi <- c()
  for(qti in query_term) {
    # (0) - Build uniprot entry URL
    urli <- paste0("https://rest.uniprot.org/uniprotkb/search?compressed=false&format=fasta&query=%28%28accession%3A",
                   qti,"%29%20AND%20",taxonomy,"%29")

    ri <- GET_retries(urli)

    seqi <- c(seqi, content(ri, as = "text", encoding = "UTF-8"))
  }

  if(any(grepl("^Error", seqi))) {
    cat(NULL, file = output_fasta_file)
    quit(save = "no", status = 1)
  }

  cat(seqi, sep = "", file = output_fasta_file)
  quit(save = "no", status = 0)
}
