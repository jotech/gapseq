# Tools required
# * clipkit (intall via pip)
# * iqtree3 (apt)
# * seqkit (apt)
# * mafft (apt)


library(httr)
library(stringr)
suppressMessages(library(Biostrings))

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

get_uniref90IDs <- function(accvec) {
  cat("\r",k,"/",length(acc_batches))
  acc_concat <- paste0("%28uniprot_id%3A",accvec,"%29",
                       collapse = "%20OR%20")
  urlc <- paste0("https://rest.uniprot.org/uniref/search?compressed=false&fields=id&format=tsv&query=%28%28",
                 acc_concat,"%29%20AND%20%28identity%3A0.9%29%29&size=500")
  ri <- GET_retries(urlc)
  uniref_acc <- content(ri, as = "text", encoding = "UTF-8")
  uniref_acc <- unlist(str_split(uniref_acc, "\n"))[-1]
  uniref_acc <- uniref_acc[uniref_acc != ""]
  k <<- k + 1
  return(uniref_acc)
}

get_uniref_seqs <- function(clustvec) {
  require(stringr)
  acc_concat <- paste0("%28id%3A",clustvec,"%29",
                       collapse = "%20OR%20")
  urlc <- paste0("https://rest.uniprot.org/uniref/search?format=fasta&query=%28",
                 acc_concat,"%29&size=500")
  ri <- GET_retries(urlc)
  seqs <- content(ri, as = "text", encoding = "UTF-8")
}

# Uniprot query settings
batch_size_clusters <- 100 # 200 should work: https://github.com/ebi-uniprot/uniprot-rest-api/issues/275#issuecomment-1173888616
max_attempts <- 10 # Maximum number of attempts per query to receive a status==2xx response from uniprot website.

# TrpA
# Query: (ec:4.2.1.20) AND ((taxonomy_id:2) OR (taxonomy_id:2157)) AND (reviewed:true) AND ((gene:trpA) OR (gene:trpA1) OR (gene:trpA2))

query_url_trpA <- "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28ec%3A4.2.1.20%29+AND+%28%28taxonomy_id%3A2%29+OR+%28taxonomy_id%3A2157%29%29+AND+%28reviewed%3Atrue%29+AND+%28%28gene%3AtrpA%29+OR+%28gene%3AtrpA1%29+OR+%28gene%3AtrpA2%29%29%29"

ri <- GET(query_url_trpA)
trpAacc <- content(ri, as = "text", encoding = "UTF-8")
trpAacc <- unlist(str_split(trpAacc, "\n"))
trpAacc <- trpAacc[trpAacc != ""]

acc_batches <- split(trpAacc, ceiling(seq_along(trpAacc)/batch_size_clusters))
k <- 1
cluster_ids <- lapply(acc_batches, FUN = get_uniref90IDs)
cluster_ids <- unique(unlist(cluster_ids))

uniref_batches <- split(cluster_ids, ceiling(seq_along(cluster_ids)/batch_size_clusters))
cluster_seqs <- lapply(uniref_batches, get_uniref_seqs)
cluster_seqs <- unlist(cluster_seqs)
cat(cluster_seqs, sep = "", file = "tmp_trpA.fasta")

# TrpB
# Query: (ec:4.2.1.20) AND ((taxonomy_id:2) OR (taxonomy_id:2157)) AND (reviewed:true) AND ((gene:trpB) OR (gene:trpB1) OR (gene:trpB2))

query_url_trpB <- "https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28ec%3A4.2.1.20%29+AND+%28%28taxonomy_id%3A2%29+OR+%28taxonomy_id%3A2157%29%29+AND+%28reviewed%3Atrue%29+AND+%28%28gene%3AtrpB%29+OR+%28gene%3AtrpB1%29+OR+%28gene%3AtrpB2%29%29%29"

ri <- GET(query_url_trpB)
trpBacc <- content(ri, as = "text", encoding = "UTF-8")
trpBacc <- unlist(str_split(trpBacc, "\n"))
trpBacc <- trpBacc[trpBacc != ""]

acc_batches <- split(trpBacc, ceiling(seq_along(trpBacc)/batch_size_clusters))
k <- 1
cluster_ids <- lapply(acc_batches, FUN = get_uniref90IDs)
cluster_ids <- unique(unlist(cluster_ids))

uniref_batches <- split(cluster_ids, ceiling(seq_along(cluster_ids)/batch_size_clusters))
cluster_seqs <- lapply(uniref_batches, get_uniref_seqs)
cluster_seqs <- unlist(cluster_seqs)
cat(cluster_seqs, sep = "", file = "tmp_trpB.fasta")

#-------------------------------------------------------------------------------
# Reconstruct pyhlogenetic tree for trpB sequences (to identify trpB1 and trpB2
# clusters)
#-------------------------------------------------------------------------------
system("src/./trpB_phylo.sh tmp_trpB.fasta trpBclust 16", intern = TRUE)





