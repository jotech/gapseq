suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))

prefix <- "trpBclust"

tree_file <- paste0(prefix, ".iqtree.treefile")
output_file <- paste0(prefix, ".cluster_assignments.tsv")

tree <- read.tree(tree_file)


# Pairwise branch-length distances between tree tips
tree_distances <- cophenetic.phylo(tree)

# Hierarchical clustering of the patristic-distance matrix
hc <- hclust(
  as.dist(tree_distances),
  method = "average"
)

# Force a two-cluster solution
membership <- cutree(hc, k = 2)
trpb2_idx <- which.min(table(membership))

assignments <- data.frame(
  sequence_id = names(membership),
  cluster = ifelse(membership==trpb2_idx,"TrpB2","TrpB1"),
  stringsAsFactors = FALSE
)
assignments <- data.table(assignments)

ggtree(tree, layout = "daylight") %<+% assignments +
  geom_tippoint(aes(color = cluster), size = 1) +
  scale_color_manual(values = c(
    `TrpB1` = "red",
    `TrpB2` = "black"
  ))

# split sequences into 4.2.1.20.fasta for trpA+trpB1 and 4.2.1.122 for trpB2
allTrpA <- readAAStringSet("tmp_trpA.fasta")
allTrpB <- readAAStringSet("tmp_trpB.fasta")

trpB1seqs <- allTrpB[which(grepl(
  assignments[cluster == "TrpB1", paste0("^",sequence_id, collapse = "|")],
  names(allTrpB)
))]
trpB2seqs <- allTrpB[which(grepl(
  assignments[cluster == "TrpB2", paste0("^",sequence_id, collapse = "|")],
  names(allTrpB)
))]

trpAB <- c(allTrpA,
           trpB1seqs)
trpB2 <- trpB2seqs

writeXStringSet(trpAB, "dat/seq/Bacteria/user/4.2.1.20.fasta")
writeXStringSet(trpAB, "dat/seq/Archaea/user/4.2.1.20.fasta")

writeXStringSet(trpB2, "dat/seq/Bacteria/user/4.2.1.122.fasta")
writeXStringSet(trpB2, "dat/seq/Archaea/user/4.2.1.122.fasta")


# clean up
file.remove(dir(".", pattern = "^trpB|^tmp_trp"))
