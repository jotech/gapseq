library(XML)
library(xml2)
library(data.table)
library(stringr)


source("src/update_ecs_via_IUBMB_functions.R")

#------------------------#
# Get latest EC database #
#------------------------#
ecdb <- download_ecdb()
setkey(ecdb, "ec_num")

# some ECs were transferred to ECs that wer later again transferred to a yet-new
# IDs and so on. This look assigns EC to the final/current valid EC.
lookup_ecs <- ecdb[!is.na(ec_num_new), ec_num_new]
names(lookup_ecs) <- ecdb[!is.na(ec_num_new), ec_num]
max_rounds <- 10
i <- 1
for(i in 1:max_rounds) {
  if(ecdb[ec_num_new %in% names(lookup_ecs), .N] > 0) {
    ecdb[ec_num_new %in% names(lookup_ecs), ec_num_new := lookup_ecs[ec_num_new]]
  } else {
    break
  }
}
if(i == max_rounds)
  warning("Potential loop in EC transfers (e.g., ECA -> ECB -> ECA). This should not happen.")
rm(i)

#--------------------------#
# Update ECs in data files #
#--------------------------#

# Pathways (dat/meta_pwy.tbl and dat/custom_pwy.tbl)
updpwy <- correctEC_pathways(ecdb)

# ec conflicts (dat/ec_conflicts.tsv)
updecc <- correctEC_conflicts(ecdb)

# ec to gaspeq/ModelSEED reactions
updrxn <- correctEC_seedrxn(ecdb)

# ec exceptions?
updexc <- correctEC_exceptions(ecdb)

