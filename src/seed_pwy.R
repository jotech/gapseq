library(data.table)
library(stringr)

seedDB <- fread("~/uni/dat/db/seed/ModelSEED_Subsystems.tsv")

length()

seedPwy <- unique(seedDB$Name)
seedPwy.dt <- data.table()
for(i in seq_along(seedPwy)){
  pwy           <- seedPwy[i]
  pwy.id        <- paste0("seed.pwy",i)
  pwy.dat       <- seedDB[Name==pwy]
  rea.name.full <- unique(pwy.dat$Role)
  rea.id        <- pwy.dat$Reaction[match(rea.name.full, pwy.dat$Role)] # same reaction name can have several IDs (take only first here)
  rea.ec        <- str_extract(rea.name.full, "(?<=\\(EC ).*?(?=\\))")
  rea.noec      <- length(na.omit(rea.ec))
  rea.ec        <- ifelse(is.na(rea.ec),"", rea.ec)
  rea.name      <- str_remove(rea.name.full, " \\(EC.*\\)")
  hierarchy     <- paste0(pwy.dat$Class[1], ";",pwy.dat$`Sub-class`[1])
  seedPwy.dt    <- rbind(seedPwy.dt,
                         data.table(id=pwy.id, name=pwy, altname="", hierarchy=hierarchy, taxrange="", reaId=paste(rea.id, collapse = ","), 
                         reaEc=paste(rea.ec, collapse = ","), keyRea="", reaName=paste(rea.name, collapse = ";"), 
                         reaNr=length(rea.id), ecNr= rea.noec, superpathway="False", status="True"))
}
seedPwy.dt[reaNr!=ecNr]

seedPwy.dt$name <- gsub("_"," ", seedPwy.dt$name)
seedPwy.dt$hierarchy <- gsub("_"," ", seedPwy.dt$hierarchy)

fwrite(seedPwy.dt, file = "~/uni/gapseq/dat/seed_pwy.tbl", sep="\t", quote = F)
