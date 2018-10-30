library(getopt)
library(stringr)
library(stringi)

# test data
#blast.res <- "../sw.locale/GCA_000020425.1-all-blast.tbl"
#gram <- "auto"
#genome.seq <- "../sw.locale/GCA_000020425.1.fna.gz"
#high.evi.rxn.BS <- 200
#topology.evidence <- "../sw.locale/GCA_000020425.1-all-Reactions.lst,../sw.locale/GCA_000020425.1-Transporter.lst"

# Please note: If Game-Staining is set to "auto" it requires the external programms barrnap, usearch and bedtools
build_draft_model_from_blast_results <- function(blast.res, topo.evi = NA, gram = "auto", model.name = NA, genome.seq = NA, script.dir, high.evi.rxn.BS = 200) {
  require(data.table)
  require(stringr)
  require(sybil)
  
  if(is.na(model.name))
    model.name <- gsub("-all-blast.tbl","",basename(blast.res), fixed = T)
  
  if(is.na(topo.evi[1]))
    topo.evi <- c()
  
  if(gram=="auto") {
    if(grepl("\\.gz$", genome.seq)) {
      require(R.utils)
      genome.seq <- gunzip(genome.seq, remove = F, temporary = T, overwrite = T)
    } else {
      
    }
    system(paste0("barrnap --quiet ",genome.seq, " | grep 16S > ", genome.seq, ".gff"))
    system(paste0("bedtools getfasta -fi ",genome.seq," -bed ", genome.seq,".gff > ", genome.seq, ".16S.fasta"))
    system(paste0("usearch -sinaps ",genome.seq,".16S.fasta -db " ,script.dir,"/../dat/seq/Bacteria/16S_graminfo/16S_gramposinfo.fna -attr grampos -tabbedout ",genome.seq,".graminfo.tsv -strand plus"))
    
    gram.dt <- fread(paste0(genome.seq,".graminfo.tsv"), header = F)
    gram.dt[, max.score := max(V3)]
    gram.dt <- gram.dt[V3 > max.score * .95]
    gram <- names(sort(table(gram.dt$V2),decreasing=TRUE)[1])
    cat("\nPredicted gram staining: ",gram,"\n")
  }
  
  seed_x_ec <- fread("../dat/seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv", header=T)
  seed_x_name <- fread("../dat/seed_reactions_corrected.tsv", header=T, stringsAsFactors = F)
  seed_x_metCyc <- fread("../dat/mnxref_seed-other.tsv", header = T)
  seed_x_aliases <- fread("../dat/seed_Enzyme_Name_Reactions_Aliases.tsv", header=T)
  
  
  dt <- fread(blast.res, header=T)
  dt <- dt[,.(rxn, name, ec, qseqid, pident, evalue, bitscore, qcovs, stitle, sstart, send)]
  dt[str_count(ec,"\\.")==1, ec := paste0(ec,".-.-")]
  dt[str_count(ec,"\\.")==2, ec := paste0(ec,".-")]
  
  
  # 1. assign modelSEED ID via ec-number
  dt_1 <- merge(dt, seed_x_ec[,.(`External ID`,seed=`MS ID`)], by.x = "ec", by.y = "External ID")
  
  # 2. assign modelSEED ID via mnxref namespace
  dt_2 <- merge(dt, seed_x_metCyc[,.(seed, other)] , by.x = "rxn", by.y = "other")
  
  # for remaining unmapped MetaCyc reaction try this:
  # 3. assign modelSEED ID via reaction name
  matches <- unique(c(dt_1[,paste(rxn, name, sep = "$")], dt_2[,paste(rxn, name, sep = "$")]))
  dt_3 <- merge(dt[!(paste(rxn,name,sep="$") %in% matches)], seed_x_name[,.(seed=id,name)] , by.x = "name", by.y = "name")
  matches <- c(matches, dt_3[,paste(rxn, name, sep = "$")])
  dt_4 <- merge(dt[!(paste(rxn,name,sep="$") %in% matches)], seed_x_aliases[,.(seed=seed.rxn.id,enzyme.name)], by.x = "name",by.y = "enzyme.name")
  matches <- c(matches, dt_4[,paste(rxn, name, sep = "$")])
  
  dt_seed <- rbindlist(list(dt_1,dt_2,dt_3,dt_4))
  dt_seed <- splitcol2rows_mget(dt_seed, "seed", "|")
  
  #TODO: for developing purposes (will be removed later on.) Write file with unmapped meta-cyc reactions (to seed DB)
  dt_unmap <- dt[!(paste(rxn,name,sep="$") %in% matches), .(rxn, name,ec, bitscore)]
  fwrite(dt_unmap, file=paste0(blast.res,"_unmapps_MCrxns.csv"), quote = FALSE, sep="\t")
  
  mseed <- seed_x_name
  mseed <- mseed[gapseq.status %in% c("approved","corrected")]
  mseed <- mseed[order(id)]
  
  rxns <- unique(dt_seed[bitscore >= high.evi.rxn.BS,seed])
  
  mseed <- mseed[(id %in% rxns) | (id %in% topo.evi)]
  
  #remove duplicate reactions
  mseed[, rxn.hash := generate_rxn_stoich_hash(stoichiometry, reversibility)]
  mseed <- mseed[order(rxn.hash,id)]
  
  mod <- modelorg(name = model.name,id = model.name)
  for(i in (1:nrow(mseed))) {
    cat("\r",i,"/",nrow(mseed))
    rxn.info <- str_split(unlist(str_split(string = mseed[i,stoichiometry],pattern = ";")), pattern = ":", simplify = T)
    
    met.comp  <- rxn.info[,3]
    met.comp.n <- ifelse(met.comp==0,"c0","e0")
    met.comp.n <- ifelse(met.comp>=2,"p0",met.comp.n)
    
    met.ids   <- paste0(rxn.info[,2],"[",met.comp.n,"]")
    met.scoef <- as.numeric(rxn.info[,1])
    met.name  <- str_replace_all(rxn.info[,5],"\\\"","")
    
    met.name  <- paste0(met.name,"-",met.comp.n)
    
    ind <- which(met.name=="" | is.na(met.name))
    met.name[ind] <- met.ids[ind]
    
    is.rev <- ifelse(mseed[i,reversibility] %in% c("<","="),T,F)
    only.backwards <- ifelse(mseed[i,reversibility]=="<",T,F)
    
    ind.new.mets <- which(met.ids %in% mod@met_id)
    ind.old.mets <- which(mod@met_id %in% met.ids[ind.new.mets])
    
    met.name[ind.new.mets] <- mod@met_name[ind.old.mets]
    
    mod <- addReact(model = mod, 
                    id = paste0(mseed[i,id],"_c0"), 
                    met = met.ids,
                    Scoef = met.scoef,
                    reversible = is.rev, 
                    metComp = as.integer(met.comp)+1,
                    ub = ifelse(only.backwards, 0, 1000),
                    lb = ifelse(is.rev, -1000, 0),
                    reactName = mseed[i, name], 
                    metName = met.name)
  }
  cat("\n")
  
  # Adding Biomass reaction
  if(gram == "neg")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramNeg.tsv"))
  if(gram == "pos")
    dt.bm <- fread(paste0(script.dir, "/../dat/seed_biomass.DT_gramPos.tsv"))
  
  mod <- addReact(mod,id = "bio1", met = dt.bm$id, Scoef = dt.bm$stoich, reversible = F, lb = 0, ub = 1000, obj = 1, 
                  reactName = paste0("Biomass reaction ",ifelse(gram=="neg","(gram -)","(gram +)")), metName = dt.bm$name, metComp = dt.bm$comp)
  
  mod <- add_missing_exchanges(mod) 
  
  # get human readlible model/bacterium name
  assembly.names <- fread("/nuuk/2018/bacref/metadata/assemblyTaxo.tsv", quote = "", fill = T)
  mod@mod_desc <- assembly.names[AssemblyID == model.name, ScientificName]
  
  # construct gapfill candidate reactions and their weight
  dt.cand <- dt_seed[bitscore < high.evi.rxn.BS]
  dt.cand[, max.bs := max(bitscore), by = "seed"]
  dt.cand <- dt.cand[max.bs == bitscore]
  dt.cand <- dt.cand[!duplicated(seed)]
  dt.cand[, max.bs := NULL]
  dt.cand[, weight := 1 - (bitscore / high.evi.rxn.BS)]
  
  return(list(mod=mod, cand.rxns=dt.cand))
}

splitcol2rows_mget <- function(dtInput, col2split, sep){
  dtInput <- dtInput[, .(tmp.add.col = unlist(strsplit(get(col2split),sep,T))), by=names(dtInput)]
  
  dtInput[, c(col2split):=NULL];
  setnames(dtInput, 'tmp.add.col', col2split); 
  return(dtInput);
}

# get current script path
if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
  # RStudio specific code
  script.dir    <- dirname(rstudioapi::getSourceEditorContext()$path)
} else{
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  script.dir  <- dirname(script.name)
}

source(paste0(script.dir,"/add_missing_exRxns.R"))
source(paste0(script.dir,"/generate_rxn_stoich_hash.R"))

# get options first
spec <- matrix(c(
  'blast.res', 'r', 1, "character", "Blast-results table generated by gapseq.sh.",
  'topology.evidence', 't', 2, "character", "File with space-separated list of reactions (SEED-namespace), that have gapseq-topology evidence.",
  'gram', 'g', 2, "character", "Gram \"pos\" OR \"neg\" OR \"auto\"? Default: \"auto\". Please note: if set to \"auto\", the external programms barrnap, usearch, and bedtools are required.",
  'model.name', 'n', 2, "character", "Name of draft model network. Default: the basename of \"blast.res\"",
  'genome.seq', 'c', 2, "character", "If gram is set to \"auto\", the genome sequence is required to search for 16S genes, which are used to predict gram-staining.",
  'high.evi.rxn.BS', "b", 2, "numeric", "Reactions with an associated blast-hit with a bitscore above this value will be added to the draft model as core reactions (i.e. high-sequence-evidence reactions)",
  'output.dir', 'o', 2, "character", "Directory to store results. Default: \".\" (alternatives not yet implemented)",
  'sbml.output', 's', 2, "logical", "Should the gapfilled model be saved as sbml? Default: FALSE (export not yet implemented)"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( is.null(opt$blast.res) | (is.null(opt$gram) & is.null(opt$genome.seq)) | (!is.null(opt$gram) && opt$gram=="auto" && is.null(opt$genome.seq))) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# Setting defaults if required
if ( is.null(opt$model.name) ) { opt$model.name = NA_character_ }
if ( is.null(opt$output.dir) ) { opt$output.dir = "." }
if ( is.null(opt$sbml.output) ) { opt$sbml.output = F }
if ( is.null(opt$high.evi.rxn.BS) ) { opt$high.evi.rxn.BS = 200 }
if ( is.null(opt$gram) ) { opt$gram = "auto" }
if ( is.null(opt$topology.evidence) ) { opt$topology.evidence = NA_character_ }

# Arguments:
blast.res         <- opt$blast.res
model.name        <- opt$model.name
gram              <- opt$gram
output.dir        <- opt$output.dir
genome.seq        <- opt$genome.seq
high.evi.rxn.BS   <- opt$high.evi.rxn.BS
topology.evidence <- opt$topology.evidence

if(is.na(model.name))
  model.name <- gsub("-all-blast.tbl","",basename(blast.res), fixed = T)

# loading transporters and reactions with topology evidence
topo.evi <- c()
for( file in unlist(str_split(topology.evidence, ",")) ){
  input.rxns <- readLines(file)
  topo.evi <- c(topo.evi, unlist(str_split(input.rxns, " ")))
}
topo.evi <- unique(topo.evi)
cat("\nLoading input reactions:", length(topo.evi), "\n")

# construct draft model
mod <- build_draft_model_from_blast_results(blast.res = blast.res, 
                                            topo.evi = topo.evi,
                                            gram = gram, 
                                            model.name = model.name, 
                                            genome.seq = genome.seq, 
                                            high.evi.rxn.BS = high.evi.rxn.BS,
                                            script.dir = script.dir)

# save draft model and reaction weights
saveRDS(mod$mod,file = paste0(model.name, ".RDS"))
saveRDS(mod$cand.rxns,file = paste0(model.name, "-rxnWeights.RDS"))

