# parse positional arguments:
# (1) Path to json file
# (2) Field to return (one of: "zenodoID", "version", "date")
args <- commandArgs(trailingOnly = TRUE)

res <- "NA"
tryCatch({
  jsf <- jsonlite::read_json(args[1])
  res <- jsf[[args[2]]][[1]][1]
  if(is.null(res))
    res <- "NA"
  if(res != "NA" & args[2] == "date")
    res <- substr(res, 1, 10)
}, error = function(e) { res <- "NA" })

cat(res)
