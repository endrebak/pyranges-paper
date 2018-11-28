source("helpers.R")
library(GenomicRanges)
library(data.table)

stop("Missing! Need to merge on strand, remember!")

f = snakemake@input[[1]]

gr = file_to_grange(f)

start.time <- Sys.time()



end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[[1]])
