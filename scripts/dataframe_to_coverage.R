source("scripts/helpers.R")
library(GenomicRanges)
library(data.table)

## chip = import()
f = snakemake@input[[1]]

gr = file_to_grange(f)

print("Starting to create Rles")

start.time <- Sys.time()
chip_list = coverage(gr)

end.time <- Sys.time()


time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[["timing"]])
capture.output(print(chip_list), file=snakemake@output[["preview"]])
