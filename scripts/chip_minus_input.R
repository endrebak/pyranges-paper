library(GenomicRanges)
library(rtracklayer)

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]
stranded = snakemake@wildcards[["stranded"]]

chip = coverage(import(fc))
background = coverage(import(fb))

start.time <- Sys.time()
res = chip - background
end.time <- Sys.time()

time.taken <- end.time - start.time

write(time.taken, file=snakemake@output[[1]])
