source("scripts/helpers.R")
library(GenomicRanges)
library(data.table)

f = snakemake@input[[1]]
df = get_df(f)

start.time <- Sys.time()

gr = GRanges(seqnames = df$Chromosome,
             ranges = IRanges(start = df$Start, end = df$End), strand = df$Strand)

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[[1]])
