library(GenomicRanges)
library(rtracklayer)
library(data.table)

## chip = import()
f = snakemake@input[[1]]
stranded = snakemake@wildcards[["stranded"]]

print("Reading data table")
cmd = paste0("zcat ", f, " | cut -f 1-3,6")
print(cmd)
df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"))

if (stranded == "stranded"){
  split_list = list(df$Chromosome, df$Strand)
} else {
  split_list = list(df$Chromosome)
}

start.time <- Sys.time()
print("Starting to create GRanges")
chip <-
  lapply(split(df, split_list), function(i){
    GRanges(seqnames = i$Chromosome,
            ranges = IRanges(start = i$Start,
                             end = i$End))
  })

print("Starting to create Rles")

chip_list = lapply(chip, coverage)

end.time <- Sys.time()


time.taken <- end.time - start.time

write(time.taken, snakemake@output[["timing"]])
capture.output(print(chip_list), file=snakemake@output[["result"]])
