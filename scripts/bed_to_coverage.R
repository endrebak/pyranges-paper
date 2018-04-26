library(GenomicRanges)
library(rtracklayer)
library(rbenchmark)

## chip = import()
f = snakemake@input[[1]],
stranded = snakemake@wildcards[["stranded"]]

df = read.table(f, sep="\t",
                col.names=c("Chromosome", "Start", "End", "Name", "Score", "Strand"),
                header=FALSE)

if (stranded == "stranded"){
  split_list = list(df$Chromosome, df$Strand)
} else {
  split_list = list(df$Chromosome)
}

start.time <- Sys.time()

chip <-
  lapply(split(df, split_list), function(i){
    GRanges(seqnames = i$Chromosome,
            ranges = IRanges(start = i$Start,
                             end = i$End))
  })

chip_list = lapply(chip, coverage)

end.time <- Sys.time()

time.taken <- end.time - start.time
