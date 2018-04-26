library(GenomicRanges)
library(rtracklayer)
library(rbenchmark)

fc = snakemake@input[["chip"]],
fb = snakemake@input[["background"]],
stranded = snakemake@wildcards[["stranded"]]

chip = import(fc)
background = import(fb)

if (stranded == "stranded"){
  split_list = list(df$Chromosome, df$Strand)
} else {
  split_list = list(df$Chromosome)
}

start.time <- Sys.time()
res = cv_plus - cv_minus
end.time <- Sys.time()

time.taken <- end.time - start.time

chip <-
  lapply(split(df, split_list), function(i){
    GRanges(seqnames = i$Chromosome,
            ranges = IRanges(start = i$Start,
                             end = i$End))
  })

chip_list = lapply(chip, coverage)

end.time <- Sys.time()

time.taken <- end.time - start.time
