library(GenomicRanges)
library(rtracklayer)
library(rbenchmark)

f = snakemake@input[[1]],
stranded = snakemake@wildcards[["stranded"]]

df = read.table(f, sep="\t",
                col.names=c("Chromosome", "Start", "End", "Name", "Score", "Strand"),
                header=FALSE)

start.time <- Sys.time()

if (stranded == "stranded"){
  gr = GRanges(seqnames = df$Chromosome,
              ranges = IRanges(start = df$Start, end = df$End), strand = df$Strand)
} else {
  gr = GRanges(seqnames = df$Chromosome,
              ranges = IRanges(start = df$Start, end = df$End))
}

values(gr) = df$Name

end.time <- Sys.time()

time.taken <- end.time - start.time
