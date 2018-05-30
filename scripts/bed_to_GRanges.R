library(GenomicRanges)
library(rtracklayer)
library(data.table)

f = snakemake@input[[1]]
stranded = snakemake@wildcards[["stranded"]]

cmd = paste0("zcat ", f, " | cut -f 1-3,6")
print(cmd)
df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"))

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

write(time.taken, snakemake@output[[1]])
