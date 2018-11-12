library(GenomicRanges)
library(data.table)

## chip = import()
f = snakemake@input[[1]]

print("Reading data table")
if (Sys.info()["sysname"] == "Darwin"){
  cmd = paste0("gzcat ", f, " | cut -f 1-3,6")
} else {
  cmd = paste0("zcat ", f, " | cut -f 1-3,6")
}

print(cmd)
df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)

start.time <- Sys.time()
print("Starting to create GRanges")
gr = GRanges(seqnames = df$Chromosome, ranges = IRanges(start = df$Start, end = df$End), strand=df$Strand)

print("Starting to create Rles")

chip_list = coverage(gr)

end.time <- Sys.time()


time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[["timing"]])
capture.output(print(chip_list), file=snakemake@output[["preview"]])
