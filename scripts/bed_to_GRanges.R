library(GenomicRanges)
library(data.table)

f = snakemake@input[[1]]

cmd = paste0("zcat ", f, " | cut -f 1-3,6")
print(cmd)
df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)

start.time <- Sys.time()

gr = GRanges(seqnames = df$Chromosome,
             ranges = IRanges(start = df$Start, end = df$End), strand = df$Strand)

## gr = GRanges(seqnames = df$V1,
##              ranges = IRanges(start = df$V2, end = df$V3), strand = df$V6)

## bg_gr = GRanges(seqnames = bg_df$V1,
##              ranges = IRanges(start = bg_df$V2, end = bg_df$V3), strand = bg_df$V6)

values(gr) = df$Name

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[[1]])
