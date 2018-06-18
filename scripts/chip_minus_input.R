library(GenomicRanges)
library(data.table)

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]


print("Reading data table")
cmd = paste0("zcat ", fc, " | cut -f 1-3,6")
cmd_bg= paste0("zcat ", fb, " | cut -f 1-3,6")
print(cmd)
chip_df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"))
input_df = fread(cmd_bg, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"))

print("creating granges")
chip = GRanges(seqnames = chip_df$Chromosome, ranges = IRanges(start = chip_df$Start, end = chip_df$End), strand=chip_df$Strand)
input = GRanges(seqnames = input_df$Chromosome, ranges = IRanges(start = input_df$Start, end = input_df$End), strand=input_df$Strand)

print("computing coverage")
chip = coverage(chip)
background = coverage(input)

start.time <- Sys.time()
res = chip - background
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[[1]])
