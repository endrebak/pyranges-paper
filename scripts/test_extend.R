

# Avoid rotating vector when using IRanges RleList

## I want to deduct the coverage in one RleList from that of another. A - B. But
## the lengths of the coverages
## aren't going to be exactly the same, which means that the vectors will be rotated. This is not what I want. If the Rle in chr1 in A is longer than that in B I want to extend B with the length difference from A.

library(GenomicRanges)
library(data.table)

fc = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/chip_1000000.bed.gz"
fb = "/mnt/scratch/endrebak/pyranges_benchmark/data/download/input_1000000.bed.gz"

print("Reading data table")
cmd = paste0("zcat ", fc, " | cut -f 1-3,6")
cmd_bg= paste0("zcat ", fb, " | cut -f 1-3,6")

chip_df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)
input_df = fread(cmd_bg, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)

print("creating granges")
chip = GRanges(seqnames = chip_df$Chromosome, ranges = IRanges(start = chip_df$Start, end = chip_df$End), strand=chip_df$Strand)
input = GRanges(seqnames = input_df$Chromosome, ranges = IRanges(start = input_df$Start, end = input_df$End), strand=input_df$Strand)

print("computing coverage")
chip = coverage(chip)
background = coverage(input)

outlist = list()
for (n in names(chip)){

  sumc = sum(runLength(chip[n]))
  sumb = sum(runLength(background[n]))

  # to avoid rotation
  if (sumc > sumb){
    print(paste0("Extending bg for chromosome ", n, " with ", sumc - sumb, " nucleotides"))
    bg = c(background[n], Rle(0, sumc - sumb))
    print("Computing chip")
    outlist[n] = chip[n] - bg
  } else if (sumb > sumc){
    print(paste0("Extending chip for chromosome ", n, " with ", sumb - sumc, " nucleotides"))
    ch = c(chip[n], Rle(0, sumb - sumc))
    outlist[n] = ch - background[n]
  }


}

print(paste0("length chr1 chip ", sum(runLength(chip["chr1"]))))
print(paste0("length chr1 bg ", sum(runLength(background["chr1"]))))

start.time <- Sys.time()
res = chip / background
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[[1]])
