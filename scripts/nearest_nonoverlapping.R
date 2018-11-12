
library(GenomicRanges)
library(data.table)

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

if (Sys.info()["sysname"] == "Darwin"){
  cmd = paste0("gzcat ", fc, " | cut -f 1-3,6")
  cmd_bg= paste0("gzcat ", fb, " | cut -f 1-3,6")
} else {
  cmd = paste0("zcat ", fc, " | cut -f 1-3,6")
  cmd_bg= paste0("zcat ", fb, " | cut -f 1-3,6")
}

print(cmd)
chip_df = fread(cmd, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)
input_df = fread(cmd_bg, header=FALSE, col.names=c("Chromosome", "Start", "End", "Strand"), stringsAsFactors=TRUE)

print("creating granges")
chip = GRanges(seqnames = chip_df$Chromosome, ranges = IRanges(start = chip_df$Start, end = chip_df$End), strand=chip_df$Strand)
input = GRanges(seqnames = input_df$Chromosome, ranges = IRanges(start = input_df$Start, end = input_df$End), strand=input_df$Strand)

print("finding nearest non-overlapping")

# https://support.bioconductor.org/p/110506/

fun2 <- function(x)
    min(x) == x & !duplicated(x)

distanceToNearest3 <- function(query, subject, fun = fun2, select=select) {
    hits <- union(
        precede(query, subject, select=select),
        follow(query, subject, select=select)
    )
    ## minimum distance by query
    distance <- distance(query[queryHits(hits)], subject[subjectHits(hits)])
    lst <- splitAsList(distance, queryHits(hits))
    keep <- unsplit(fun(lst), queryHits(hits))
    ## clean-up
    mcols(hits)$distance <- distance
    hits[keep]
}


start.time <- Sys.time()
result = distanceToNearest3(chip, input, select="all")
## result = chip[result]
end.time <- Sys.time()

time.taken <- end.time - start.time


time.taken <- as.numeric(time.taken, units="secs")

print(result)

write(time.taken, file=snakemake@output[["time"]])

## write.table(result, file=snakemake@output[["result"]], sep=" ", row.names=FALSE, quote=FALSE)
