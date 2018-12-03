source("scripts/helpers.R")

library(GenomicRanges)

code = snakemake@params[["code"]]

print("code")
print(code)

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

gr1 = file_to_grange(fc)
## print(gr1)
gr2 = file_to_grange(fb, filetype=snakemake@wildcards[["filetype"]])

start.time <- Sys.time()
eval(parse(text=code))
end.time <- Sys.time()

time.taken = end.time - start.time
time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[["time"]])



result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
