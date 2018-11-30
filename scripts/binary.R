source("scripts/helpers.R")

library(GenomicRanges)

operation = snakemake@params[["code"]]

"operation is stored in f, behind the scenes"
operation = paste0("f <- function(gr1, gr2){\n", operation, "\nreturn(result)}")
print(paste0("Performing operation ", operation, " for ", snakemake@wildcards[["operation"]]))
eval(parse(text=operation))
## f = match.fun()

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

chip = file_to_grange(fc)
background = file_to_grange(fb)

start.time <- Sys.time()
result = f(chip, background)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[[1]])


result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
