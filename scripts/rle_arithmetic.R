source("scripts/helpers.R")

f = match.fun(snakemake@params[[1]])

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

chip = file_to_coverage(fc)
background = file_to_coverage(fb)

start.time <- Sys.time()
res = f(chip, background)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[[1]])
