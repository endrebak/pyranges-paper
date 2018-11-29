source("scripts/helpers.R")

f = match.fun(snakemake@params[[1]])

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

c1p = file_to_coverage(fc, "+")
c1m = file_to_coverage(fc, "-")
c2p = file_to_coverage(fb, "+")
c2m = file_to_coverage(fb, "-")

start.time <- Sys.time()
resp = f(c1p, c2p)
resm = f(c1m, c2m)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[[1]])
