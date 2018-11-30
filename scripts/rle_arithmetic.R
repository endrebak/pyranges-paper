source("scripts/helpers.R")

operation = snakemake@params[["code"]]
operation = paste0("f <- function(c1p, c2p, c1m, c2m){\n", operation, "\nreturn(result)}")

eval(parse(text=operation))
print(operation)


fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

c1p = file_to_coverage(fc, "+")
c1m = file_to_coverage(fc, "-")
c2p = file_to_coverage(fb, "+")
c2m = file_to_coverage(fb, "-")

start.time <- Sys.time()
result = f(c1p, c2p, c1m, c2m)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[["time"]])

result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
