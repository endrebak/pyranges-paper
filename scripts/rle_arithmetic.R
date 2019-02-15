source("scripts/helpers.R")

operation = snakemake@params[["code"]]
operation = paste0("f <- function(c1p, c2p, c1m, c2m){\n", operation, "\nreturn(result)}")

eval(parse(text=operation))
print(operation)

filetype = snakemake@wildcards[["filetype"]]

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

c = file_to_grange(fc)
b = file_to_grange(fb, filetype)

cp = coverage(c[strand(c) == "+"])
cm = coverage(c[strand(c) == "-"])

bp = coverage(b[strand(b) == "+"])
bm = coverage(b[strand(b) == "-"])


start.time <- Sys.time()
result = f(cp, bp, cm, bm)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, file=snakemake@output[["time"]])

result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
