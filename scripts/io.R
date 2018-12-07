source("scripts/helpers.R")

library(rtracklayer)

operation = snakemake@params[["io_operation"]]
operation = paste0("f <- function(gr){\n", operation, "\nreturn(result)}")
print(paste0("Performing operation ", operation, " for ", snakemake@wildcards[["unary_operation"]]))

file = snakemake@input[[1]]

eval(parse(text=operation))
start.time <- Sys.time()

result = f(gr)

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[["time"]])

result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
