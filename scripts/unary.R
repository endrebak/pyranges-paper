source("scripts/helpers.R")

operation = snakemake@params[["operation"]]
operation = paste0("f <- function(gr){\n", operation, "\nreturn(result)}")
print(paste0("Performing operation ", operation, " for ", snakemake@wildcards[["unary_operation"]]))

fc = snakemake@input[[1]]
gr = file_to_grange(fc)

eval(parse(text=operation))
start.time <- Sys.time()

f(gr)

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[[1]])
