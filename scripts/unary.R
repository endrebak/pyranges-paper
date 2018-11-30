source("scripts/helpers.R")

operation = snakemake@params[["operation"]]
operation = paste0("f <- function(gr){\n", operation, "\nreturn(result)}")
print(paste0("Performing operation ", operation, " for ", snakemake@wildcards[["unary_operation"]]))

fc = snakemake@input[[1]]

if (snakemake@wildcards[["unary_operation"]] != "dataframe_to_genomicrange"){
  gr = file_to_grange(fc)
} else {
  df = get_df(fc)
}

eval(parse(text=operation))
start.time <- Sys.time()

if (snakemake@wildcards[["unary_operation"]] != "dataframe_to_genomicrange"){
  result = f(gr)
} else {
  result = f(df)
}

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")

write(time.taken, snakemake@output[["time"]])

result_string = capture.output(print(result))
write(result_string, snakemake@output[["result"]])
