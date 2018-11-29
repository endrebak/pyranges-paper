operation = snakemake@params[["operation"]]

fc = snakemake@input[["chip"]]
fb = snakemake@input[["background"]]

chip = file_to_grange(fc)
input = file_to_grange(fb)

print(paste0("Performing operation ", operation))

"operation is stored in m, behind the scenes"
eval(parse(text=operation))

start.time <- Sys.time()
result = m(chip, background)
end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken <- as.numeric(time.taken, units="secs")


write(time.taken, file=snakemake@output[["time"]])

print(head(result))
## write.table(result, file=snakemake@output[["result"]], sep=" ", row.names=FALSE, quote=FALSE)
