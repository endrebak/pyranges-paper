from helpers import read_file, file_to_grange, init_ray; init_ray(snakemake.wildcards)

from time import time
import datetime

operation = snakemake.wildcards.operation

f = snakemake.input[0]
f2 = snakemake.input[1]
gr = file_to_grange(f)
gr2 = file_to_grange(f2)

m = getattr(gr, operation)

start = time()
result = m(gr2, strandedness="same")
end = time()
total = end - start

print(result)

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
