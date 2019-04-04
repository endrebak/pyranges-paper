from helpers import read_file, file_to_grange, init_ray
init_ray(snakemake.wildcards)

from time import time
import datetime

import pyranges as pr

operation = snakemake.wildcards.write_operation
w = snakemake.wildcards

f = snakemake.input[0]

print(snakemake.params.code)
start = time()
exec(snakemake.params.code)
end = time()
total = end - start

print(result)

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output.time, "w+").write(minutes_seconds)

# result.df.to_csv("test.csv", sep="\t", index=False, header=True)

result_string = str(result)

open(snakemake.output.result, "w+").write(result_string)
