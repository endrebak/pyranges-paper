
from helpers import read_file, file_to_coverage, init_ray; init_ray(snakemake.wildcards)

from pyrle import PyRles

from time import time
import datetime
import numpy as np

c1 = file_to_coverage(snakemake.input.chip)
c2 = file_to_coverage(snakemake.input.background)

code = snakemake.params.code

# pyop = {"add": "__add__", "subtract": "__sub__", "multiply": "__mul__", "divide": "__truediv__"}[operation]

# m = getattr(chip, pyop)

start = time()

# result = m(background)
exec(code)

end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

result_string = str(result)

open(snakemake.output.result, "w+").write(result_string)
