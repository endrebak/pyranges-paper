
from helpers import read_file, file_to_coverage, init_ray; init_ray(snakemake.wildcards)

from pyrle import PyRles

from time import time
import datetime
import numpy as np

chip = file_to_coverage(snakemake.input.chip)
background = file_to_coverage(snakemake.input.background)

operation = snakemake.wildcards.operation

pyop = {"add": "__add__", "subtract": "__sub__", "multiply": "__mul__", "divide": "__truediv__"}[operation]

m = getattr(chip, pyop)

start = time()

result = m(background)

end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

print(result)
