from helpers import read_file, file_to_grange, init_ray; init_ray(snakemake.wildcards)

from pyrle import PyRles

from time import time
import datetime
import numpy as np

chip = file_to_grange(snakemake.input[0], dtype=np.int64)

print(chip)
start = time()
grles = PyRles(chip, stranded=True)
end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output["timing"], "w+").write(minutes_seconds)

print(grles)
