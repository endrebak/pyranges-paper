from helpers import read_file, init_ray; init_ray(snakemake.wildcards)
from pyranges import PyRanges

from time import time
import datetime

f = snakemake.input[0]
chip = read_file(f)

start = time()
granges = PyRanges(chip, copy_df=False)
end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
