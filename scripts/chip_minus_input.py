import pandas as pd


import ray
num_cores = int(snakemake.wildcards.num_cores)
if num_cores != 1:
    ray.init(num_cpus=num_cores)
else:
    ray.init(local_mode=True, num_cpus=1) # logging_level=logging.CRITICAL # local_mode=True


from pyranges import PyRanges
from pyrle import PyRles

from time import time
import datetime


chip_f = snakemake.input.chip
background_f = snakemake.input.background

print("reading")
nrows = None
chip = pd.read_table(chip_f, sep="\t", nrows=nrows,
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(),
                     dtype={"Chromosome": "category", "Strand": "category"})
print(chip)

grles_chip = PyRles(chip, stranded=True)
print(grles_chip)

background = pd.read_table(background_f, sep="\t", nrows=nrows,
                           usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(),
                           dtype={"Chromosome": "category", "Strand": "category"})

grles_background = PyRles(background, stranded=True)
print(grles_background)
print("done reading")

start = time()

print("starting")
result = grles_chip.div(grles_background)
print("ending")

end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

print(result)
