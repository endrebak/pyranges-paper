import pandas as pd
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

grles_chip = PyRles(chip, stranded=True)

background = pd.read_table(background_f, sep="\t", nrows=nrows,
                           usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(),
                           dtype={"Chromosome": "category", "Strand": "category"})

grles_background = PyRles(background, stranded=True)
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
