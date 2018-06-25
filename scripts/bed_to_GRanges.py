import pandas as pd
from pyranges import PyRanges

from time import time
import datetime


chip = pd.read_table(snakemake.input[0], sep="\t",
                        usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

start = time()
granges = PyRanges(chip)
end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
