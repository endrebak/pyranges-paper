import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

from time import time
import datetime


chip = pd.read_table(snakemake.input[0], sep="\t", usecols=[0, 1, 2, 5],
                     header=None, names="Chromosome Start End Strand".split(),
                     dtype={"Chromosome": "category", "Strand": "category"})
start = time()
grles = PyRles(chip, stranded=True)
end = time()
total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output["timing"], "w+").write(minutes_seconds)
open(snakemake.output["preview"], "w+").write(str(grles))
