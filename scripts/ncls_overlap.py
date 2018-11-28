
from ncls import NCLS


from time import time
import datetime
import pandas as pd
import numpy as np

c = snakemake.input.chip
b = snakemake.input.background

chip = pd.read_table(b, sep="\t",
                     usecols=[1, 2], header=None, names="Start End".split(), engine="c",
                     dtype={"Start": np.int64, "End": np.int64} )

background = pd.read_table(b, sep="\t",
                     usecols=[1, 2], header=None, names="Start End".split(), engine="c",
                     dtype={"Start": np.int64, "End": np.int64} )


tree = NCLS(background.Start.values, background.End.values, background.index.values)

start = time()
result = tree.all_overlaps_both(chip.Start.values, chip.End.values, chip.index.values)
end = time()

# print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
