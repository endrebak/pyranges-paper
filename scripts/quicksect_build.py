
from quicksect import IntervalNode, Interval, IntervalTree


from time import time
import datetime
import pandas as pd
import numpy as np

# c = snakemake.input.chip
b = snakemake.input.background


background = pd.read_table(b, sep="\t",
                     usecols=[1, 2], header=None, names="Start End".split(), engine="c",
                     dtype={"Start": np.int32, "End": np.int32} )

start = time()

tree = IntervalTree()
for start_, end_ in zip(background.Start, background.End):
    tree.add(start_, end_)

end = time()

# print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
