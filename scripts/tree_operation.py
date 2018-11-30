from ncls import NCLS
from quicksect import IntervalTree

from time import time
import datetime
import pandas as pd
import numpy as np

c = snakemake.input.chip
b = snakemake.input.background

df1 = pd.read_table(b, sep="\t",
                    usecols=[1, 2], header=None, names="Start End".split(), engine="c",
                    dtype={"Start": np.int32, "End": np.int32} )

df2 = pd.read_table(b, sep="\t",
                    usecols=[1, 2], header=None, names="Start End".split(), engine="c",
                    dtype={"Start": np.int32, "End": np.int32} )

if snakemake.wildcards.tree_operation == "tree_build":

    start = time()
    exec(snakemake.params.build_code)
    end = time()

else:

    exec(snakemake.params.build_code)

    start = time()
    exec(snakemake.params.overlaps_code)
    end = time()

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
