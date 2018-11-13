import pandas as pd

import numpy as np
# import ray
# ray.init(num_cpus=2)

from pyranges import PyRanges

from time import time
import datetime

chip_f = snakemake.input.chip
background_f = snakemake.input.background

nrows = None

chip = pd.read_table(chip_f, sep="\t", usecols=[0, 1, 2, 5], header=None, nrows=nrows,
                     names="Chromosome Start End Strand".split(),
        dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})

start_init = time()
cgr = PyRanges(chip, copy_df=False)
end_init = time()

print(end_init - start_init)

background = pd.read_table(background_f, sep="\t", usecols=[0, 1, 2, 5], nrows=nrows,
                           header=None, names="Chromosome Start End Strand".split(),
        dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})

start_init = time()
bgr = PyRanges(background, copy_df=False)
end_init = time()

print(end_init - start_init)
start = time()
result = cgr.nearest(bgr, strandedness="same")

end = time()

print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

# result.df.to_csv(snakemake.output.result, sep=" ", index=False)
