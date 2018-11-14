import numpy as np
import pandas as pd

import ray
num_cores = int(snakemake.wildcards.num_cores)
if num_cores != 1:
    ray.init(num_cpus=num_cores)
else:
    ray.init(local_mode=True, num_cpus=1) # logging_level=logging.CRITICAL # local_mode=True


from pyranges import PyRanges

from time import time
import datetime

chip_f = snakemake.input.chip
background_f = snakemake.input.background

chip = pd.read_table(chip_f, sep="\t",
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(),
                     dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})

cgr = PyRanges(chip, copy_df=False)

background = pd.read_table(background_f, sep="\t",
                            usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(),
        dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})
bgr = PyRanges(background, copy_df=False)

start = time()

result = cgr.subtraction(bgr, strandedness="same")

end = time()

print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

# result.df.to_csv(snakemake.output.result, sep=" ", index=False)
