
import numpy as np
import pandas as pd
from pyranges import PyRanges

from time import time
import datetime

chip_f = snakemake.input.chip
background_f = snakemake.input.background

chip = pd.read_table(chip_f, sep="\t",
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(), engine="c",
                     dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})


cgr = PyRanges(chip, copy_df=False)

background = pd.read_table(background_f, sep="\t",
                           usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split(), engine="c",
                           dtype= {"Start": np.int32, "End": np.int32, "Chromosome": "category", "Strand": "category"})

bgr = PyRanges(background, copy_df=False)

start = time()

result = cgr.overlap(bgr, strandedness="same")

end = time()

print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

# result.df.to_csv(snakemake.output.result, sep=" ", index=False)
