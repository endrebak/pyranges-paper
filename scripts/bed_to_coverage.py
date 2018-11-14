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

print(grles)
