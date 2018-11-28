import ray

def init_ray(w):
    num_cores = int(w.num_cores)
    if num_cores != 1:
        ray.init(num_cpus=num_cores)
    else:
        ray.init(local_mode=True, num_cpus=1) # logging_level=logging.CRITICAL # local_mode=True

import pandas as pd
import numpy as np

d = {"bed": [0, 1, 2, 5], "gtf": [0, 3, 4, 6]}


def read_file(f):

    ext = f.split(".")[-2]

    cols = d[ext]

    df = pd.read_table(f, sep="\t", usecols=cols,
                         header=None, names="Chromosome Start End Strand".split(),
                         dtype={"Chromosome": "category", "Strand": "category", "Start": np.int32, "End": np.int32})

    return df


def file_to_grange(f):

    from pyranges import PyRanges
    df = read_file(f)

    return PyRanges(df)
