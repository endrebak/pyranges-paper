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


def read_file(f, dtype=np.int32):

    ext = f.split(".")[-2]

    cols = d[ext]

    # print("using dtype: ", dtype)
    df = pd.read_table(f, sep="\t", usecols=cols,
                         header=None, names="Chromosome Start End Strand".split(),
                         dtype={"Chromosome": "category", "Strand": "category", "Start": dtype, "End": dtype})

    return df


def file_to_grange(f, dtype=np.int32, filetype="reads"):

    from pyranges import PyRanges, read_gtf

    if dtype == np.int64:
        extended = True
    else:
        extended = False

    if filetype == "reads":
        df = read_file(f, dtype)
        gr = PyRanges(df, extended=extended)
    elif filetype == "annotation":
        gr = read_gtf(f)

    return gr



def file_to_coverage(f, dtype=np.int64):

    from pyranges import PyRles
    gr = file_to_grange(f, dtype)
    return PyRles(gr, stranded=True)
