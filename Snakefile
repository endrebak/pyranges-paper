"Snakemake workflow to benchmark pyranges vs bioconductor genomicranges/s4vectors."

from time import time
import datetime
import glob

import numpy as np
import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

from os import environ

ss = pd.read_table("sample_sheet.txt", sep=" ", header=0)
print(ss)

if not environ.get("TMUX"):
    raise Exception("Not using TMUX!")

import platform

prefix = "/mnt/scratch/endrebak/pyranges_benchmark"
iterations = range(3)
libraries = "bioconductor pyranges_1 pyranges_2 pyranges_8 pyranges_24 pyranges_48 bx-python".split()
sizes = [int(f) for f in [1e6, 1e7, 1e8]] #, 1e9, 1e10]]
pybedtool_sizes = [int(f) for f in [1e6, 1e7]] #, 1e9, 1e10]]

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))

sort = ["sorted"]
num_cores = [1, 2, 4, 8, 24, 48]

wildcard_constraints:
    filetype = regex("reads annotation".split()),
    # size = regex(sizes),
    # iteration = regex(iterations),
    # libraries = regex(libraries + ["ncls", "bx-python", "pybedtools"]),
    # num_cores = regex(num_cores),
    # subset = "(''|_subset)"

# genomicranges_pyranges_only = "pyranges bioconductor".split()
filetypes = "reads annotation".split()


def _expand(functions):

    path = "{prefix}/benchmark/{function}/{library}/{filetype}/{iteration}_{size}_time.txt"

    outfiles = []
    for function in functions:
        libraries = list(ss[ss["Function"] == function].Library)

        multicore = ss[(ss.Function == function) & (ss.Library == "pyranges")].Multicore
        print(multicore)
        assert len(multicore) < 2
        if len(multicore) == 1:
            multicore = multicore.iloc[0]
            category = ss[(ss.Function == function) & (ss.Library == "pyranges")].Category.iloc[0]
        else:
            multicore = 0

        libraries = [l for l in libraries if l != "pyranges"]
        if multicore:
            libraries += ["pyranges_{}".format(n) for n in num_cores]
        else:
            libraries += ["pyranges_1"]


        smaller_libs = set("bx-python ncls pybedtools".split())
        larger_libs = set(libraries) - smaller_libs
        smaller_libs = set(libraries) & smaller_libs

        if category == "tree":
            _sizes = sizes[:-1]
        else:
            _sizes = sizes

        outfiles = expand(path, function=functions, prefix=prefix, iteration=iterations, size=_sizes, library=larger_libs, filetype=filetypes) + \
            expand(path, function=functions, prefix=prefix, iteration=iterations, size=pybedtool_sizes, library=smaller_libs, filetype=filetypes)

        return outfiles


single_pyranges_functions = ss[ss.Category == "single"].Function
single_pyranges_files = _expand(single_pyranges_functions)

binary_pyranges_functions = ss[ss.Category == "binary"].Function
binary_pyranges_files = _expand(binary_pyranges_functions)

rle_functions = ss[ss.Category == "rle"].Function
rle_files = _expand(rle_functions)

tree_functions = ss[ss.Category == "tree"].Function
tree_files = _expand(tree_functions)

category_dict = {"single": single_pyranges_files,
                 "binary": binary_pyranges_files,
                 "rle": rle_files,
                 "tree": tree_files}


def correct_file(w):

    ft = w.filetype

    if ft == "reads":
        return "{prefix}/data/download/input_{size}.bed.gz".format(**w)
    else:
        return "{prefix}/data/download/annotation_{size}.gtf.gz".format(**w)


for rule in glob.glob("rules/*.smk"):
    print("including: " + rule, file=sys.stderr)
    include: rule
