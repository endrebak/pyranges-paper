"Snakemake workflow to benchmark pyranges vs bioconductor genomicranges/s4vectors."

from time import time
import datetime
import glob

import yaml
import numpy as np
import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

from os import environ

ss = pd.read_table("sample_sheet.txt", sep=" ", header=0)

binary_map = yaml.load(open("supplementaries/binary.yaml"))
rle_map = yaml.load(open("supplementaries/rle.yaml"))
unary_map = yaml.load(open("supplementaries/unary.yaml"))
tree_map = yaml.load(open("supplementaries/tree.yaml"))

category_code = {"tree": tree_map, "unary": unary_map, "rle": rle_map, "binary": binary_map}

if not environ.get("TMUX"):
    raise Exception("Not using TMUX!")

import platform

prefix = "/mnt/scratch/endrebak/pyranges_benchmark"
iterations = [0, 1] # range(3)
libraries = "bioconductor pyranges_1 pyranges_2 pyranges_8 pyranges_24 pyranges_48 bx-python".split()
sizes = [int(f) for f in [1e5, 1e6, 1e7, 1e8]] #, 1e9, 1e10]]
pybedtool_sizes = [int(f) for f in [1e5, 1e6, 1e7]] #, 1e9, 1e10]]

# print(binary_map["pybedtools"])

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))

sort = ["sorted"]
num_cores = [1, 2, 4, 8, 24, 48]

small_run = True
if small_run:
    filetypes = "annotation"
    num_cores = [1]
    sizes = [int(f) for f in [1e5]] #, 1e9, 1e10]]
    pybedtool_sizes = [int(f) for f in [1e5]] #, 1e9, 1e10]]
    iterations = [0]

# print(list(unary_map["pybedtools"]) + list(unary_map["bioconductor"]))
# raise

# print(list(unary_map["pybedtools"]) + list(unary_map["bioconductor"]))
wildcard_constraints:
    filetype = regex("reads annotation".split()),
    unary_operation = regex(list(unary_map["pybedtools"]) + list(unary_map["bioconductor"])),
    operation = regex(list(binary_map["pybedtools"]) + list(binary_map["bioconductor"])),
    tree_operation = regex(list(tree_map["ncls"]) + list(tree_map["bx-python"])),
    rle_operation = regex(rle_map["pyranges"])
    # size = regex(sizes),
    # iteration = regex(iterations),
    # libraries = regex(libraries + ["ncls", "bx-python", "pybedtools"]),
    # num_cores = regex(num_cores),
    # subset = "(''|_subset)"

# genomicranges_pyranges_only = "pyranges bioconductor".split()
filetypes = "reads annotation".split()


def _expand(functions, path="{prefix}/benchmark/{function}/{library}/{filetype}/{iteration}_{size}_time.txt",
            sizes=sizes, iterations=iterations, filetypes=filetypes, num_cores=num_cores):

    outfiles = []
    for function in functions:
        # print("----" * 5)
        # print("function", function)
        libraries = list(ss[ss["Function"] == function].Library)
        category = ss[ss["Function"] == function].Category.iloc[0]
        # print(libraries)

        multicore = ss[(ss.Function == function) & (ss.Library == "pyranges")].Multicore
        assert len(multicore) < 2
        if len(multicore) == 1:
            multicore = multicore.iloc[0]
            # category = ss[(ss.Function == function) & (ss.Library == "pyranges")].Category.iloc[0]
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

        # if "dataframe" in function:
        #     print(function)
        #     print("larger_libs", larger_libs)
        #     print("smaller_libs", smaller_libs)
        # print(function)
        # print(libraries)
        # print(filetypes)

        _outfiles = expand(path, function=function, prefix=prefix, iteration=iterations, size=_sizes, library=larger_libs, filetype=filetypes) + \
            expand(path, function=function, prefix=prefix, iteration=iterations, size=pybedtool_sizes, library=smaller_libs, filetype=filetypes)

        # for f in _outfiles:
        #     if "dataframe" in f and "pybedtools" in f:
        #         print(f)
        outfiles.extend(_outfiles)

    # print(outfiles)
    return outfiles


single_pyranges_functions = ss[ss.Category == "single"].Function.drop_duplicates()
single_pyranges_files = _expand(single_pyranges_functions)

# for f in single_pyranges_files:
#     if "dataframe" in f and "pybedtools" in f:
#         print(f)

binary_pyranges_functions = ss[ss.Category == "binary"].Function.drop_duplicates()
binary_pyranges_files = _expand(binary_pyranges_functions)

rle_functions = ss[ss.Category == "rle"].Function.drop_duplicates()
rle_files = _expand(rle_functions)

tree_functions = ss[ss.Category == "tree"].Function.drop_duplicates()
tree_files = _expand(tree_functions)

category_dict = {"single": single_pyranges_files,
                 "binary": binary_pyranges_files,
                 "rle": rle_files,
                 "tree": tree_files}


# print(rle_map)
# raise
def correct_file(w):

    ft = w.filetype

    if ft == "reads":
        return "{prefix}/data/download/input_{size}.bed.gz".format(**w)
    else:
        return "{prefix}/data/download/annotation_{size}.gtf.gz".format(**w)


for rule in glob.glob("rules/*.smk"):
    print("including: " + rule, file=sys.stderr)

    # line to comment in if you want to rerun analyses with -F without regenerating the data
    # if "generate_data" in rule: continue

    # never want to redownload data
    if "download" in rule: continue

    include: rule


rule all:
    input:
        expand("{prefix}/benchmark/graphs/time_{filetype}_{category}.pdf", prefix=prefix, filetype=filetypes, category=category_dict)


rule generate:
    input:
        expand("{prefix}/data/download/annotation_{size}.gtf.gz", prefix=prefix, size=sizes),
        expand("{prefix}/data/download/chip_{size}.bed.gz", prefix=prefix, size=sizes),
        expand("{prefix}/data/download/background_{size}.bed.gz", prefix=prefix, size=sizes),
