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

if platform.system() == "Darwin":
    prefix = "/Users/endrebakkenstovner/large_data/pyranges_paper"
    # iterations = [0] # range(10)
    iterations = [1] # range(10)
    sizes = [int(f) for f in [1e5, 1e6, 1e7]] # , 1e6, 1e7]]
    pybedtool_sizes = [int(f) for f in [1e5, 1e6, 1e7]] # , 1e6, 1e7]]
    libraries = "bioconductor pyranges_1 pyranges_2 pyranges_4 pybedtools".split()
else:
    prefix = "/mnt/scratch/endrebak/pyranges_benchmark"
    iterations = range(3)
    libraries = "bioconductor pyranges_1 pyranges_2 pyranges_8 pyranges_24 pyranges_48".split()
    sizes = [int(f) for f in [1e6, 1e7, 1e8]] #, 1e9, 1e10]]
    pybedtool_sizes = [int(f) for f in [1e6, 1e7]] #, 1e9, 1e10]]


no_pybedtools_libraries = "bioconductor pyranges_1 pyranges_2 pyranges_4".split()

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))

# sort = "sorted unsorted".split()
sort = ["sorted"]
num_cores = [1, 2, 4, 8, 24, 48]
wildcard_constraints:
    # chip = "(chip|input)", # regex(test_files.keys()),
    filetype = regex("reads annotation".split()),
    size = regex(sizes),
    iteration = regex(iterations),
    libraries = regex(libraries + ["ncls", "bx-python", "pybedtools"]),
    num_cores = regex(num_cores),
    sorted = regex(sort),
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

        outfiles = expand(path, function=functions, prefix=prefix, iteration=iterations, size=sizes, library=larger_libs, filetype=filetypes) + \
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

rule all:
    input:
        cluster_files
# cluster_files_pybedtools = expand("{prefix}/benchmark/cluster/{library}/{filetype}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools")
# file_to_granges_files = expand("{prefix}/benchmark/file_to_granges/{library}/{filetype}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries),
# single_pyranges_files = [bed_to_granges_files, cluster_files, cluster_files_pybedtools] # cluster, sort

# bed_to_coverage_files = expand("{prefix}/benchmark/bed_to_coverage/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries),


# chip_minus_input_files = expand("{prefix}/benchmark/chip_minus_input/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries),



# intersection_files = expand("{prefix}/benchmark/intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries)
# intersection_files_pybedtools = expand("{prefix}/benchmark/intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools")

# overlap_files = expand("{prefix}/benchmark/overlap/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries)
# overlap_files_pybedtools = expand("{prefix}/benchmark/overlap/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools")

# set_intersection_files = expand("{prefix}/benchmark/set_intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries)
# set_intersection_files_pybedtools = expand("{prefix}/benchmark/set_intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools")

# nearest_files = expand("{prefix}/benchmark/nearest{type}/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, type=["", "_nonoverlapping"])
# nearest_files_pybedtools = expand("{prefix}/benchmark/nearest{type}/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", type=["", "_nonoverlapping"])

# subtract_files = expand("{prefix}/benchmark/subtract/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries)
# subtract_files_pybedtools = expand("{prefix}/benchmark/subtract/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools")


# rle_files = [bed_to_coverage_files, chip_minus_input_files] # +, -, %


# tree_build = expand("{prefix}/benchmark/tree_build/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library=["ncls", "bx-python"])
# tree_overlap = expand("{prefix}/benchmark/tree_overlap/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library=["ncls", "bx-python"])

# tree_files = tree_build + tree_overlap

# binary_pyranges_files = sum([intersection_files, intersection_files_pybedtools,
#                              overlap_files, overlap_files_pybedtools, set_intersection_files,
#                              set_intersection_files_pybedtools, nearest_files, nearest_files_pybedtools,
#                              subtract_files, subtract_files_pybedtools], [])

# io_files = []

# category_dict = {}
# category_dict["rle"] = rle_files
# category_dict["single_pyranges"] = single_pyranges_files
# category_dict["binary_pyranges"] = binary_pyranges_files
# category_dict["tree"] = tree_files
# category_dict["io"] = io_files

# all_files = [bed_to_granges_files, bed_to_coverage_files, chip_minus_input_files, intersection_files, intersection_files_pybedtools, overlap_files, overlap_files_pybedtools, set_intersection_files, set_intersection_files_pybedtools, nearest_files, nearest_files_pybedtools, subtract_files, subtract_files_pybedtools]

# graph_files = [f"{prefix}/benchmark/graphs/time.pdf", f"{prefix}/benchmark/graphs/memory.pdf"]


def correct_file(w):

    ft = w.filetype

    if ft == "reads":
        return "{prefix}/data/download/input_{size}.bed.gz".format(**w)
    else:
        return "{prefix}/data/download/annotation_{size}.gtf.gz".format(**w)


for rule in glob.glob("rules/*.smk"):
    include: rule

# # rule all:
# #     input:
# #         all_files,
# #         expand("{prefix}/benchmark/graphs/time.pdf", prefix=prefix),
# #         expand("{prefix}/benchmark/graphs/memory.pdf", prefix=prefix)



# rule graphs:
#     input:
#         graph_files


# rule set_intersection:
#     input:
#         set_intersection_files


# rule rle:
#     input:
#         chip_minus_input_files


# rule subtract:
#     input:
#         subtract_files


# rule bed_to_coverage:
#     input:
#         bed_to_coverage_files


# rule nearest:
#     input:
#         nearest_files

# rule tree:
#     input:
#         tree_files
