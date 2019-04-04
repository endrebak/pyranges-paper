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

ss = pd.read_csv("sample_sheet.txt", sep="\t", header=0, comment="#")

binary_map = yaml.load(open("supplementaries/binary.yaml"))
rle_map = yaml.load(open("supplementaries/rle.yaml"))
unary_map = yaml.load(open("supplementaries/unary.yaml"))
io_map = yaml.load(open("supplementaries/io.yaml"))
write_map = yaml.load(open("supplementaries/write.yaml"))
tree_map = yaml.load(open("supplementaries/tree.yaml"))

category_code = {"tree": tree_map, "unary": unary_map, "rle": rle_map, "binary": binary_map, "io": io_map}

if not environ.get("TMUX"):
    raise Exception("Not using TMUX!")

# import platform

prefix = "/mnt/scratch/endrebak/pyranges_benchmark"
iterations = [0, 1] # range(3)
libraries = "bioconductor pyranges_1 pyranges_4 pyranges_8 bx-python".split() #pyranges_2  pyranges_24 pyranges_48
sizes = [int(f) for f in [1e5, 1e6, 1e7]] #, 1e8]] #, 1e9, 1e10]]
pybedtool_sizes = [int(f) for f in [1e5, 1e6, 1e7]] #, 1e9, 1e10]]

# print(binary_map["pybedtools"])

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))

sort = ["sorted"]
num_cores = [1, 4, 8] # 2, 24, 48]

filetypes = "annotation reads".split()

small_run = False
if small_run:
    filetypes = "annotation"
    num_cores = [1]
    sizes = [int(f) for f in [1e5]] #, 1e9, 1e10]]
    pybedtool_sizes = [int(f) for f in [1e5]] #, 1e9, 1e10]]
    iterations = [0]

unary_ops = list(unary_map["pybedtools"]) + list(unary_map["bioconductor"])
binary_ops = list(binary_map["pybedtools"]) + list(binary_map["bioconductor"])
tree_ops = list(tree_map["ncls"]) + list(tree_map["bx-python"])
rle_ops = list(rle_map["pyranges"])
io_ops = list(io_map["pyranges"])
write_ops = list(write_map["pyranges"])


wildcard_constraints:
    filetype = regex("reads annotation".split()),
    unary_operation = regex(unary_ops),
    operation = regex(binary_ops),
    tree_operation = regex(tree_ops),
    rle_operation = regex(rle_ops),
    io_operation = regex(io_ops),
    libraries = regex(libraries + ["ncls", "bx-python", "pybedtools"]),
    category = regex(ss.Category.drop_duplicates().tolist()),
    num_cores = regex(num_cores),
    function = regex(ss.Function.drop_duplicates()),
    measure = "(time|memory)"


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


        if "read_" not in function:
            _outfiles = expand(path, function=function, prefix=prefix, iteration=iterations, size=_sizes, library=larger_libs, filetype=filetypes) + \
                expand(path, function=function, prefix=prefix, iteration=iterations, size=pybedtool_sizes, library=smaller_libs, filetype=filetypes)
        elif function == "read_gtf":
            _outfiles = expand(path, function=function, prefix=prefix, iteration=iterations, size=_sizes, library=larger_libs, filetype="annotation")
        elif function in ["read_bed", "read_bam"]:
            _outfiles = expand(path, function=function, prefix=prefix, iteration=iterations, size=_sizes, library=larger_libs, filetype="reads")

        outfiles.extend(_outfiles)

    return outfiles


single_pyranges_functions = ss[ss.Category == "unary"].Function.drop_duplicates()
single_pyranges_files = _expand(single_pyranges_functions)
# print(single_pyranges_files)
# raise

binary_pyranges_functions = ss[ss.Category == "binary"].Function.drop_duplicates()
binary_pyranges_files = _expand(binary_pyranges_functions)

rle_functions = ss[ss.Category == "rle"].Function.drop_duplicates()
rle_files = _expand(rle_functions)

tree_functions = ss[ss.Category == "tree"].Function.drop_duplicates()
tree_files = _expand(tree_functions)

io_functions = ss[ss.Category == "io"].Function.drop_duplicates()
io_files = _expand(io_functions)

category_dict = {"unary": single_pyranges_files,
                 "binary": binary_pyranges_files,
                 "rle": rle_files,
                 "tree": tree_files,
                 "io": io_files}


def correct_file(w):

    ft = w.filetype

    if ft == "reads":
        return "{prefix}/data/download/input_{size}.bed.gz".format(**w)
    else:
        return "{prefix}/data/download/annotation_{size}.gtf.gz".format(**w)

extensions = "pdf png".split()
time_files = expand("{prefix}/benchmark/graphs/time_{filetype}_{category}.{extension}", prefix=prefix, filetype=filetypes, category=category_dict, extension=extensions)
memory_files = expand("{prefix}/benchmark/graphs/memory_{filetype}_{category}.{extension}", prefix=prefix, filetype=filetypes, category=category_dict, extension=extensions)
main_paper_graphs = expand("{prefix}/benchmark/graphs/main_paper_{filetype}_{category}.{extension}", measure="time memory".split(), prefix=prefix, filetype=filetypes, category="binary", extension=extensions)
time_mem_together_graphs = expand("supplementary_paper/time_memory_together_{function}.{extension}",
       measure="time memory".split(), prefix=prefix, filetype=filetypes,
       function=ss.Function.drop_duplicates(), extension="png")

all_mds = expand("supplementary_paper/{function}_all.md", function=ss.Function.drop_duplicates())

for rule in glob.glob("rules/*.smk"):
    # print("including: " + rule, file=sys.stderr)

    # line to comment in if you want to rerun analyses with -F without regenerating the data
    # if "generate_data" in rule: continue

    # never want to redownload data
    if "download" in rule: continue

    include: rule



rule all:
    input:
        time_files, memory_files, main_paper_graphs, time_mem_together_graphs, all_mds

rule differences:
    input:
        expand("{prefix}/benchmark/differences/{num_cores}_{filetype}_{category}_differences.txt",
               num_cores = sorted(num_cores), filetype=filetypes, category="unary binary rle".split(), prefix=prefix),
        expand("{prefix}/benchmark/differences/{num_cores}_{filetype}_{category}_differences.txt",
               num_cores = [1], filetype=filetypes, category="tree io".split(), prefix=prefix)

# f = "supplementary_paper/{category}_all.md"
rule supplementary:
    input:
        "supplementary_paper/README.pdf"

        # expand(f, category=category_dict)

# rule supplementary_code:
#     input:
#         expand("{prefix}/benchmark/supplementaries/{category}.pdf", prefix=prefix, category=category_dict)


# rule supplementary_graphs:
#     input:
#         expand("supplementary_paper/{measure}.md", measure="time memory".split())



rule generate:
    input:
        expand("{prefix}/data/download/annotation_{size}.gtf.gz", prefix=prefix, size=sizes),
        expand("{prefix}/data/download/chip_{size}.bed.gz", prefix=prefix, size=sizes),
        expand("{prefix}/data/download/background_{size}.bed.gz", prefix=prefix, size=sizes),
