"Snakemake workflow to benchmark pyranges vs bioconductor genomicranges/s4vectors."

from time import time
import datetime

import numpy as np
import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

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

wildcard_constraints:
    chip = "(chip|input)", # regex(test_files.keys()),
    size = regex(sizes),
    iteration = regex(iterations),
    libraries = regex(libraries),
    num_cores = regex([1, 2, 4, 8, 24, 48]),
    sorted = regex(sort)

# genomicranges_pyranges_only = "pyranges bioconductor".split()


bed_to_coverage_files = expand("{prefix}/benchmark/bed_to_coverage/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries, sorted=sort),


bed_to_granges_files = expand("{prefix}/benchmark/bed_to_granges/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries, sorted=sort),

chip_minus_input_files = expand("{prefix}/benchmark/chip_minus_input/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=no_pybedtools_libraries, sorted=sort),


intersection_files = expand("{prefix}/benchmark/intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, sorted=sort)
intersection_files_pybedtools = expand("{prefix}/benchmark/intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", sorted=sort)

overlap_files = expand("{prefix}/benchmark/overlap/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, sorted=sort)
overlap_files_pybedtools = expand("{prefix}/benchmark/overlap/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", sorted=sort)

set_intersection_files = expand("{prefix}/benchmark/set_intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, sorted=sort)
set_intersection_files_pybedtools = expand("{prefix}/benchmark/set_intersection/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", sorted=sort)

nearest_files = expand("{prefix}/benchmark/nearest{type}/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, sorted=sort, type=["", "_nonoverlapping"])
nearest_files_pybedtools = expand("{prefix}/benchmark/nearest{type}/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", sorted=sort, type=["", "_nonoverlapping"])

subtract_files = expand("{prefix}/benchmark/subtract/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=sizes, library=libraries, sorted=sort)
subtract_files_pybedtools = expand("{prefix}/benchmark/subtract/{library}/{iteration}_{size}_time.txt", prefix=prefix, iteration=iterations, size=pybedtool_sizes, library="pybedtools", sorted=sort)


graph_files = [f"{prefix}/benchmark/graphs/time.pdf", f"{prefix}/benchmark/graphs/memory.pdf"]


rule all:
    input:
        bed_to_granges_files,
        bed_to_coverage_files,
        chip_minus_input_files,
        intersection_files,
        intersection_files_pybedtools,
        overlap_files,
        overlap_files_pybedtools,
        set_intersection_files,
        set_intersection_files_pybedtools,
        nearest_files,
        nearest_files_pybedtools,
        subtract_files,
        subtract_files_pybedtools


rule graphs:
    input:
        graph_files


rule set_intersection:
    input:
        set_intersection_files


rule rle:
    input:
        chip_minus_input_files


rule subtract:
    input:
        subtract_files


rule bed_to_coverage:
    input:
        bed_to_coverage_files


rule nearest:
    input:
        nearest_files


rule fetch_chromsizes:
    output:
        "{prefix}/data/download/chromsizes.txt"
    shell:
        "fetchChromSizes hg38 | grep -v '_' > {output[0]}"


rule generate_data:
    input:
        "{prefix}/data/download/chromsizes.txt"
    output:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    shell:
        "bedtools random -n {wildcards.size} -g {input[0]} | gzip -9 > {output[0]}"


rule sort:
    input:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    output:
        "{prefix}/data/download/{chip}_{size}_sorted.bed.gz"
    shell:
        "zcat {input[0]} | sort -k1,1 -k2,2n | gzip > {output[0]}"


rule pyranges_bed_to_coverage:
    "How long it takes to turn a bed-file into a PyRles-object."
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/bed_to_coverage.py"


rule bioconductor_bed_to_coverage:
    "How long it takes to turn a bed-file into an RleList."
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_time.txt",
        preview = "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_preview.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_coverage.R"


rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/pyranges_{num_cores}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_GRanges.py"



rule bioconductor_bed_to_GRanges:
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_GRanges.R"


rule pyranges_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/pyranges_{num_cores}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_minus_input/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/chip_minus_input.py"


rule bioconductor_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/bioconductor/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_minus_input/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/chip_minus_input.R"




rule pybedtools_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/cluster_pybedtools.py"


rule pyranges_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/cluster.py"


rule bioconductor_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/cluster.R"



rule pybedtools_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/pybedtools_overlap.py"


rule pyranges_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/overlap.py"


rule bioconductor_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/overlap.R"


rule pybedtools_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/intersection_pybedtools.py"


rule pyranges_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/intersection.py"


rule bioconductor_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/intersection.R"


rule pybedtools_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/set_intersection_pybedtools.py"


rule pyranges_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/set_intersection.py"


rule bioconductor_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/set_intersection.R"


rule pyranges_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/subtract.py"


rule pybedtools_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/subtract_pybedtools.py"


rule bioconductor_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/subtraction.R"


rule pyranges_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/nearest.py"


rule pybedtools_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/pybedtools_nearest.py"


rule bioconductor_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/nearest.R"


rule pybedtools_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/pybedtools_nearest_nonoverlapping.py"


rule pyranges_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "scripts/nearest_nonoverlapping.py"


rule bioconductor_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/nearest_nonoverlapping.R"


rule collect_times:
    input:
        bed_to_coverage_files, chip_minus_input_files, set_intersection_files, nearest_files, intersection_files, subtract_files
    output:
        f"{prefix}/benchmark/collected_timings.txt"
    run:
        rowdicts = []
        for f in input:
            bmark_f = f.replace("time.txt", "benchmark.txt")
            function, library, timingfile = f.split("/")[-3:]
            iteration, size = timingfile.split("_")[:2]

            timing = open(f).readlines()[0].strip()

            if library in ["pyranges", "pybedtools"]:
                minutes, seconds, fraction = timing.split(".")
                minutes, seconds = int(minutes), int(seconds)
                seconds += minutes * 60

                timing = ".".join(str(s) for s in [seconds, fraction])

            max_rss = pd.read_table(bmark_f, sep="\t", usecols=[2], skiprows=1, squeeze=True, header=None).values[0] / 1024

            rowdict = {"Iteration": iteration, "MaxRSSGB": max_rss,
                       "Seconds": timing, "Function": function, "Library": library, "Log10NBIntervals":
                       np.log10(size)}

            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts).sort_values("Function Library NBIntervals".split())
        column_order = "Function Library NBIntervals MaxRSSGB Seconds Iteration".split()

        df[column_order].to_csv(output[0], index=False, sep="\t")


rule graph_results:
    input:
        "{prefix}/benchmark/collected_timings.txt"
    output:
        "{prefix}/benchmark/graphs/time.pdf"
    script:
        "scripts/graph_time.R"


rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings.txt"
    output:
        "{prefix}/benchmark/graphs/memory.pdf"
    script:
        "scripts/graph_memory.R"
