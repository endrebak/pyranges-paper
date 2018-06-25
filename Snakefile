"Snakemake workflow to benchmark pyranges vs bioconductor genomicranges/s4vectors."

from time import time
import datetime

import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

prefix = "/mnt/scratch/endrebak/pyranges_benchmark"

test_files = {"chip": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/Histone_H3K27me3/Aorta/UCSD.Aorta.H3K27me3.STL003.bed.gz",
              "input": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/ChIP-Seq_Input/Aorta/UCSD.Aorta.Input.STL002.bed.gz"}

sizes = [int(f) for f in [1e6, 5e6, 1e7, 1.5e7]]

iterations = range(10)

libraries = "bioconductor pyranges".split()

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))


wildcard_constraints:
    chip = regex(test_files.keys()),
    size = regex(sizes),
    iteration = regex(iterations),
    libraries = regex(libraries)


bed_to_coverage_files = expand("{prefix}/benchmark/bed_to_coverage/{library}/{iteration}_{size}_time.txt",
                               prefix=prefix, iteration=iterations, size=sizes, library=libraries),


bed_to_granges_files = expand("{prefix}/benchmark/bed_to_granges/{library}/{iteration}_{size}_time.txt",
                              prefix=prefix, iteration=iterations, size=sizes, library=libraries),

chip_minus_input_files = expand("{prefix}/benchmark/chip_minus_input/{library}/{iteration}_{size}_time.txt",
                                prefix=prefix, iteration=iterations, size=sizes, library=libraries),

intersection_files = expand("{prefix}/benchmark/intersection/{library}/{iteration}_{size}_time.txt",
                            prefix=prefix, iteration=iterations, size=sizes, library=libraries)

nearest_files = expand("{prefix}/benchmark/nearest/{library}/{iteration}_{size}_time.txt",
                       prefix=prefix, iteration=iterations, size=sizes, library=libraries)

graph_files = [f"{prefix}/benchmark/graphs/time.png", f"{prefix}/benchmark/graphs/memory.png"]


rule all:
    input:
        graph_files


rule graphs:
    input:
        graph_files


rule intersection:
    input:
        intersection_files


rule rle:
    input:
        chip_minus_input_files


rule bed_to_granges:
    input:
        bed_to_granges_files


rule bed_to_coverage:
    input:
        bed_to_coverage_files


rule nearest:
    input:
        nearest_files



rule download_data:
    output:
        "{prefix}/data/download/{chip}.bed.gz"
    params:
        lambda w: test_files[w.chip]
    shell:
        "curl {params[0]} > {output[0]}"


# to test if sorting affects the results
rule shuffle:
    input:
        "{prefix}/data/download/{chip}.bed.gz"
    output:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    shell:
        "zcat {input[0]} | shuf -n {wildcards.size} | gzip > {output[0]}"


rule pyranges_bed_to_coverage:
    "How long it takes to turn a bed-file into a PyRles-object."
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/pyranges/{iteration}_{size}_time.txt",
        preview = "{prefix}/benchmark/bed_to_coverage/pyranges/{iteration}_{size}_preview.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/pyranges/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_coverage.py"


rule bioconductor_bed_to_coverage:
    "How long it takes to turn a bed-file into an RleList."
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_preview.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_coverage.R"


rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/pyranges/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/pyranges/{iteration}_{size}_benchmark.txt"
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
        "{prefix}/benchmark/chip_minus_input/pyranges/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_minus_input/pyranges/{iteration}_{size}_benchmark.txt"
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


rule pyranges_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pyranges/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/intersection/pyranges/{iteration}_{size}_result.txt"
    benchmark:
        "{prefix}/benchmark/intersection/pyranges/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/intersection.py"


rule bioconductor_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_result.txt"
    benchmark:
        "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/intersection.R"


rule pyranges_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/pyranges/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/nearest/pyranges/{iteration}_{size}_result.txt"
    benchmark:
        "{prefix}/benchmark/nearest/pyranges/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/nearest.py"


rule bioconductor_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_result.txt"
    benchmark:
        "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "scripts/nearest.R"


rule collect_times:
    input:
        bed_to_granges_files, bed_to_coverage_files, chip_minus_input_files, intersection_files, nearest_files
    output:
        "{prefix}/benchmark/collected_timings.txt"
    run:
        rowdicts = []
        for f in input:
            bmark_f = f.replace("time.txt", "benchmark.txt")
            function, library, timingfile = f.split("/")[-3:]
            iteration, size = timingfile.split("_")[:2]

            timing = open(f).readlines()[0].strip()

            if library == "pyranges":
                minutes, seconds, fraction = timing.split(".")
                minutes, seconds = int(minutes), int(seconds)
                seconds += minutes * 60

                timing = ".".join(str(s) for s in [seconds, fraction])

            max_rss = pd.read_table(bmark_f, sep="\t", usecols=[2], skiprows=1, squeeze=True, header=None).values[0] / 1024

            rowdict = {"Iteration": iteration, "MaxRSSGB": max_rss,
                       "Seconds": timing, "Function": function, "Library": library, "NBIntervals":
                       int(size)/int(1e6)}

            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts).sort_values("Function Library NBIntervals".split())
        column_order = "Function Library NBIntervals MaxRSSGB Seconds Iteration".split()

        df[column_order].to_csv(output[0], index=False, sep="\t")


rule graph_results:
    input:
        "{prefix}/benchmark/collected_timings.txt"
    output:
        "{prefix}/benchmark/graphs/time.png"
    script:
        "scripts/graph_time.R"


rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings.txt"
    output:
        "{prefix}/benchmark/graphs/memory.png"
    script:
        "scripts/graph_memory.R"
