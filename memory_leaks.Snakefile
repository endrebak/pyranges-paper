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

runs = [1, 10, 100, 1000]

def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))


wildcard_constraints:
    chip = regex(test_files.keys()),
    size = regex(sizes),
    runs = regex(runs)

functions = "bed_to_granges bed_to_coverage nearest intersection chip_minus_input".split()

all_files = expand("{prefix}/memcheck/{function}/{runs}.txt",
                   prefix=prefix, runs=runs, function=functions)

rule all:
    input:
        f"{prefix}/benchmark/graphs/memleak.png",
        all_files,
        f"{prefix}/memcheck/collected_results.txt"


rule shuffle:
    input:
        "{prefix}/data/download/{chip}.bed.gz"
    output:
        "{prefix}/data/download/{chip}_100000.bed.gz"
    shell:
        "zcat {input[0]} | shuf -n 100000 | gzip > {output[0]}"


rule pandas_baseline:
    input:
        chip = "{prefix}/data/download/chip_100000.bed.gz",
        background = "{prefix}/data/download/input_100000.bed.gz",
    output:
        "{prefix}/memcheck/pandas_baseline/1.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/pandas_mem.py {input[0]} {input[1]}"


rule pyranges_bed_to_coverage:
    input:
        "{prefix}/data/download/chip_100000.bed.gz"
    output:
        "{prefix}/memcheck/bed_to_coverage/{runs}.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/bed_to_GRanges_mem.py {input[0]} {wildcards.runs}"


rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/chip_100000.bed.gz"
    output:
        "{prefix}/memcheck/bed_to_granges/{runs}.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/bed_to_GRanges_mem.py {input[0]} {wildcards.runs}"


rule pyranges_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_100000.bed.gz",
        background = "{prefix}/data/download/input_100000.bed.gz",
    output:
        "{prefix}/memcheck/chip_minus_input/{runs}.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/chip_minus_input_mem.py {input[0]} {input[1]} {wildcards.runs}"


rule pyranges_intersection:
    input:
        chip = "{prefix}/data/download/chip_100000.bed.gz",
        background = "{prefix}/data/download/input_100000.bed.gz",
    output:
        "{prefix}/memcheck/intersection/{runs}.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/intersection_mem.py {input[0]} {input[1]} {wildcards.runs}"

rule pyranges_nearest:
    input:
        chip = "{prefix}/data/download/chip_100000.bed.gz",
        background = "{prefix}/data/download/input_100000.bed.gz",
    output:
        "{prefix}/memcheck/nearest/{runs}.txt"
    shell:
        "valgrind --log-file={output[0]} --tool=memcheck --leak-check=full --num-callers=30 --suppressions=valgrind-python.supp python scripts/nearest_mem.py {input[0]} {input[1]} {wildcards.runs}"


rule collect_files:
    input:
        pandas = "{prefix}/memcheck/pandas_baseline/1.txt",
        all_files = all_files
    output:
        "{prefix}/memcheck/collected_results.txt"
    run:
        def parse_valgrind(f):
            return [int(n.split(":")[1].split("by")[0].strip().replace(",", "")) for n in open(f).readlines()[-11:-8]]

        pdd, pdi, pdp = parse_valgrind(input.pandas)

        rowdicts = []
        for f in input.all_files:

            func, run = f.split("/")[-2:]
            run = run.split(".")[0]

            definite, indirectly, possibly = parse_valgrind(f)
            definite -= pdd
            indirectly -= pdi
            possibly -= pdp

            definite = max(definite, 0)
            indirectly = max(indirectly, 0)
            possibly = max(possibly, 0)

            for s, v in zip("Definite Possibly Indirectly".split(), [definite, possibly, indirectly]):
                rowdicts.append({"Function": func, "Iterations": run, "LeakType": s, "MegaBytes": v / (1024 * 1024)})

        df = pd.DataFrame.from_dict(rowdicts).sort_values("Function Iterations".split())
        column_order = "Function Iterations LeakType MegaBytes".split()

        df[column_order].to_csv(output[0], index=False, sep="\t")


rule graph_memleaks:
    input:
        "{prefix}/memcheck/collected_results.txt"
    output:
        "{prefix}/benchmark/graphs/memleak.png"
    script:
        "scripts/graph_memleaks.R"
