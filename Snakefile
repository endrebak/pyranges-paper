"Snakemake workflow to benchmark pyranges vs genomicranges."

from time import time
import datetime

import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

prefix = "/mnt/scratch/endrebak/pyranges_benchmark"

test_files = {"chip": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/Histone_H3K27me3/Aorta/UCSD.Aorta.H3K27me3.STL003.bed.gz",
              "input": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/ChIP-Seq_Input/Aorta/UCSD.Aorta.Input.STL002.bed.gz"}

sizes = [int(f) for f in [1e6, 5e6, 1e7, 1.5e7]] #, 23547846]]


def regex(lst):

    return "({})".format("|".join([str(e) for e in lst]))


rule all:
    input:
        expand("{prefix}/benchmark/bed_to_coverage/{library}/{chip}_{size}_time.txt",
               prefix=prefix, chip=test_files, size=sizes, library="bioconductor pyranges".split()),
        expand("{prefix}/benchmark/bed_to_granges/{library}/{chip}_{size}_time.txt",
               prefix=prefix, chip=test_files, size=sizes, library="bioconductor pyranges".split()),
        expand("{prefix}/benchmark/chip_minus_input/{library}/{size}_time.txt",
               prefix=prefix, size=sizes, library="bioconductor pyranges".split()),
        expand("{prefix}/benchmark/intersection/{library}/{size}_time.txt",
               prefix=prefix, size=sizes, library="bioconductor pyranges".split())


rule intersection:
    input:
        expand("{prefix}/benchmark/intersection/{library}/{size}_time.txt",
               prefix=prefix, size=sizes, library="bioconductor pyranges".split())


rule rle:
    input:
        expand("{prefix}/benchmark/chip_minus_input/{library}/{size}_time.txt",
               prefix=prefix, size=sizes, library="bioconductor pyranges".split()),


rule bed_to_granges:
    input:
        expand("{prefix}/benchmark/bed_to_granges/{library}/{chip}_{size}_time.txt",
                chip=test_files,prefix=prefix, size=sizes, library="bioconductor pyranges".split()),

rule bed_to_coverage:
    input:
        expand("{prefix}/benchmark/bed_to_coverage/{library}/{chip}_{size}_time.txt",
               chip=test_files,prefix=prefix, size=sizes, library="bioconductor pyranges".split()),



wildcard_constraints:
    chip = regex(test_files.keys()),
    size = regex(sizes)


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
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/pyranges/{chip}_{size}_time.txt",
        preview = "{prefix}/benchmark/bed_to_coverage/pyranges/{chip}_{size}_preview.txt"
    threads:
        25
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/pyranges/{chip}_{size}_benchmark.txt"
    run:
        chip = pd.read_table(input[0], sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        start = time()
        grles = PyRles(chip, stranded=True)
        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%-M.%-S.%f')

        open(output["timing"], "w+").write(minutes_seconds)
        open(output["preview"], "w+").write(str(grles))



rule bioconductor_bed_to_coverage:
    "How long it takes to turn a bed-file into an RleList."
    input:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    output:
        timing = "{prefix}/benchmark/bed_to_coverage/bioconductor/{chip}_{size}_time.txt",
        result = "{prefix}/benchmark/bed_to_coverage/bioconductor/{chip}_{size}_preview.txt"
    threads:
        25
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/bioconductor/{chip}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_coverage.R"


rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/pyranges/{chip}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/pyranges/{chip}_{size}_benchmark.txt"
    run:
        chip = pd.read_table(input[0], sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

        start = time()
        granges = PyRanges(chip)
        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%-M.%-S.%f')

        open(output[0], "w+").write(minutes_seconds)



rule bioconductor_bed_to_GRanges:
    input:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{chip}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{chip}_{size}_benchmark.txt"
    script:
        "scripts/bed_to_GRanges.R"



rule pyranges_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/pyranges/{size}_time.txt"
    run:
        chip_f = input.chip
        background_f = input.background

        print("reading")
        nrows = None
        chip = pd.read_table(chip_f, sep="\t", nrows=nrows,
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        grles_chip = PyRles(chip, stranded=True)

        background = pd.read_table(background_f, sep="\t", nrows=nrows,
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        grles_background = PyRles(background, stranded=True)
        print("done reading")

        start = time()

        print("starting")
        result = grles_chip.sub(grles_background)
        print("ending")

        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%-M.%-S.%f')

        open(output[0], "w+").write(minutes_seconds)



rule bioconductor_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/bioconductor/{size}_time.txt"
    script:
        "scripts/chip_minus_input.R"


rule pyranges_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pyranges/{size}_time.txt",
        result = "{prefix}/benchmark/intersection/pyranges/{size}_result.txt"
    run:
        chip_f = input.chip
        background_f = input.background

        chip = pd.read_table(chip_f, sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        cgr = PyRanges(chip)

        background = pd.read_table(background_f, sep="\t",
                                   usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        bgr = PyRanges(background)

        start = time()
        print("starting")
        result = cgr.set_intersection(bgr, strandedness="same")

        end = time()

        print("ending")
        print(result)

        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%-M.%-S.%f')

        open(output.time, "w+").write(minutes_seconds)
        result.df.to_csv(output.result, sep=" ", index=False)





rule bioconductor_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/bioconductor/{size}_time.txt",
        result = "{prefix}/benchmark/intersection/bioconductor/{size}_result.txt"
    script:
        "scripts/intersection.R"


rule pyranges_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/nearest/pyranges/{size}_time.txt"
    run:
        chip_f = input.chip
        background_f = input.background

        chip = pd.read_table(chip_f, sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        cgr = PyRanges(chip)

        background = pd.read_table(background_f, sep="\t",
                                   usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        bgr = PyRanges(background)

        start = time()
        result = cgr.nearest(bgr)

        end = time()

        print(result)

        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%-M.%-S.%f')

        open(output[0], "w+").write(minutes_seconds)


# rule bioconductor_nearest:
#     input:
#         chip = "{prefix}/data/download/chip_{size}.bed.gz",
#         background = "{prefix}/data/download/input_{size}.bed.gz",
#     output:
#         "{prefix}/benchmark/intersection/bioconductor/{size}_time.txt"
#     script:
#         "scripts/chip_minus_input.R"
