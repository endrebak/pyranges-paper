"Snakemake workflow to benchmark pyranges vs genomicranges."

from time import time
import datetime

from pyranges import GRanges
from pyrle import GRles

prefix = "/mnt/scratch/endrebak/pyranges_benchmark/"

test_files = {"chip": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/Histone_H3K27me3/Aorta/UCSD.Aorta.H3K27me3.STL003.bed.gz",
              "input": "ftp://ftp.genboree.org/EpigenomeAtlas/Current-Release/experiment-sample/ChIP-Seq_Input/Aorta/UCSD.Aorta.Input.STL002.bed.gz"}

def regex(lst):

    return "({})".format("|".join(lst))


rule wildcard_constraints:
    chip = regex(test_files.keys()),
    stranded = regex("stranded unstranded".split())


rule download_data:
    output:
        "{prefix}/data/download/{chip}.bed.gz"
    params:
        lambda w: test_files[w.chip]
    shell:
        "curl {params[0]} > {output[0]}"


rule pyranges_bed_to_coverage:
    "How long it takes to turn a bed-file into a GRles-object."
    input:
        "{prefix}/data/download/{chip}_{stranded}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_coverage/pyranges/{chip}_{stranded}_time.txt"
    threads:
        25
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/pyranges/{chip}_{stranded}_benchmark.txt"
    run:
        chip = pd.read_table(input[0], sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        stranded = True if wildcards.stranded == "stranded" else False
        start = time()
        grles = GRles(chip, n_jobs=25, stranded=stranded)
        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%M\t%S\n')

        open(output[0], "w+").write(minutes_seconds)


rule bioconductor_bed_to_coverage:
    "How long it takes to turn a bed-file into an RleList."
    input:
        "{prefix}/data/download/{chip}_{stranded}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_coverage/bioconductor/{chip}_{stranded}_time.txt"
    threads:
        25
    benchmark:
        "{prefix}/benchmark/bed_to_coverage/bioconductor/{chip}_{stranded}_benchmark.txt"
    script:
        "scripts/bed_to_coverage.R"


rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/{chip}_{stranded}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_pyranges/pyranges/{chip}_{stranded}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_pyranges/bioconductor/{chip}_{stranded}_benchmark.txt"
    run:
        chip = pd.read_table(input[0], sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

        stranded = True if wildcards.stranded == "stranded" else False
        start = time()
        granges = GRanges(chip, stranded=stranded)
        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%M\t%S\n')

        open(output[0], "w+").write(minutes_seconds)



rule bioconductor_bed_to_GRanges:
    input:
        "{prefix}/data/download/{chip}_{stranded}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_bioconductor/bioconductor/{chip}_{stranded}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_bioconductor/bioconductor/{chip}_{stranded}_benchmark.txt"
    script:
        "scripts/bed_to_GRanges.R"



rule pyranges_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip.bed.gz",
        background = "{prefix}/data/download/input.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/pyranges/{chip}_{stranded}_time.txt"
    run:
        chip_f = input.chip
        background_f = input.background
        stranded = wildcards.stranded

        chip = pd.read_table(chip_f, sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        stranded = True if wildcards.stranded == "stranded" else False
        grles_chip = GRles(chip, n_jobs=25, stranded=stranded)

        background = pd.read_table(background_f, sep="\t",
                             usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
        grles_background = GRles(background, n_jobs=25, stranded=stranded)

        start = time()

        result = grles_chip.sub(grles_background)

        end = time()
        total = end - start

        total_dt = datetime.datetime.fromtimestamp(total)

        minutes_seconds = total_dt.strftime('%M\t%S\n')

        open(output[0], "w+").write(minutes_seconds)
