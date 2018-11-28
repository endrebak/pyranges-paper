


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
        "../scripts/chip_minus_input.py"


rule bioconductor_chip_minus_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_minus_input/bioconductor/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_minus_input/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/chip_minus_input.R"


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
        "../scripts/bed_to_coverage.py"


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
        "../scripts/bed_to_coverage.R"
