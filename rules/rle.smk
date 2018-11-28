


rule pyranges_chip_divide_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_divide_input/pyranges_{num_cores}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_divide_input/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/chip_divide_input.py"


rule bioconductor_chip_divide_input:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/chip_divide_input/bioconductor/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/chip_divide_input/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/chip_divide_input.R"
