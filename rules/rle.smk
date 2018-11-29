


rule pyranges_rle_arithmetic:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/rle_{operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/rle_{operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/rle_arithmetic.py"


rule bioconductor_rle_arithmetic:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        "{prefix}/benchmark/rle_{operation}/bioconductor/{filetype}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/rle_{operation}/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        lambda w: {"add": "+", "subtract": "-", "divide": "/", "multiply": "*"}[w.operation]
    script:
        "../scripts/rle_arithmetic.R"
