


rule pyranges_rle_arithmetic:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = correct_file,
    output:
        time = "{prefix}/benchmark/rle_{rle_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/rle_{rle_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/rle_{rle_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        code = lambda w: rle_map["pyranges"][w.rle_operation]
    script:
        "../scripts/rle_arithmetic.py"


rule bioconductor_rle_arithmetic:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = correct_file,
    output:
        time = "{prefix}/benchmark/rle_{rle_operation}/bioconductor/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/rle_{rle_operation}/bioconductor/{filetype}/{iteration}_{size}.result"
    benchmark:
        "{prefix}/benchmark/rle_{rle_operation}/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        code = lambda w: rle_map["bioconductor"][w.rle_operation]
    script:
        "../scripts/rle_arithmetic.R"
