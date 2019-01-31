rule pyranges_io:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{io_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{io_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{io_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        code = lambda w: io_map["pyranges"][w.io_operation]
    wildcard_constraints:
        io_operation = regex(io_ops)
    script:
        "../scripts/io.py"


rule bioconductor_io:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{io_operation}/bioconductor/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{io_operation}/bioconductor/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{io_operation}/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        operation = lambda w: io_map["bioconductor"][w.io_operation]
    wildcard_constraints:
        io_operation = regex(io_ops)
    script:
        "../scripts/io.R"
