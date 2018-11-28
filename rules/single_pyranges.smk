


rule pyranges_dataframe_to_pyranges:
    input:
        correct_file
    output:
        "{prefix}/benchmark/dataframe_to_granges/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/dataframe_to_granges/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/dataframe_to_GRanges.py"



rule bioconductor_dataframe_to_GRanges:
    input:
        correct_file
    output:
        "{prefix}/benchmark/dataframe_to_granges/bioconductor/{filetype}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/dataframe_to_granges/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/dataframe_to_GRanges.R"


rule pybedtools_cluster:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/cluster/pybedtools/{filetype}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pybedtools/{filetype}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/pybedtools_cluster.py"


rule pyranges_cluster:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/cluster/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/cluster.py"


rule bioconductor_cluster:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/cluster/bioconductor/{filetype}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/cluster.R"


rule pyranges_start_end_sort:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/sort/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
    script:
        "../scripts/sort.py"


rule pybedtools_start_end_sort:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/sort/pybedtools/{filetype}/{iteration}_{size}_time.txt",
    script:
        "../scripts/pybedtools_sort.py"



rule pyranges_dataframe_to_coverage:
    "How long it takes to turn a dataframe-file into a PyRles-object."
    input:
        correct_file
    output:
        timing = "{prefix}/benchmark/dataframe_to_coverage/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/dataframe_to_coverage/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/dataframe_to_coverage.py"


rule bioconductor_dataframe_to_coverage:
    "How long it takes to turn a dataframe-file into an RleList."
    input:
        correct_file
    output:
        timing = "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_time.txt",
        preview = "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_preview.txt"
    benchmark:
        "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/dataframe_to_coverage.R"
