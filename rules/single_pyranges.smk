

rule pyranges_bed_to_pyranges:
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/pyranges_{num_cores}/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/bed_to_GRanges.py"



rule bioconductor_bed_to_GRanges:
    input:
        "{prefix}/data/download/chip_{size}.bed.gz"
    output:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{iteration}_{size}_time.txt"
    benchmark:
        "{prefix}/benchmark/bed_to_granges/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/bed_to_GRanges.R"


rule pybedtools_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/cluster_pybedtools.py"


rule pyranges_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/cluster.py"


rule bioconductor_cluster:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/cluster/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/cluster/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/cluster.R"
