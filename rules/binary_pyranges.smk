rule pybedtools_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/pybedtools_overlap.py"


rule pyranges_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/overlap.py"


rule bioconductor_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/overlap/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/overlap/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/overlap.R"


rule pybedtools_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/intersection_pybedtools.py"


rule pyranges_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/intersection.py"


rule bioconductor_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/intersection/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/intersection.R"


rule pybedtools_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/set_intersection_pybedtools.py"


rule pyranges_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/set_intersection.py"


rule bioconductor_set_intersection:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/set_intersection/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/set_intersection/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/set_intersection.R"


rule pyranges_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/subtract.py"


rule pybedtools_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/subtract_pybedtools.py"


rule bioconductor_subtract:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/subtract/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/subtract/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/subtraction.R"


rule pyranges_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/nearest.py"


rule pybedtools_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/pybedtools_nearest.py"


rule bioconductor_nearest:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/nearest.R"


rule pybedtools_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/pybedtools/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/pybedtools/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/pybedtools_nearest_nonoverlapping.py"


rule pyranges_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/pyranges_{num_cores}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/pyranges_{num_cores}/{iteration}_{size}_benchmark.txt"
    threads:
        4
    script:
        "../scripts/nearest_nonoverlapping.py"


rule bioconductor_nearest_nonoverlapping:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/nearest_nonoverlapping/bioconductor/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/nearest_nonoverlapping/bioconductor/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/nearest_nonoverlapping.R"
