rule quicksect_build:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/tree_build/bx-python/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/tree_build/bx-python/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/quicksect_build.py"


rule quicksect_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/tree_overlap/bx-python/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/tree_overlap/bx-python/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/quicksect_overlap.py"



rule ncls_build:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/tree_build/pyranges_1/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/tree_build/pyranges_1/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/ncls_build.py"


rule ncls_overlap:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/tree_overlap/pyranges_1/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/tree_overlap/pyranges_1/{iteration}_{size}_benchmark.txt"
    script:
        "../scripts/ncls_overlap.py"
