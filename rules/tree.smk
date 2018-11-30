rule tree_operation:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/{tree_operation}/{library}/{filetype}/{iteration}_{size}_time.txt",
    benchmark:
        "{prefix}/benchmark/{tree_operation}/{library}/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        build_code = lambda w: tree_map[w.library]["tree_build"] if w.library == "bx-python" else tree_map["pyranges"]["tree_build"],
        overlaps_code = lambda w: tree_map[w.library]["tree_overlap"] if w.library == "bx-python" else tree_map["pyranges"]["tree_overlap"]
    script:
        "../scripts/tree_operation.py"
