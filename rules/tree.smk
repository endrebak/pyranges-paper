rule tree_operation:
    input:
        chip = "{prefix}/data/download/chip_{size}.bed.gz",
        background = "{prefix}/data/download/input_{size}.bed.gz",
    output:
        time = "{prefix}/benchmark/{tree_operation}/{library}/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{tree_operation}/{library}/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{tree_operation}/{library}/{filetype}/{iteration}_{size}_benchmark.txt"
    wildcard_constraints:
        tree_operation = regex(tree_ops)
    params:
        build_code = lambda w: tree_map[w.library.replace("pyranges_1", "ncls")]["tree_build"],
        overlaps_code = lambda w: tree_map[w.library.replace("pyranges_1", "ncls")]["tree_overlap"]
    script:
        "../scripts/tree_operation.py"
