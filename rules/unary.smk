


# rule pyranges_dataframe_to_pyranges:
#     input:
#         correct_file
#     output:
#         "{prefix}/benchmark/dataframe_to_granges/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt"
#     benchmark:
#         "{prefix}/benchmark/dataframe_to_granges/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
#     script:
#         "../scripts/dataframe_to_GRanges.py"



# rule bioconductor_dataframe_to_GRanges:
#     input:
#         correct_file
#     output:
#         "{prefix}/benchmark/dataframe_to_granges/bioconductor/{filetype}/{iteration}_{size}_time.txt"
#     benchmark:
#         "{prefix}/benchmark/dataframe_to_granges/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
#     script:
#         "../scripts/dataframe_to_GRanges.R"

rule pyranges_unary:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{unary_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{unary_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{unary_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        code = lambda w: unary_map["pyranges"][w.unary_operation]
    wildcard_constraints:
        unary_operation = regex(unary_ops)
    script:
        "../scripts/unary.py"


rule bioconductor_unary:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{unary_operation}/bioconductor/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{unary_operation}/bioconductor/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{unary_operation}/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        operation = lambda w: unary_map["bioconductor"][w.unary_operation]
    wildcard_constraints:
        unary_operation = regex(unary_ops)
    script:
        "../scripts/unary.R"


rule pybedtools_unary:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{unary_operation}/pybedtools/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{unary_operation}/pybedtools/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{unary_operation}/pybedtools/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        operation = lambda w: unary_map["pybedtools"][w.unary_operation]
    wildcard_constraints:
        unary_operation = regex(unary_ops)
    script:
        "../scripts/pybedtools_unary.py"



# rule pyranges_dataframe_to_coverage:
#     "How long it takes to turn a dataframe-file into a PyRles-object."
#     input:
#         correct_file
#     output:
#         timing = "{prefix}/benchmark/dataframe_to_coverage/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
#     benchmark:
#         "{prefix}/benchmark/dataframe_to_coverage/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
#     threads:
#         4
#     script:
#         "../scripts/dataframe_to_coverage.py"


# rule bioconductor_dataframe_to_coverage:
#     "How long it takes to turn a dataframe-file into an RleList."
#     input:
#         correct_file
#     output:
#         timing = "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_time.txt",
#         preview = "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_preview.txt"
#     benchmark:
#         "{prefix}/benchmark/dataframe_to_coverage/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
#     script:
#         "../scripts/dataframe_to_coverage.R"
