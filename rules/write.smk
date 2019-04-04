
rule bed_to_bam:
    input:
        bed = "{prefix}/data/download/{chip}_{size}.bed.gz",
        genome = "{prefix}/data/download/chromsizes.txt"
    output:
        "{prefix}/data/download/{chip}_{size}.bam"
    shell:
        "bedToBam -i {input.bed} -g {input.genome} > {output[0]}"


def correct_file_io(w):

    op = w.write_operation

    if op == "read_bed":
        return "{prefix}/data/download/input_{size}.bed.gz".format(**w)
    elif op == "read_gtf":
        return "{prefix}/data/download/annotation_{size}.gtf.gz".format(**w)
    else:
        return "{prefix}/data/download/input_{size}.bam".format(**w)




rule pyranges_io:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{write_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{write_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{write_operation}/pyranges_{num_cores}/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        code = lambda w: io_map["pyranges"][w.write_operation]
    wildcard_constraints:
        write_operation = regex(write_ops)
    script:
        "../scripts/io.py"


rule bioconductor:
    input:
        correct_file
    output:
        time = "{prefix}/benchmark/{write_operation}/bioconductor/{filetype}/{iteration}_{size}_time.txt",
        result = "{prefix}/benchmark/{write_operation}/bioconductor/{filetype}/{iteration}_{size}.result",
    benchmark:
        "{prefix}/benchmark/{write_operation}/bioconductor/{filetype}/{iteration}_{size}_benchmark.txt"
    params:
        operation = lambda w: io_map["bioconductor"][w.write_operation]
    wildcard_constraints:
        write_operation = regex(write_ops)
    script:
        "../scripts/io.R"
