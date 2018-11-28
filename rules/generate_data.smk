

rule fetch_chromsizes:
    output:
        "{prefix}/data/download/chromsizes.txt"
    shell:
        "fetchChromSizes hg38 | grep -v '_' > {output[0]}"


rule generate_data:
    input:
        "{prefix}/data/download/chromsizes.txt"
    output:
        "{prefix}/data/download/{chip}_{size}.bed.gz"
    shell:
        "bedtools random -n {wildcards.size} -g {input[0]} | gzip -9 > {output[0]}"



rule download_genode:
    output:
        "{prefix}/data/download/annotation.gtf.gz"
    shell:
        "curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz > {output[0]}"


rule unzip_gencode:
    input:
        "{prefix}/data/download/annotation.gtf.gz"
    output:
        temp("{prefix}/data/download/annotation.gtf")
    shell:
        "zcat {input[0]} | grep -Pv '^#' > {output[0]}"


rule generate_random_annotation:
    input:
        "{prefix}/data/download/annotation.gtf"
    output:
        "{prefix}/data/download/annotation_{size}.gtf.gz"
    shell:
        "shuf -r -n {wildcards.size} {input[0]} | gzip > {output[0]}"
