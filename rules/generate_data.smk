

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
