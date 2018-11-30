

rule download_genode:
    output:
        "{prefix}/data/download/annotation.gtf.gz"
    shell:
        "curl ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz > {output[0]}"
