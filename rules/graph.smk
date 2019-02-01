titles = {"tree": "bx-python vs. PyRanges NCLS",
          "rle": "BioConductor S4Vectors vs. PyRanges RLEs",
          "io": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "unary": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "binary": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges"}

method_names = {"tree": ", Tree Methods: ",
                "rle": ", RLE Methods: ",
                "io": ", IO Methods: ",
                "unary": ", Unary Methods: ",
                "binary": ", Binary Methods: "}


def get_title(w):

    c = w.category
    f = " (Filetype" + w.filetype + ")"

    return "Time Usage" + method_names[c] + f


rule graph_time:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/time_{filetype}_{category}.pdf"
    params:
        title = get_title,
        subset = False
    script:
        "../scripts/graph_time.R"



rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/memory_{filetype}_{category}.pdf"
    params:
        title = lambda w: get_title(w).replace("Time", "Memory"),
        subset = False
    script:
        "../scripts/graph_memory.R"


rule subset_collected_timings:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/main_paper_collected_timings_mean_{filetype}_{category}.txt"
    shell:
        'head -1 {input[0]} > {output[0]}; grep -P "(\\bintersect|nearest\\b|join|subtract)" {input[0]} | grep -v -P "(_2|_24|_48)" >> {output[0]}'



# rule graph_time_subset:
#     input:
#         "{prefix}/benchmark/main_paper_collected_timings_mean_{filetype}_{category}.txt"
#     output:
#         "{prefix}/benchmark/graphs/main_paper_time_{filetype}_{category}.pdf"
#     params:
#         title = get_title,
#         subset = True
#     script:
#         "../scripts/graph_time.R"


rule graph_paper:
    input:
        "{prefix}/benchmark/main_paper_collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/main_paper_{filetype}_{category}.pdf"
    params:
        title = get_title,
        subset = True
    script:
        "../scripts/graph_paper.R"
