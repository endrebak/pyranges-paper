titles = {"tree": "bx-python vs. PyRanges NCLS",
          "rle": "BioConductor S4Vectors vs. PyRanges RLEs",
          "I/O": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "single": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "binary": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges"}

method_names = {"tree": ", Tree Methods: ",
                "rle": ", RLE Methods: ",
                "I/O": ", IO Methods: ",
                "single": ", Unary Methods: ",
                "binary": ", Binary Methods: "}


def get_title(w):

    c = w.category
    f = " (Filetype" + w.filetype + ")"

    return "Time Usage" + method_names[c] + f # + titles[c]


rule graph_results:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/time_{filetype}_{category}.pdf"
    params:
        title = get_title
    script:
        "../scripts/graph_time.R"


rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/memory_{filetype}_{category}.pdf"
    params:
        title = lambda w: get_title[w].replace("Time", "Memory")
    script:
        "../scripts/graph_memory.R"
