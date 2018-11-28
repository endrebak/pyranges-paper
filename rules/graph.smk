d = {"tree": "bx-python vs. ncls",
     "rle": "BioConductor S4Vectors vs. PyRles",
     "I/O": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
     "single_pyranges": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
     "binary_pyranges": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges"}


def get_title(w):

    c = w.category

    return "Time Usage: " + d[c]


rule graph_results:
    input:
        "{prefix}/benchmark/collected_timings_mean_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/time_{category}.pdf"
    params:
        title = get_title
    script:
        "../scripts/graph_time.R"


rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings_mean_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/memory_{category}.pdf"
    params:
        title = lambda w: get_title[w].replace("Time", "Memory")
    script:
        "../scripts/graph_memory.R"
