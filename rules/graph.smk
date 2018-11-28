
rule graph_results:
    input:
        "{prefix}/benchmark/collected_timings_mean_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/time.pdf"
    script:
        "../scripts/graph_time.R"


rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings_mean_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/memory.pdf"
    script:
        "../scripts/graph_memory.R"
