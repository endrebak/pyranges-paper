titles = {"tree": "bx-python vs. PyRanges NCLS",
          "rle": "BioConductor S4Vectors vs. PyRanges RLEs",
          "io": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "unary": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges",
          "binary": "BioConductor GenomicRanges vs. pybedtools vs. PyRanges"}

method_names = {"tree": "Tree Methods",
                "rle": "RLE Methods",
                "io": "IO Methods",
                "unary": "Unary Methods",
                "binary": "Binary Methods"}


def get_title(w):

    c = w.category
    f = " (" + w.filetype.capitalize() + ")"

    return "Time Usage " + method_names[c] + f


rule graph_time:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/time_{filetype}_{category}.{extension}"
    params:
        title = get_title,
        subset = False
    script:
        "../scripts/graph_time.R"



rule graph_memory:
    input:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/memory_{filetype}_{category}.{extension}"
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




rule graph_paper:
    input:
        "{prefix}/benchmark/main_paper_collected_timings_mean_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/graphs/main_paper_{filetype}_{category}.{extension}"
    params:
        title = get_title,
        subset = True
    script:
        "../scripts/graph_paper.R"


rule collect_filetype_split_function:
    input:
        lambda w: expand("{{prefix}}/benchmark/collected_timings_mean_{filetype}_{category}.txt", filetype=filetypes, category=ss[ss.Function == w.function].Category)
    output:
        "{prefix}/benchmark/collected_timings_mean_{function}.txt"
    run:
        dfs = []
        for f in input:
            filetype = f.split("/")[-1].split("_")[-2]
            df = pd.read_csv(f, sep="\t")
            df.insert(df.shape[1], "Filetype", filetype)
            dfs.append(df)

        df = pd.concat(dfs)
        df[df.Function == wildcards.function].to_csv(output[0], sep="\t")


def fix_description(desc):

    from textwrap import wrap
    assert len(desc) < 160, "Description too long (>= 160 chars.)"
    return "\n".join(wrap(desc, width=80)) + "\n"

descriptions = pd.read_csv("supplementaries/descriptions.yaml", sep="\t", header=0)
print(descriptions)

rule graph_time_memory_together:
    input:
        "{prefix}/benchmark/collected_timings_mean_{function}.txt"
    output:
        "{prefix}/benchmark/graphs/time_memory_together_{function}.{extension}"
    params:
        function = lambda w: w.function.capitalize(),
        description = lambda w: fix_description(descriptions[descriptions.Function == w.function].Description.iloc[0])
    script:
        "../scripts/graph_time_mem_together.R"
