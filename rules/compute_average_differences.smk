import pandas as pd

rule find_average_differences:
    input:
        "{prefix}/benchmark/collected_timings_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/differences/{num_cores}_{filetype}_{category}_differences.txt"
    run:
        f = input[0]
        df = pd.read_csv(f, sep="\t", header=0)

        cores = str(wildcards.num_cores)

        rowdicts = []
        for g, gdf in df.groupby("Function"):
            if "bioconductor" not in list(gdf.Library):
                continue

            gdf = gdf[(gdf.Library == "bioconductor") | (gdf.Library == "pyranges_" + cores)]

            for size, gdf2 in gdf.groupby("Log10NBIntervals"):

                gl = gdf2.groupby("Library")
                mean_seconds = gl.Seconds.mean()
                mean_mem = gl.MaxRSSGB.mean()

                sdiff = 10 ** mean_seconds["bioconductor"] / (10 ** mean_seconds["pyranges_" + cores])
                memdiff = mean_mem["bioconductor"] / mean_mem["pyranges_" + cores]

                rowdicts.append({"Mem": memdiff, "Seconds": sdiff, "Function": g, "Size": size})

        c = pd.DataFrame.from_dict(rowdicts)

        c.to_csv(output[0], index=False, sep="\t")


# rule mean_average_difference:
#     input:
#         expand("{{prefix}}/benchmark/differences/{{num_cores}}_{filetype}_{{category}}_differences.txt", filetype="annotation reads".split())
#     output:
#         "{prefix}/benchmark/differences/{num_cores}_{category}_mean_differences.txt"
#     run:
#         dfs = []
#         for f in input:
#             dfs.append(pd.read_csv(f, sep="\t", header=0))

#         df = pd.concat(dfs)

#         rowdicts = []
#         for g, gdf in df.groupby("Size"):
#             mem = gdf.Mem.mean()
#             seconds = gdf.Seconds.mean()

#             rowdicts.append({"Seconds": seconds, "Mem": mem, "Size": g})

#         odf = pd.DataFrame.from_dict(rowdicts)
#         print(odf)
