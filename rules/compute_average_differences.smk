import pandas as pd

rule find_average_differences:
    input:
        "{prefix}/benchmark/collected_timings_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/differences/{filetype}_{category}_differences.txt"
    run:
        f = input[0]
        df = pd.read_csv(f, sep="\t", header=0)

        rowdicts = []
        for g, gdf in df.groupby("Function"):
            if "bioconductor" not in list(gdf.Library):
                continue

            gdf = gdf[(gdf.Library == "bioconductor") | (gdf.Library == "pyranges_1")]
            # gdf = gdf[gdf.Log10NBIntervals == float(wildcards.size)]
            for size, gdf2 in gdf.groupby("Log10NBIntervals"):

                gl = gdf2.groupby("Library")
                mean_seconds = gl.Seconds.mean()
                mean_mem = gl.MaxRSSGB.mean()

                sdiff = mean_seconds["bioconductor"] - mean_seconds["pyranges_1"]
                memdiff = mean_mem["bioconductor"] - mean_mem["pyranges_1"]

                rowdicts.append({"Mem": memdiff, "Seconds": sdiff, "Function": g, "Size": size})

        c = pd.DataFrame.from_dict(rowdicts)

        c.to_csv(output[0], index=False, sep="\t")
