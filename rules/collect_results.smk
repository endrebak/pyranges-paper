
def get_files(w):

    return [f for f in category_dict[w.category] if "/" + w.filetype + "/" in f]


rule collect_times:
    input:
        get_files
    output:
        "{prefix}/benchmark/collected_timings_{filetype}_{category}.txt"
    run:
        rowdicts = []
        for f in input:
            # print(f)
            bmark_f = f.replace("time.txt", "benchmark.txt")

            function, library, filetype, timingfile = f.split("/")[-4:]
            iteration, size = timingfile.split("_")[:2]
            size = int(size)
            # print(function, library, timingfile )
            # print(iteration, size)

            timing = open(f).readlines()[0].strip()

            if "bioconductor" not in library:
                minutes, seconds, fraction = timing.split(".")
                minutes, seconds = int(minutes), int(seconds)
                seconds += minutes * 60

                timing = ".".join(str(s) for s in [seconds, fraction])

            timing = np.log10(float(timing))

            max_rss = pd.read_csv(bmark_f, sep="\t", usecols=[2], skiprows=1, squeeze=True, header=None).values[0] / 1024

            rowdict = {"Iteration": iteration, "MaxRSSGB": max_rss,
                        "Seconds": timing, "Function": function, "Library": library, "Log10NBIntervals":
                        np.log10(size)}

            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts)
        print(df.head())
        df = df.sort_values("Function Library Log10NBIntervals".split())

        column_order = "Function Library Log10NBIntervals MaxRSSGB Seconds Iteration".split()
        df[column_order].to_csv(output[0], index=False, sep="\t")


rule find_mean_times:
    input:
        "{prefix}/benchmark/collected_timings_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    run:
        df = pd.read_csv(input[0], sep="\t", header=0)

        rowdicts = []
        for g, gdf in df.groupby("Function Library Log10NBIntervals".split()):

            seconds = gdf.Seconds.mean()
            seconds_sd = gdf.Seconds.std()
            memory = gdf.MaxRSSGB.mean()
            memory_sd = gdf.MaxRSSGB.std()

            head = gdf.head(1)
            rowdict = {}
            for row in gdf:
                rowdict[row] = gdf[row].iloc[0]

            rowdict.update({"Seconds": seconds, "MaxRSSGB": memory, "SecondsSD": seconds_sd, "MemorySD": memory_sd})

            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts)["Function Library Log10NBIntervals MaxRSSGB Seconds SecondsSD MemorySD".split()]
        print(df.head())
        df.to_csv(output[0], sep="\t", index=False)
