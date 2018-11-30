
def get_files(w):
    # print("category", w.category)


    return [f for f in category_dict[w.category]]


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

            if library != "bioconductor":
                minutes, seconds, fraction = timing.split(".")
                minutes, seconds = int(minutes), int(seconds)
                seconds += minutes * 60

                timing = ".".join(str(s) for s in [seconds, fraction])

            timing = np.log10(float(timing))
            try:
                max_rss = pd.read_table(bmark_f, sep="\t", usecols=[2], skiprows=1, squeeze=True, header=None).values[0] / 1024

                rowdict = {"Iteration": iteration, "MaxRSSGB": max_rss,
                           "Seconds": timing, "Function": function, "Library": library, "Log10NBIntervals":
                           np.log10(size)}
            # because maxrssgb does not work on macOS
            except:
                rowdict = {"Iteration": iteration, "Seconds": timing, "Function": function, "Library": library, "Log10NBIntervals":
                           np.log10(size)}


            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts)
        print(df.head())
        df = df.sort_values("Function Library Log10NBIntervals".split())

        try:
            column_order = "Function Library Log10NBIntervals MaxRSSGB Seconds Iteration".split()
            df[column_order].to_csv(output[0], index=False, sep="\t")
        except:
            column_order = "Function Library Log10NBIntervals Seconds Iteration".split()
            df[column_order].to_csv(output[0], index=False, sep="\t")


rule find_mean_times:
    input:
        "{prefix}/benchmark/collected_timings_{filetype}_{category}.txt"
    output:
        "{prefix}/benchmark/collected_timings_mean_{filetype}_{category}.txt"
    run:
        df = pd.read_table(input[0], sep="\t", header=0)

        rowdicts = []
        for g, gdf in df.groupby("Function Library Log10NBIntervals".split()):

            seconds = gdf.Seconds.mean()
            seconds_sd = gdf.Seconds.std()
            memory = gdf.MaxRSSGB.mean()
            print(g)
            print(gdf)
            print(gdf.MaxRSSGB)
            memory_sd = gdf.MaxRSSGB.std()
            print(memory_sd)
            print("-------")

            head = gdf.head(1)
            rowdict = {}
            for row in gdf:
                rowdict[row] = gdf[row].iloc[0]

            rowdict.update({"Seconds": seconds, "MaxRSSGB": memory, "SecondsSD": seconds_sd, "MemorySD": memory_sd})

            rowdicts.append(rowdict)

        df = pd.DataFrame.from_dict(rowdicts)["Function Library Log10NBIntervals MaxRSSGB Seconds SecondsSD MemorySD".split()]
        print(df.head())
        df.to_csv(output[0], sep="\t", index=False)
