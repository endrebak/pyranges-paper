from helpers import read_file, file_to_grange, init_ray; init_ray(snakemake.wildcards)

from time import time
import datetime

operation = snakemake.wildcards.unary_operation
w = snakemake.wildcards

f = snakemake.input[0]

if operation != "dataframe_to_genomicrange":
    gr = file_to_grange(f, filetype=snakemake.wildcards.filetype)
else:
    from pyranges import PyRanges
    df = read_file(f, filetype=snakemake.wildcards.filetype)

print(snakemake.params.code)
start = time()
exec(snakemake.params.code)
end = time()
total = end - start

print(result)

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output.time, "w+").write(minutes_seconds)

result_string = str(result)

open(snakemake.output.result, "w+").write(result_string)

if int(w.size) == int(1e5) and int(w.iteration) == 0:
    # print("writing df! " * 10)
    path = "{prefix}/benchmark/actual_results/{function}/{filetype}/pyranges.txt".format(prefix=w.prefix, function=operation, filetype=w.filetype)
    import subprocess, os
    subprocess.check_output("mkdir -p " + os.path.dirname(path), shell=True)
    # print(path)
    try:
        result.df.sort_values("Chromosome Start End".split()).to_csv(path, sep="\t", index=False)
    except:
        pass
