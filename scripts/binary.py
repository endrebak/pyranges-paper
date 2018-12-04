from helpers import read_file, file_to_grange, init_ray; init_ray(snakemake.wildcards)

from time import time
import datetime

w = snakemake.wildcards

operation = w.operation
code = snakemake.params.code

print(code)

f = snakemake.input.chip
f2 = snakemake.input.background
gr = file_to_grange(f)
gr2 = file_to_grange(f2, filetype=snakemake.wildcards.filetype)

start = time()
exec(snakemake.params.code)
end = time()
total = end - start

print(result)

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

result_string = str(result)

open(snakemake.output.result, "w+").write(result_string)

# print(w)
if int(w.size) == int(1e5) and int(w.iteration) == 0:
    # print("writing df! " * 10)
    path = "{prefix}/benchmark/actual_results/{function}/{filetype}/pyranges.txt".format(prefix=w.prefix, function=w.operation, filetype=w.filetype)
    import subprocess, os
    subprocess.check_output("mkdir -p " + os.path.dirname(path), shell=True)
    # print(path)
    try:
        result.df.sort_values("Chromosome Start End".split()).to_csv(path, sep="\t", index=False)
    except:
        pass
