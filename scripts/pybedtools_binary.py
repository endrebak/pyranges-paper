from pybedtools import BedTool
import subprocess
import tempfile

from time import time
import datetime

c = snakemake.input.chip
b = snakemake.input.background

w = snakemake.wildcards
if w.filetype == "annotation":
    extension = "gtf"
else:
    extension = "bed"

pb1 = BedTool(c)
pb2 = BedTool(b)

op = snakemake.params.code
print(op)

start = time()

exec(op)

end = time()

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output.time, "w+").write(minutes_seconds)

with tempfile.NamedTemporaryFile() as f:

    try:
        result.saveas(f.name)
    except:
        open(f.name, "w+").write(str(result))

    cmd = "head {} > {}".format(f.name, snakemake.output.result)
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = """wc -l {} | py -x '"Number of lines: " + x.strip().split()[0]' >> {}""".format(f.name, snakemake.output.result)
    subprocess.check_output(cmd, shell=True)


if int(w.size) == int(1e5) and int(w.iteration) == 0:
    path = "{prefix}/benchmark/actual_results/{function}/{filetype}/pybedtools.txt".format(prefix=w.prefix, function=w.operation, filetype=w.filetype)
    # print(path)
    import subprocess, os
    result.head()
    subprocess.check_output("mkdir -p " + os.path.dirname(path), shell=True)
    try:
        print("saving to path")
        result.saveas(path)
    except:
        print("an error occured in " + w.operation)
        pass
    # print(path)
    # subprocess("sort -k1,1 -k2,3n {}")
