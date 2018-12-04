
import tempfile
import subprocess
from pybedtools import BedTool

from time import time
import datetime


w = snakemake.wildcards
if w.filetype == "annotation":
    extension = "gtf"
else:
    extension = "bed"

c = snakemake.input[0]

pb1 = BedTool(c)

op = snakemake.params.operation
print(op)

start = time()

exec(op)

end = time()

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output.time, "w+").write(minutes_seconds)


with tempfile.NamedTemporaryFile() as f:

    result.saveas(f.name)

    cmd = "head {} > {}".format(f.name, snakemake.output.result)
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    cmd = """wc -l {} | py -x '"Number of lines: " + x.strip().split()[0]' >> {}""".format(f.name, snakemake.output.result)
    subprocess.check_output(cmd, shell=True)

if int(w.size) == int(1e5) and int(w.iteration) == 0:
    path = "{prefix}/benchmark/actual_results/{function}/{filetype}/pybedtools.txt".format(prefix=w.prefix, function=w.unary_operation, filetype=w.filetype)
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
