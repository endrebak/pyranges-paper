from pybedtools import BedTool
import subprocess
import tempfile

from time import time
import datetime

c = snakemake.input.chip
b = snakemake.input.background

# print("c =", c)
# print("b =", b)

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
