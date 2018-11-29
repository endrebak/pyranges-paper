
from pybedtools import BedTool

from time import time
import datetime

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

open(snakemake.output[0], "w+").write(minutes_seconds)
