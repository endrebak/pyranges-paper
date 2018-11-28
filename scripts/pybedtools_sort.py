
from pybedtools import BedTool

from time import time
import datetime

c = snakemake.input[0]

chip = BedTool(c)

start = time()

result = chip.sort()

end = time()

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
