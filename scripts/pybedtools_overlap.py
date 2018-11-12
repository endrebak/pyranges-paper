from pybedtools import BedTool

from time import time
import datetime

c = snakemake.input.chip
b = snakemake.input.background

chip = BedTool(c)
bg = BedTool(b)

start = time()

# is_sorted = True if snakemake.wildcards.sorted == "sorted" else False

# if is_sorted:
#     result = chip.intersect(bg, s=True, sorted=is_sorted, wa=True)
# else:
result = chip.intersect(bg, s=True, wa=True)

end = time()

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)
