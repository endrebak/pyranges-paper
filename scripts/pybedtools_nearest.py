
from pybedtools import BedTool

from time import time
import datetime

c = snakemake.input.chip
b = snakemake.input.background

chip = BedTool(c)
bg = BedTool(b)

start = time()

result = chip.sort().closest(bg.sort(), s=True)

end = time()

# print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

# result.df.to_csv(snakemake.output.result, sep=" ", index=False)
