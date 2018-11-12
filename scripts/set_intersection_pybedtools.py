
from pybedtools import BedTool

from time import time
import datetime

c = snakemake.input.chip
b = snakemake.input.background

chip = BedTool(c)
bg = BedTool(b)

start = time()

# is_sorted = True if snakemake.wildcards.sorted == "sorted" else False

# if not is_sorted:
#     sc = chip.sort(stream=True).merge(s=True, c=[4, 5, 6], o="distinct")
#     sb = bg.sort(stream=True).merge(s=True, c=[4, 5, 6], o="distinct")
# else:
sc = chip.merge(s=True, c=[4, 5, 6], o="distinct")
sb = bg.merge(s=True, c=[4, 5, 6], o="distinct")

result = sc.intersect(sb, s=True)

end = time()

# print(result)

total = end - start

total_dt = datetime.datetime.fromtimestamp(total)

minutes_seconds = total_dt.strftime('%-M.%-S.%f')

open(snakemake.output[0], "w+").write(minutes_seconds)

# result.df.to_csv(snakemake.output.result, sep=" ", index=False)
