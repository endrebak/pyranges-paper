import pandas as pd
from pyranges import PyRanges

from sys import argv

infile = argv[1]
runs = int(argv[2])

chip = pd.read_table(infile, sep="\t",
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

for i in range(runs):
    print("runs, iteration:", runs, i)
    granges = PyRanges(chip)
