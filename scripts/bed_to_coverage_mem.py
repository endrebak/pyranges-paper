import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles
from sys import argv

infile = argv[1]
runs = int(argv[2])

chip = pd.read_table(infile, sep="\t",
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

for i in range(runs):

    print("bed_to_coverage_mem: runs, iteration:", runs, i)

    grles = PyRles(chip, stranded=True)
