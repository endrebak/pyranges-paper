
import pandas as pd
from pyranges import PyRanges

from sys import argv

infile = argv[1]
infile_bg = argv[2]
runs = int(argv[3])

chip_f = infile
background_f = infile_bg

nrows = None
chip = pd.read_table(chip_f, sep="\t", nrows=nrows,
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

background = pd.read_table(background_f, sep="\t", nrows=nrows,
                           usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
cgr = PyRanges(chip)
bgr = PyRanges(background)

for i in range(runs):

    print("nearest", i, runs)

    result = cgr.nearest(bgr, strandedness="same")
