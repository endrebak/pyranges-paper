import pandas as pd
from pyranges import PyRanges
from pyrle import PyRles

from sys import argv

infile = argv[1]
infile_bg = argv[2]
runs = int(argv[3])

chip_f = infile
background_f = infile_bg

nrows = None
chip = pd.read_table(chip_f, sep="\t", nrows=nrows,
                     usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())

grles_chip = PyRles(chip, stranded=True)

background = pd.read_table(background_f, sep="\t", nrows=nrows,
                           usecols=[0, 1, 2, 5], header=None, names="Chromosome Start End Strand".split())
grles_background = PyRles(background, stranded=True)

for i in range(runs):

    grles_chip.div(grles_background)
