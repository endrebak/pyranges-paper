from pandas import read_table

from sys import argv

f1 = argv[1]
f2 = argv[2]

nrows = None
chip = read_table(f1, nrows=nrows)
bg = read_table(f2, nrows=nrows)
