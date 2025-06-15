import sys
sys.path.insert(1, '../../')
import matplotlib.pyplot as plt
import numpy as np
from read_csv import read_csv
from math import log, gamma, pi
from cyclotomics import Cyclotomic

seeds = 5

cs = [1, 3, 4, 5, 8, 15, 16]

for c in cs:
	print(c)
	fins = ["prof_c%d_seed%d.csv"%(c, s) for s in range(seeds)]
	Ds = [read_csv(f) for f in fins]
	fout = open("prof_c%d_avg.csv"%c, 'w')
	l = min([len(D["i"]) for D in Ds])
	print("  i, ellQ, ellK", file=fout)
	for i in range(l):
		try:
			avgK = f"{sum([D["ellK"][i] for D in Ds])/seeds:.6f}"
		except:
			avgK = ""
		avgQ = f"{sum([D["ellQ"][i] for D in Ds])/seeds:.6f}"
		print(f"{i:3}, {avgQ}, {avgK}", file=fout)
	fout.close()