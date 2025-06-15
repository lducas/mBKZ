import sys
sys.path.insert(1, '../../')
import matplotlib.pyplot as plt
import numpy as np
from read_csv import read_csv
from math import log, gamma, pi
from predictions import slope_Q, slope_K
from cyclotomics import Cyclotomic

seeds = 5

cs = [1, 3, 4, 5, 8, 15, 16]

for c in cs:
	K = Cyclotomic(c)
	d = K.deg

	fins = ["slope_c%d_seed%d.csv"%(c, s) for s in range(seeds)]
	Ds = [read_csv(f) for f in fins]
	fout = open("slope_c%d_avg.csv"%c, 'w')
	l = min([len(D["beta"]) for D in Ds])
	print("beta, betaK, slope_Q, pred_Q_low, pred_Q_high", file=fout)
	for i in range(l):
		avgQ = sum([D["slopeQ"][i] for D in Ds])/seeds		
		# avgQ /= 2 # Because the old version used to report slope of squared norms
		beta = Ds[0]["beta"][i]
		if beta < 10:
			continue
		assert(not beta % d)
		betaK = beta//d
		predQ_low = slope_Q(K, betaK, "lower estimate")
		predQ_high = slope_Q(K, betaK, "upper bound")

		print(f"{beta:3}, {betaK:3}, {avgQ:.6f}, {predQ_low:.6f}, {predQ_high:.6f}", file=fout)
	fout.close()
