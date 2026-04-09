#!/usr/bin/env python3
import sys
sys.path.insert(1, './g6k')

from cyclotomics import Cyclotomic
from math import log, exp, ceil
from modlatred import random_qary_cyclotomic, ModuleLatticeReduction
from g6k.algorithms.workout import workout
from g6k.utils.stats import dummy_tracer
import scipy.special as sc
from numpy.linalg import slogdet
from numpy import array
from predictions import avg_lskewness, avg_lindex


c = 3
q = 16317
ft = "ld"
samples = 1000
maxdim = 80
mindim = 5

for s in sys.argv[1:]:
	exec(s)


def lskewness(K, v):
	# 3rd term of equation (3)
	eucl_norm = K.ip_Q(v, v)**.5
	alg_norm = K.algebraic_norm(K.ip_K(v, v))**.5
	log_gain = log(K.deg)/2 + log(alg_norm)/K.deg - log(eucl_norm)
	return log_gain


def experiment_skew_index(K, r, samples, q=19997):
	L_lskew, L_lindex, L_nzindex = [], [], []

	if r < 2:
		raise ValueError("Sampling rank-1 modules not implemented \
						  as there is no clear way of doing so...")

	while len(L_lskew) < samples:
		B = random_qary_cyclotomic(K, r, r//2, q)
		rZ = r * K.deg
		try:
			mlr = ModuleLatticeReduction(B, K, restructure_delta_prog=.03 if K.cond<16 else .01)
			mlr.svp_reduce()
		except:
			print("failed to restructure, trying a new lattice", file=sys.stderr)
			print(f"{sys.argv[0]}: q={q}, cond={K.cond}, deg={K.deg}, ft={ft}, r={r}, samples={samples}", file=sys.stderr)
			continue

		v = mlr.M.B[0]
		vOK = K.vOK_Zbasis(v)
		I = array(list(mlr.M.B[0:K.deg]))
		_, vOK_ldet = slogdet(vOK @ vOK.transpose())
		_, I_ldet = slogdet(I @ I.transpose())

		l_lindex = (vOK_ldet - I_ldet)/2
		index = int(round(exp(l_lindex)))
		assert(abs(index/exp(l_lindex) - 1) < 1e-5)
		assert(index >= 1)

		L_lskew.append(lskewness(K, v))
		L_lindex.append(-log(index)/K.deg)
			
	return sum(L_lskew)/samples, sum(L_lindex)/samples

K = Cyclotomic(c)

print("   r,   rd, expe_skew_gap, pred_skew_gap, expe_index_gap, pred_index_gap")
for r in range(max(2, ceil(mindim/K.deg)), maxdim//K.deg):
    B = avg_lskewness(K, r)
    D = avg_lindex(K, r)
    A, C = experiment_skew_index(K, r, samples)
    print(f"{r:3d}, {r*K.deg:3d},{A:.5f}, {B:.5f}, {C:.5f}, {D:.5f}")

