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
from fpylll.util import gaussian_heuristic


c = 8
q = 16317
ft = "ld"
r = 2

K = Cyclotomic(c)
B = random_qary_cyclotomic(K, r, r//2, q)
rZ = r * K.deg
mlr = ModuleLatticeReduction(B, K)

GH = gaussian_heuristic([mlr.M.get_r(i, i) for i in range(mlr.M.d)])

svs = mlr.enumerate_all_short_vectors(GH * 4)
print(len(svs))

data = []

for x in svs:
	eucl_norm = K.ip_Q(x, x)
	alg_norm = float(K.algebraic_norm(K.ip_K(x, x))**.5)
	data.append((eucl_norm, alg_norm))

print("Shortest vector for the Euclidean Norm")
data.sort(key=lambda x : x[0])
print(data[0])

print("Shortest vector for the Algebraic Norm")
data.sort(key=lambda x : x[1])
print(data[0])
