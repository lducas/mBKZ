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
from g6k.utils.stats import dummy_tracer
from g6k.algorithms.workout import workout
from g6k import Siever
from fpylll import IntegerMatrix, FPLLL, BKZ, GSO, LLL
from fpylll import Enumeration, EvaluatorStrategy
from fpylll.fplll.enumeration import EnumerationError


def enumerate_all_short_vectors(self, sq_rad):
    """
    Output a list of (almost) all vectors of squared euclidean
    norm less than sq_rad.
    """
    n = self.M.d
    workout(self, dummy_tracer, 0, n, 0)

    nr_solutions = 100
    res = []
    while nr_solutions < 101 or len(res) >= nr_solutions:
        nr_solutions *= 2
        enum = Enumeration(self.M, 
            strategy=EvaluatorStrategy.FIRST_N_SOLUTIONS, 
            nr_solutions=nr_solutions)
        try:
            res = enum.enumerate(0, n, sq_rad, 0)
        except EnumerationError:
            pass
    svs = [self.M.B.multiply_left(x) for (_, x) in res]
    return svs



c = 16
q = 16317
ft = "ld"
r = 2

K = Cyclotomic(c)
B = random_qary_cyclotomic(K, r, r//2, q)
M = Siever.MatGSO(IntegerMatrix.from_matrix(B))
g6k = Siever(M)
n = g6k.M.d
print("Created a lattice of dimension %d"%n)
GH = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])
f = 2
print("Enumerating up to %.4f * GH"%f)

svs = enumerate_all_short_vectors(g6k, GH * f**2)
print("Found %d short vectors"%len(svs))

data = []

for x in svs:
	K.check_cyclic_embedding(x)
	eucl_norm = K.ip_Q(x, x)
	alg_norm = float(K.algebraic_norm(K.ip_K(x, x))**.5)
	data.append((eucl_norm, alg_norm))

print("Shortest vector for the Euclidean Norm")
data.sort(key=lambda x : x[0])
print(data[0])

print("Shortest vector for the Algebraic Norm")
data.sort(key=lambda x : x[1])
print(data[0])
