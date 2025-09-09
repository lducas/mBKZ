#!/usr/bin/env python3
import sys
sys.path.insert(1, './g6k')

from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix, FPLLL, BKZ, GSO, LLL
from time import time
from math import ceil
from numpy import array
from cyclotomics import Cyclotomic
from modlatred import random_qary_cyclotomic, ModuleLatticeReduction, slope


ld = "ld"
dd = "dd"
m = 160
c = 3
n = None
q = 16317
ft = "ld"
struct=True

for s in sys.argv[1:]:
    exec(s)

K = Cyclotomic(c)
n = n or K.deg * (m//(2*(K.deg)))
assert(m % K.deg == 0)
assert(n % K.deg == 0)

print(f"m={m}, n={n}, q={q}, cond={K.cond}, deg={K.deg}, ft={ft}", file=sys.stderr)

print(f"beta,\ttours,\tslopeQ, \tslopeK,  \ttime,\ttotal-time")

T0 = time()
Tlast = time()

def barf(beta, rcount=0):
    global Tlast, T0, mlr
    T = time()
    sQ = slope(mlr.profile_Q())
    sK = slope(mlr.profile_K())
    print(f"{beta:3d},\t{rcount:3d},\t{sQ:.5f},\t{sK:.5f},\t{T-Tlast:4.2f},\t{T-T0:4.2f}")
    Tlast = T

trials = 0
while True:
    trials += 1
    try:
        B = random_qary_cyclotomic(K, m//K.deg, n//K.deg, q)
        if struct:
            mlr = ModuleLatticeReduction(B, K, float_type=ft, restructure_delta_prog=.03)
        else:
            Q = Cyclotomic(1)
            mlr = ModuleLatticeReduction(B, Q, float_type=ft)
        for t in range(5):
            mlr.lll(0, m, delta=.99)
            mlr.restructure()
        barf(1, trials)
        barf(2, 5)
        break
    except:
        print("failed to restructure, trying a new lattice", file=sys.stderr)
        print(f"{sys.argv[0]}: m={m}, n={n}, q={q}, cond={K.cond}, deg={K.deg}, ft={ft}", file=sys.stderr)


tours = 5*mlr.K.deg
for beta in range(ceil(5/mlr.K.deg)*mlr.K.deg, 81, mlr.K.deg):
    mlr.bkz(beta, tours)
    barf(beta, tours)
