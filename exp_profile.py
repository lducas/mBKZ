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
m = 60
c = 3
n = None
q = 16317
ft = "ld"
beta = 20
tours = 30

for s in sys.argv[1:]:
    exec(s)

K = Cyclotomic(c)
n = n or K.deg * (m//(2*(K.deg)))
assert(m % K.deg == 0)
assert(n % K.deg == 0)
assert(beta % K.deg == 0)

print(f"m={m}, n={n}, q={q}, cond={K.cond}, deg={K.deg}, ft={ft}, beta={beta}, tours={tours}", file=sys.stderr)
T0 = time()
Tlast = time()

trials = 0
while True:
    trials += 1
    try:
        B = random_qary_cyclotomic(K, m//K.deg, n//K.deg, q)
        mlr = ModuleLatticeReduction(B, K, float_type=ft, restructure_delta_prog=.03)
        for t in range(5):
            mlr.lll(0, m)
            print("pre LLL tour=%d (cond=%d)"%(t, K.cond), file=sys.stderr)
            mlr.restructure()
        break
    except:
        print("failed to restructure, trying a new lattice", file=sys.stderr)
        print(f"{sys.argv[0]}: m={m}, n={n}, q={q}, cond={K.cond}, deg={K.deg}, ft={ft}", file=sys.stderr)


for b in range(0, beta, K.deg):
  if b<10:
    continue
  print("pre progBKZ-%d (cond=%d)"%(b, K.cond), file=sys.stderr)
  mlr.bkz(b, 1)

for t in range(tours):
    print("BKZ-%d tour=%d (cond=%d)"%(beta, t, K.cond), file=sys.stderr)
    mlr.bkz(beta, 1)
    mlr.restructure()

lQ = mlr.profile_Q(True)
lK = mlr.profile_K(True)
d = K.deg

print("  i, ellQ, ellK ")
for i in range(m):
    lKi = f"{lK[i//d]/d:.5f}" if not i%d else ""
    print(f"{i:3d}, {lQ[i]:.5f}, "+lKi)
