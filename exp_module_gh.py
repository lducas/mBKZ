#!/usr/bin/env python3
import sys
sys.path.insert(1, './g6k')

from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix, FPLLL, BKZ, GSO, LLL
from time import time
from math import ceil
from numpy import array
from numpy.linalg import slogdet
from cyclotomics import Cyclotomic
from modlatred import random_qary_cyclotomic, ModuleLatticeReduction, slope
from predictions import lghK, lghZ

c = 3
q = 16317
ft = "ld"
samples = 1000
maxdim = 80
mindim = 5

for s in sys.argv[1:]:
    exec(s)

K = Cyclotomic(c)

print(f"q={q}, cond={K.cond}, deg={K.deg}, ft={ft}, maxdim={maxdim}, samples={samples}", file=sys.stderr)

print("   r,   rd, expe_gap, pred_gap")
for r in range(ceil(mindim/K.deg), maxdim//K.deg):
    l_lambda1 = []
    while len(l_lambda1) < samples:
        B = random_qary_cyclotomic(K, r, r//2, q)
        try:
            mlr = ModuleLatticeReduction(B, K, restructure_delta_prog=.03 if K.cond<16 else .01)
            mlr.svp_reduce()
        except:
            print("Failed to restructure, trying a new lattice", file=sys.stderr)
            print(f"{sys.argv[0]}: q={q}, cond={K.cond}, deg={K.deg}, ft={ft}, r={r}, samples={samples}", file=sys.stderr)
            continue

        prof = mlr.profile_Q()
        l_lambda1.append((prof[0] - sum(prof)/len(prof)))

    avg_l_lambda1 = sum(l_lambda1)/samples
    pred_lghZ = lghZ(r * K.deg)
    pred_lghK = lghK(K, r)
    print(f"{r:3d}, {r*K.deg:3d},{avg_l_lambda1 - pred_lghZ:.5f},{pred_lghK - pred_lghZ:.5f}")


