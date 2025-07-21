#!/usr/bin/env python3
import sys
sys.path.insert(1, './g6k')

from fpylll.util import gaussian_heuristic
from fpylll import IntegerMatrix, FPLLL, BKZ, GSO, LLL
from time import time
from math import ceil, sqrt
from numpy import array, apply_along_axis, reshape, concatenate, block, zeros
from cyclotomics import Cyclotomic
from modlatred import random_qary_cyclotomic, ModuleLatticeReduction, slope, Z


ld = "ld"
dd = "dd"
c = 3
n = 20
q = 16317
ft = "ld"
mBKZ = False
mKannan = False
verb = True

for s in sys.argv[1:]:
    exec(s)

if mBKZ:
    assert(mKannan)

K = Cyclotomic(c)
m = 2*n
assert(m % K.deg == 0)
assert(n % K.deg == 0)
d = K.deg

if verb:
    print(f"m={m}, n={n}, q={q}, cond={K.cond}, deg={K.deg}, ft={ft}, mBKZ={mBKZ}, , mBKZ={mKannan}", file=sys.stderr)
T0 = time()
Tlast = time()

tot = 0
sigma = 2
s = K.spherical_sample(n//d, sigma=sigma)
e = K.spherical_sample(n//d, sigma=sigma)

B = random_qary_cyclotomic(K, m//d, n//d, q)
A = B[n:, :c*n//d]

s_ = reshape(s, (len(s)//c, c))
s_ = apply_along_axis(K.reduce_mod_defpoly, 1, s_)//c
s_ = s_.flatten()

t = (s_ @ A + e)
t = reshape(t, (len(t)//c, c))
t = apply_along_axis(K.reduce_mod_defpoly, 1, t)//c
t = t % q
t = apply_along_axis(K.cyclic_embedding, 1, t)
t = t.flatten()



t0 = concatenate((t, zeros(c*n//d, dtype=int)))
B_ = block([[B, Z(2*n, c)], [K.vOK_Zbasis(t0), K.vOK_Zbasis(K.one)]])

if not mKannan:
    B_ = B[:m+1]

# print(B_)
norm2 = (s@s) + (e@e)

if verb:
    print("e|s :", e, s)
    print("norm^2:", norm2)

if mBKZ:
    mlr = ModuleLatticeReduction(B_, K, float_type=ft, restructure_delta_prog=.03)
else:
    Q = Cyclotomic(1)
    mlr = ModuleLatticeReduction(B_, Q, float_type=ft, restructure_delta_prog=.03)

for t in range(5):
    mlr.lll(0, mlr.M.d)
    mlr.restructure()

ad = d if mBKZ else 1

tours = 5*ad
for beta in range(ceil(3/ad)*ad, 81, ad):
    v = array(mlr.M.B[0])[:c*m//d]
    if (v @ v) == norm2:
        break
    print(f"Running mBKZ_{mlr.K.cond} beta={beta} (mBKZ={mBKZ}, mBKZ={mKannan})", file=sys.stderr)
    mlr.bkz(beta, tours)
    mlr.lll(0, mlr.M.d)
    mlr.restructure()

#print("done")
print(beta - ad)
# print("stopped before beta=", beta)
if verb:
    print("found norm", v @ v)
    print("vector", v)
