from cyclotomics import Cyclotomic
from predictions import slope_Q
import sys


Q = Cyclotomic(1)
cs = [3, 4, 5, 7, 8, 9, 15, 16]
Ks = [Cyclotomic(c) for c in cs]

def sign(num):
    return -1 if num < 0 else 1

ub = "upper bound"
le = "lower estimate"
estimate = ub

for s in sys.argv[1:]:
    exec(s)


def find_beta_eq(beta, K):
    # Search beta equivalent
    if beta < 60:
        raise ValueError("Not confident that predicted slope is decreasing below that point")
    sl = slope_Q(Q, beta)
    betaK = beta//K.deg
    slK = slope_Q(K, betaK)
    sign = 1 if slK < sl else -1
    while sign * slK < sign * sl:
        betaK += sign
        slK = slope_Q(K, betaK, estimate=estimate)

    # Interpolate to non-integral value for smooth plots
    slK_ = slope_Q(K, betaK - sign, estimate=estimate)
    betaKeq = betaK - sign * (sl - slK)/(slK_ - slK)
    return K.deg * betaKeq


print("beta, " + ", ".join([(" " if c<10 else "") + f"beta_eq{c}" for c in cs]))
for beta in range(70, 1001, 10):
    print(f"{beta:4d}, " + ", ".join([f"{find_beta_eq(beta, K) - beta:9.3f}" for K in Ks]))

