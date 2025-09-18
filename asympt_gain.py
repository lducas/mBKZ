from cyclotomics import Cyclotomic
from math import log, pi
import csv

Q = Cyclotomic(1)
cs = [1, 3, 4, 5, 7, 8, 9, 11, 13, 15, 16, 17, 32]
Ks = [Cyclotomic(c) for c in cs]


def asympt_beta_eq(beta, K):
    d = K.deg
    return beta + (d - 1 if abs(K.disc) == d**d else (log(abs(K.disc)) - d*log(d))*beta/(d*(log(beta))))

writer = csv.writer(open("data/asympt_gain.csv", 'w'))
writer.writerow(['beta', 'beta_eq1' , 'beta_eq3' , 'beta_eq4' , 'beta_eq5' , 'beta_eq7' , 'beta_eq8' , 'beta_eq9' , 'beta_eq11', 'beta_eq13', 'beta_eq15', 'beta_eq16', 'beta_eq17', 'beta_eq32'])

print("beta, " + ", ".join([f"beta_eq{c}"+(" " if c<10 else "") for c in cs]))
for beta in range(60, 1001, 10):
    print(f"{beta:4d}, " + ", ".join([f"{asympt_beta_eq(beta, K) - beta:9.3f}" for K in Ks]))
    writer.writerow([beta]+[asympt_beta_eq(beta, K) - beta for K in Ks])
