import matplotlib.pyplot as plt
import numpy as np
from math import log, lgamma, pi
from cyclotomics import Cyclotomic
from scipy.special import digamma

def lVball(n): 
    """
    Logarithmic volume of a Euclidean n-ball
    :param n: dimension n
    "returns: log(Voll(n-ball))
    """    
    return (n/2)*log(pi) - lgamma(n/2+1)

def lghZ(n):
    """
    Logarithmic Gaussian Heuristic for Z-lattices (Section 3.1)
    Formula from Theorem 1 of https://eprint.iacr.org/2018/856
    :param n: dimension n of the lattice
    "returns: lgh_Q(n)
    """    
    return log(2)/n + lgamma(1+1/n) - lVball(n)/n

def lghK(K, r):
    """
    Logarithmic Gaussian Heuristic for module lattices over a cyclotomic field K (Section 4.2)
    Based on Theorem 3 of https://arxiv.org/abs/2402.10305, with adjustment discussed
    below Heuristic 4
    :param K: cyclotomic number field K
    :param r: K-rank r of the module lattice
    "returns: lgh_K(n)
    """    
    return lghZ(K.deg*r) + log(K.nb_rootsofunity/2)/(K.deg*r)
    # This version below, based on Theorem 3 appears to give less precise results
    # for even conductors.
    #return lghZ(K.deg*r) + log(K.cond)/(K.deg*r) 


def avg_lskewness(K, r, estimate="spherical model"):
    """
    Estimate the average log skewness term in the slope equation (Term t_3)
    Namely, log(sqrt(d) N(b)^(1/d) / ||b||) following Section 4.4
    :param K: cyclotomic number field K
    :param r: K-rank r of the module lattice
    :param estimate: Estimation method. Default = "spherical model", alternative = "upper bound"
    "returns: the estimated skewness
    """    

    if estimate == "upper bound":
        return 0

    if estimate == "spherical model":
        d = K.deg
        dR = K.nb_real_embeddings
        dC = K.nb_complex_embeddings//2
        logNbb = (dR * digamma(r/2) + 2*dC * (digamma(r) - log(2)) - d*digamma(d*r/2))
        return logNbb/(2*d) + log(d)/2

    raise ValueError(f"skewness estimate method '{estimate}' not recognized")




def avg_lindex(K, r, estimate="density estimate"):
    """
    Estimate the average log index term in the slope equation (Term 4)
    Namely log N(I) following Section 4.5
    :param K: cyclotomic number field K
    :param r: K-rank r of the module lattice
    :param estimate: Estimation method. Default = "density estimate", alternative = "upper bound"
    "returns: the estimated log index
    """    

    if estimate == "upper bound":
        return 0
    if estimate == "density estimate":
        return K.logderivative_zeta(r)/K.deg

skewness_estimate={"upper bound": "upper bound",
                   "lower estimate": "spherical model"}
index_estimate   ={"upper bound": "upper bound",
                   "lower estimate": "density estimate"}

def slope_K(K, beta, estimate="upper bound"):
    """
    Estimate the average slope of the K-profile of module-BKZ over K (Section 4, Equation 4)
    :param K: cyclotomic number field K
    :param beta: blocksize (as a K-rank)
    :param estimate: Estimation method. Default = "upper bound", alternative = "lower estimate"
    "returns: the estimated slope_K
    """    
    d = K.deg
    terms = lghK(K, beta)                       # t_1
    terms += log(abs(K.disc))/(2*d) - log(d)/2  # t_2
    terms += avg_lskewness(K, beta, skewness_estimate[estimate]) # t_3
    terms += avg_lindex(K, beta, index_estimate[estimate]) # t_4
    return - (2*d/(beta-1)) * terms

def slope_Q(K, beta, estimate="upper bound"):
    """
    Estimate the average slope of the Q-profile of module-BKZ over K (Section 4, Equation 3)
    :param K: cyclotomic number field K
    :param beta: blocksize (as a K-rank)
    :param estimate: Estimation method. Default = "upper bound", alternative = "lower estimate"
    "returns: the estimated slope_Q
    """    
    return slope_K(K, beta, estimate) / (K.deg**2)
