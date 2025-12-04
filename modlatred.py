import sys
sys.path.insert(1, './g6k')

from numpy import array, block, identity, zeros, roll, reshape, where
from numpy.random import randint
from fpylll import IntegerMatrix, FPLLL, BKZ, GSO, LLL
from math import log
from g6k import Siever
from cyclotomics import Cyclotomic
from g6k.utils.stats import dummy_tracer
from g6k.algorithms.workout import workout
from g6k.siever import Siever

def Z(a, b):
    """
    Returns an integral zero matrix of dimension a * b
    """
    return zeros((a, b), dtype=int)

def I(a):
    """
    Returns an integral identity matrix of dimension a * a
    """
    return identity(a, dtype=int)

def random_qary(m, n, q):
    """
    Returns the basis of a random q-ary lattice with m variables and n equations mod q
    """    
    assert(m >= n)
    A = randint(0, q-1, (m-n, n))
    B = block([[q*I(n), Z(n, m-n)], [A, I(m-n)]])
    return B

def random_qary_cyclotomic(K, m, n, q):
    """
    Returns the basis of a random q-ary module-lattice over K 
    with m variables and n equations mod q*O_K
    """    
    A = array([random_qary(m, n, q) for i in range(K.deg)])
    A[1:,:n,:] = 0; A[1:,:,n:] = 0
    R = block([[K.vOK_Zbasis(K.cyclic_embedding(A[:,i,j])) for j in range(m)] for i in range(m)])
    return R

def slope(l):
    """
    Returns the experimental slope of a given profile l
    """    
    n = len(l)
    i_mean = (n - 1) * 0.5
    l_mean = sum(l)/n
    v1, v2 = 0.0, 0.0

    for i in range(n):
        v1 += (i - i_mean) * (l[i] - l_mean)
        v2 += (i - i_mean) * (i - i_mean)
    return v1 / v2


class ModuleLatticeReduction(Siever):
    """
    A class derived from G6K Siever class, implementing the restructuring strategy
    from Section 3.3
    """

    def __init__(self, B, K, params=None, seed=None, debug=False, restructure_delta_prog=.05, float_type="ld"):
        """
        Extend the __init__ of G6K siever with extra parameters K, debug, 
        restructure_delta_prog, float_type, and apply restructuring.
        """
        self.debug = debug
        self.K = K
        self.float_type=float_type
        self.restructure_delta_prog = restructure_delta_prog
        M = self.MatGSO(IntegerMatrix.from_matrix(B), float_type=float_type)
        if (M.d % self.K.deg):
            raise ValueError("Rank of the structured lattice must be a multiple of K.deg")
        Siever.__init__(self, M, params, seed)
        self.restructure()

    def lll(self, l, r, delta=.99):
        """
        Override G6K Siever class lll method to give control of the parameter delta.
        Runs LLL from position l to r.
        :param l: Left end of the block to run LLL on (inclusive)
        :param r: Right end of the block to run LLL on (exclusive)
        :param delta: Reduction strength parameter delta in Lovasz condition
        """
        lll = LLL.Reduction(self.M, delta=delta)
        if not self.params.dual_mode:
            lll(l, l, r)
        else:
            m = self.full_n
            lll(m-r, m-r, m-l)

        self.initialized=False

    def position_check(self, kappa):
        """
        Check that position kappa is a legal position for a module operation.
        Raise an error if not.
        """
        if kappa%self.K.deg:
            raise ValueError("Structure check must be at positions multiples of K.deg")

    def is_structured_at(self, kappa, end=None, ret_w=False):
        """
        Test whether the basis has the adequate module structure at position kappa,
        that is check whether, after projection, the vectors B[kappa .. kappa+deg-1]
        indeed form the basis of a rank-1 module. (Section 3.3)

        This is done by checking that B[kappa] * O_K is a sublattice of 
        B[kappa .. kappa+deg-1]. We call Babai merely for expressing vectors in
        basis B. 
        :param kappa: Position at which to test the module structure
        :param end: Ignore the basis beyond that index; an optimization that can 
                    be used when we have been working in a limited range
        :param ret_w: Option to return a counterexample
        :return: A counterexample iff ret_w, otherwise a boolean
        """
        if self.K.deg==1:
            return True
        M, deg = self.M, self.K.deg
        end = end or M.d
        self.position_check(kappa)

        # Consider each basis vector v of B[kappa] * O_K
        # As an optimization, the first vector B[kappa] itself can be skipped
        for v in self.K.vOK_Zbasis(M.B[kappa])[1:]:
            # write v in base B
            w = array(M.babai(v, 0, end, gso=False))
            # test that w is indeed a linear combination of only B[kappa:kappa+deg]
            if (not any(w[kappa+1:kappa+deg])) or any(w[kappa+deg:end]):
                return w if ret_w else False
        return None if ret_w else True

    def is_structured(self):
        """
        Test whether the basis is structured at all positions
        :return: a boolean
        """
        if self.K.deg==1:
            return True
        for kappa in range(0, self.M.d, self.K.deg):
            if not self.is_structured_at(kappa):
                return False
        return True


    def restructure_at(self, kappa, end=None):
        """
        Enforce the module structure at position kappa (Section 3.3)
        Assumes structure is enforced outside of the range [kappa:end]
        :param kappa: Position at which to enforce the module structure
        :param end: Ignore the basis beyond that index; an optimization that can 
                    be used when we have been working in a limited range
        """
        if self.K.deg==1:
            return 

        end = end or self.M.d
        self.position_check(kappa)
        M, deg = self.M, self.K.deg

        # Assumes structure is enforced at prior positions
        for i in range(0, kappa, deg) if self.debug else []:
            if not self.is_structured_at(i):
                raise ValueError("Basis is not structured at positions %d < %d"%(i, kappa))

        n, m = M.B.nrows, M.B.ncols

        # Because of g6k+fpylll API (spec. the use of U and UInv in Mat.GSO)
        # we can not insert extra vectors to the basis, so we recreate a new 
        # matgso object without U/Uinv for that. We eliminate linear dependencies
        # with LLL, and reconstruct at the end.

        # In this process, LLL may break the structure again! Repeating and
        # decreasing delta as we go seems to eventually lead to the desired
        # structure when the conductor is small (say <= 16). To avoid this hack,
        # we would require fplll API to allow setting delta < 1/4 so as to 
        # do linear elimination while avoiding reduction side effects on the basis.
        trials = 0
        while not self.is_structured_at(kappa, end):
            B, X = Z(n, m), Z(end + self.K.deg - 1, m) 
            self.M.B.to_matrix(B)
            delta = .99 - self.restructure_delta_prog * trials
            if delta <= .3+ self.restructure_delta_prog and delta > .3:
                print("struggling with deg=%d at pos %d"%(self.K.deg, kappa), file=sys.stderr)

            # Create a set of generators of the current basis starting
            # with a basis of B[kappa] * O_K at positions B[kappa:kappa+deg]
            X[:kappa] = B[:kappa]
            X[kappa:kappa+deg] = self.K.vOK_Zbasis(B[kappa])
            X[kappa+deg:] = B[kappa+1:end]

            # Apply LLL to find the linear dependencies
            MX = GSO.Mat(IntegerMatrix.from_matrix(X), float_type=self.float_type)
            MX.update_gso()
            lll = LLL.Reduction(MX, delta=delta)
            lll(kappa, kappa, end+deg-1, 0)
            MX.B.to_matrix(X)
            # pull back the new basis, ignoring the first zero-vectors found by LLL            
            B[:end] = X[:end]

            # Reinitialize the Structured Siever with the new basis
            del self.M
            self.M = self.MatGSO(IntegerMatrix.from_matrix(B), float_type=self.float_type)
            self.lll(0, M.d, delta=.26)
            self.initialized = False
            del B, X
            trials += 1


    def restructure(self, start=0, end=None):
        """
        Enforce the module structure at all positions
        """
        if self.K.deg==1:
            return 
        end = end or self.M.d
        for kappa in range(start, end, self.K.deg):
            self.restructure_at(kappa, end)
        if self.debug:
            assert(self.is_structured())


    def svp_reduce(self, kappa=0, beta=None, full_restructure=True):
        """
        Run SVP reduction at position kappa while maintaining the module structure
        :param kappa: position
        :param beta: blocksize (expressed as a Q-rank, not a K-rank)
        :param full_restructure: Restructure all affected positions. The alternative
        is to restructure only at kappa, which is a viable optimization inside BKZ
        as the next position will be broken and repaired again by the next SVP call 
        """
        self.position_check(kappa)
        beta_ = min(beta or self.M.d, self.M.d - kappa)
        self.position_check(beta_)
        # G6K fails for absurdingly small dimensions, but LLL should suffice there
        if beta_ > 4:
            # Being less aggressive than in G6K to be quite sure we really solve SVP
            # We do a workout rather than just a pump, as typically done inside BKZ (pnj-BKZ)
            # and we take a bit less dims for free.
            d4f = int(min(max(0,(beta_ - 40)/2), 10 + 0.075*beta_))
            workout(self, dummy_tracer, kappa, beta_, d4f)

        self.lll(kappa, kappa+beta_)
        restruct_indices = range(kappa, kappa+beta_, self.K.deg) if full_restructure else [kappa] 
        for i in restruct_indices:
            self.restructure_at(kappa, kappa+beta_)

    def bkz(self, beta, tours=1, verbose=False):
        """
        Run tours of module BKZ
        :param beta: blocksize (expressed as a Q-rank, not a K-rank)
        :param tours: number of tours
        :param verbose: 
        """
        for t in range(tours):
            if verbose:
                print("tour %d / %d, cond=%d"%(t, tours, self.K.cond), file=sys.stderr)
            for kappa in range(0, self.M.d, self.K.deg):
                self.svp_reduce(kappa, beta, full_restructure=False)
            assert(self.is_structured())

    def hkz(self):
        """
        Module HKZ reduction
        """
        self.bkz(self.M.d, 1)

    def profile_Q(self, rescale=False):
        """
        Output the Q-profile of the current basis (Section 3.1)
        :param rescale: rescale the lattice to volume 1
        """        
        self.M.update_gso()
        prof = [log(self.M.get_r(i,i))/2 for i in range(self.M.d)]
        if not rescale:
            return prof
        avg = sum(prof)/len(prof)
        for i in range(len(prof)):
            prof[i] -= avg
        return prof


    def profile_K(self, rescale=False):
        """
        Output the K-profile of the current basis (Section 4.1)
        :param rescale: rescale the lattice to volume 1
        """        
        prof_Q = self.profile_Q(rescale)
        return [sum(prof_Q[i:i+self.K.deg]) for i in range(0, self.M.d, self.K.deg)]
