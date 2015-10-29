from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib
#from . import linalg

import numpy as np
import itertools

def numericalRank(M, V, weightExp = 1):
    """Estimate the numerical rank of the matrix ``M`` with element-wise
    variances ``V``.

    Parameters
    ----------
    M : np.array
        The estimated matrix estimator(Y,X). This is recomputed if not given
    V : np.array
        The element-wise variances of the entries in M
    
    Returns
    -------
    dim : int
        The numerical rank of ``M``.
    """
    
    if weightExp == 0:
        s = np.linalg.svd(M, compute_uv=False)
        err = np.sum( np.sqrt(V) ) / np.sqrt( V.size )
    else:
        tV = np.maximum( 1e-15, V)
        wI = ( 1/( np.average(tV, 1)[:,np.newaxis] ) )**(0.5 * weightExp)
        wJ = ( 1/( np.average(tV, 0)[np.newaxis,:] ) )**(0.5 * weightExp)
        s = np.linalg.svd(wI * M * wJ, compute_uv=False)
        err = np.sum( wI * np.sqrt(V) * wJ ) / np.sqrt( V.size )
    dim = 0
    while dim < s.size and s[dim] > err: dim += 1
    return dim


def estimateDimension(estimator, X, Y, eF = None, weightExp = 1):
    """Estimate the model dimension for the underlying matrix F.

    Parameters
    ----------
    estimator : tom.Estimator (or equivalent object)
        An object that provides estimates for the function values f(x).
    X : tom.Sequences (or list of tom.Sequence ?)
        The set of indicative sequences, each of type tom.Sequence.
    Y : tom.Sequences (or list of tom.Sequence ?)
        The set of characteristic sequences, each of type tom.Sequence.
    eF : np.array, optional
        The estimated matrix estimator(Y,X). This is recomputed if not given
        
    Returns
    -------
    dim : int
    The estimated model dimension.
    """

    if eF is None: eF = estimator.f(Y,X)
    varianceParams = { p: getattr(estimator, p) for p in EstimatorVarianceParams }
    [setattr(estimator, p, val) for p, val in EstimatorExactVarianceParams.items()];
    V = estimator.v(Y, X)
    [setattr(estimator, p, val) for p, val in varianceParams.items()];
    return numericalRank(eF, V, weightExp)


def identifySubspace(eF, dim, W = None):
    if W is None: # We don't use any weights
        U,s,VT = np.linalg.svd(eF)
        s = np.sqrt(s[:dim])
        B = U[:,:dim] * s[None,:]; A = s[:,None] * VT[:dim,:]
    elif type(W) in [list, tuple]: # row/column weights
        wI = np.sqrt(W[0]); wJ = np.sqrt(W[1])
        if wI.ndim == 1: wI = wI[:,None]
        if wJ.ndim == 1: wJ = wJ[None,:]
        U,s,VT = np.linalg.svd( wI * eF * wJ )
        s = np.sqrt(s[:dim])
        B = 1/wI * U[:,:dim] * s[None,:]
        A = s[:,None] * VT[:dim,:] * 1/wJ
    else: # element-wise weights
        wI = ( 1 / np.sqrt( np.average(1/W, 1) ) )[:,None]
        wJ = ( 1 / np.sqrt( np.average(1/W, 0) ) )[None,:]
        B,A = identifySubspace(eF, dim, [wI, wJ])
        B = B.copy(order='F'); A = A.copy(order='F')
        improveWLRA(B, A, eF, W, 1e-7, 1000, "LLT")
    return B,A


def learnSpectral(estimator, X, Y, dimension = 0, weightExp = 0):
    """Spectral learning algorithm for (IO)-OOMs.

    This function returns an (IO)-OOM that is estimated using the spectral
    learning algorithm. The training data is provided in the form of an
    "estimator" together with sets of indicative and characteristic
    sequences given by "indicatives" and "characteristics". The target
    dimension may be specified by "dimemsion", otherwise it will be
    selected numerically. Row and column weights will be taken into account
    if the "weightExp" parameter is set to e.g. 1.

    Parameters
    ----------
    estimator : tom.Estimator (or equivalent object)
        An object that provides estimates for the function values f(x).
    X : tom.Sequences (or list of tom.Sequence ?)
        The set of indicative sequences, each of type tom.Sequence.
    Y : tom.Sequences (or list of tom.Sequence ?)
        The set of characteristic sequences, each of type tom.Sequence.
    dimension : int, optional
        The model target dimension. Determined numerically by default.
    weightExp : {0, >0}, optional
        The exponent to apply to the row and column weights. If zero, no row
        and column weights are used.
        
    Returns
    -------
    oom : tom.Oom
    The estimated (IO)-OOM.
 
    Notes
    -----
    The data is provided in the form of an estimator object. This may be
    a tom.Estimator, or an object that provides the following interface:
        nU() -> returns the size of the input alphabet
        nO() -> returns the size of the observation alphabet
        f(Y,X) -> returns the matrix of estimates f(xy)
        v(Y,X) -> returns the corresponding variance estimates  
        f(Y,X,o,u=0) -> returns the matrix of estimates f(x[u]oy)
        v(Y,X,o,u=0) -> returns the corresponding variance estimates

    If the target dimension is not specified (``dimension = 0``), an
    appropriate target dimension is selected based on the numerical rank
    of the matrix F = f(Y,X).

    Simple row and column weights can be used, and are computed as the inverse
    of the row-/column-wise average variance estimates of F, taken to the power
    given in ``weightExp``. The same row and column weights will also be used in
    the second step of the algorithm when solving the equations
    F * tau_z = F_z for tau_z.
    
    References
    ----------
    .. [1] M. Thon, H. Jaeger, *Links Between Multiplicity Automata,
    Observable Operator Models and Predictive State Representations -- a
    Unified Learning Framework*, JMLR, 2014, pg. 15"""
    
    nU = estimator.nU(); nO = estimator.nO()
    e = _tomlib.Sequences(1)
    threshold = 1e-12

    dim = dimension
    if dim == 0: dim = estimateDimension(estimator, X, Y)

    # Obtain weights:
    if weightExp == 0:
        wI = np.ones((len(Y), 1))
        wJ = np.ones((1, len(X)))
    else:
        V = estimator.v(Y, X)
        wI = ( np.sum(V, 1)[:,None] / len(Y) )**-(0.5 * weightExp)
        wJ = ( np.sum(V, 0)[None,:] / len(X) )**-(0.5 * weightExp)

    U,s,V_T = np.linalg.svd(wI * estimator.f(Y, X) * wJ)
    Ud = U[:,:dim]; Sd_inv = s[:dim]; Vd_T = V_T[:dim,:]
    for i in range(dim):
        Sd_inv[i] = 1 / Sd_inv[i] if Sd_inv[i] > threshold else 0
    del U, s, V_T
    C = Ud.transpose() * wI.transpose()
    Q = wJ.transpose() * Vd_T.transpose() * Sd_inv

    # Get OOM:
    oom = _tomlib.Oom()
    oom.setSize(dim, nO, nU)
    for u, o in itertools.product(range(max(1,nU)), range(nO)):
        oom.tau(o,u, C.dot(estimator.f(Y, X, o, u)).dot(Q) )
    oom.sig( estimator.f(e, X).dot(Q) )
    oom.w0( C.dot(estimator.f(Y, e)) )
    # oom.w0( oom.stationaryState() )
    oom.init()
    return oom


def learnWeightedSpectral(estimator, X, Y, dimension = 0, method = 'GLS', subspace = None):
    nU = estimator.nU(); nO = estimator.nO()
    e = _tomlib.Sequences(1)
    dim = dimension
    if dim == 0: dim = estimateDimension(estimator, X, Y)

    if subspace is None:
        B,A = identifySubspace(estimator.f(Y,X), dim, 1/estimator.v(Y,X))
    else:
        [B,A] = subspace

    oom = _tomlib.Oom()
    oom.setSize(dim, nO, nU)
    
    for u, o in itertools.product(range(max(1,nU)), range(nO)):
        Wz = 1/estimator.v(Y,X,o,u)
        Az = solveFastWLS(B, estimator.f(Y,X,o,u), Wz, True)
        if method == "WLS":
            # 3.1) compute new weights WAz for Az by [WAz]j = diag(B^\top * [Wz]j * B)
            WAz = np.zeros(Az.shape);
            for j in range(Az.shape[1]):
                WAz[:,j] = np.diag((B.transpose() * Wz[:,j]).dot(B))
            # 3.2) Tz <-- WLS solution to Tz * A = Az with weights WAz
            oom.tau(o,u, solveFastWLS(A, Az, WAz) )
        else: # GLS
            # 3.1) compute new block-diagonal weight matrix for Az
            #      with blocks Waz_j = B^\top * [Wz]j * B
            WAz = np.zeros((Az.shape[0], Az.shape[0]*Az.shape[1]))
            for j in range(Az.shape[1]):
                WAz[:, j*dim:(j+1)*dim] = (B.transpose() * Wz[:,j]).dot(B)
            # 3.2) Tz <-- GLS solution to Tz * A = Az with weights WAz
            oom.tau(o,u, solveFastGLS(A, Az, WAz) )

    oom.sig( solveFastWLS(A, estimator.f(e,X), 1/estimator.v(e,X)) )
    oom.w0( solveFastWLS(B, estimator.f(Y,e), 1/estimator.v(Y,e), True) )
        
    oom.init()
    return oom


# Old ??

# r/c weighted spectal learning:
# input: sets of characteristic and indicative sequences, estimator
# output: model

# 1. Make sure ind[0] (= char[0]) = empty word
# 2. Get F and row wX and column weights wY, as well as wz for each z.
#    These will be 1/counts.
# 3. Set DX = diag(wX)**0.5 and DY = diag(wY)**0.5
# 5. Ud, Sd, Vd_T = d-truncated SVD of DY*F*DX -- best weighted rank-d approx.
# 5b. Alternatively, compute Ud, Sd, Vd_T = d-truncted SVD of DY*[F F_]*DX
# 6. Set C = Ud^T * DY -- maps cols to d-dimensional representation
# 7. Set A = Sd*Vd_T (first cols) = C*F*DX, and set Az = C*F_z*DX
# 8. Set B = [A; wz**0.5 * Azs; (1/N)**0.5 * fX * DX]
# 9. Set Ud, Sd, Vd_T = d-truncated SVD of B. Ud = [U; Uz1; ...; Uzn; Usig]
# 10. Set tau_z = Uzd * Ud^{-1} / wz**0.5
    
def simpleWeightedSpectral(est, Y, X, d, weights = 'average'):
    nO = est.nO_
    e = _tomlib.Sequences(1)
    threshold = 1e-15

        
    # Weights:
    wY = est.v(Y, e)**-0.5
    wX = est.v(e, X)**-0.5

    # Best weighted rank-d approx. to F
    U,S,V_T = np.linalg.svd(wY * est.f(Y, X) * wX)
    Ud = U[:,:d]; Sd_inv = S[:d]; Vd_T = V_T[:d,:]
    for i in range(d):
        Sd_inv[i] = 1 / Sd_inv[i] if Sd_inv[i] > threshold else 0
    del U, S, V_T
    
    C = Ud.transpose() * wY.transpose()
    Q = wX.transpose() * Vd_T.transpose() * Sd_inv

    # Get OOM:
    oom = _tomlib.Oom()
    oom.setSize(d, nO, 0)
    for z in range(nO):
        oom.tau(z, 0, C.dot(est.f(Y, X, z)).dot(Q))
    oom.sig(est.f(e,X).dot(Q))
    # oom.w0(C.dot(est.f(Y,e)))
    oom.w0(oom.stationaryState())
    oom.reset()
    return oom


def weightedSpectral(est, Y, X, d, valueErrorBias = 1):
    nO = est.nO_
    e = _tomlib.Sequences(1)
    
    # Weights:
    wY = est.v(Y, e)**-0.5
    wX = est.v(e, X)**-0.5
    wz = [(est.f(e, e, z)[0,0]+0.00001)/1.00001 for z in range(nO)]
    wz = [(wz[z] * (1-wz[z]))**-0.5 for z in range(nO)]
    # Use same column weights for Fz
    # wz = [1 for z in range(nO)]

    # Best weighted rank-d approx. to F
    U,S,V_T = np.linalg.svd(wY * est.f(Y, X) * wX)
    Ud = U[:,:d]; Sd = np.diag(S[:d]); Vd_T = V_T[:d,:]
    del U, S, V_T
    
    A = np.zeros(( (nO+1)*d+1, X.size() ))
    # A[:d,:] = Ud.transpose().dot(wY * est.f(Y,X) * wX)
    A[:d,:] = Sd.dot(Vd_T) * valueErrorBias
    for z in range(nO):
        A[(z+1)*d:(z+2)*d,:] = Ud.transpose().dot(wY * est.f(Y, X, z) * wX * wz[z] )
    A[-1:,:] = 1 * est.f(e, X) * wX # Note that the weight of f(e) should be 1

    U,S,V_T = np.linalg.svd(A)
    del A
    Ud = U[:,:d]; Sd = np.diag(S[:d]); Vd_T = V_T[:d,:]
    del U, S, V_T
    Uinv = np.linalg.inv(Ud[:d,:] / valueErrorBias)

    # Apinv = pinv(A[:d,:])
    
    # Get OOM:
    oom = _tomlib.Oom()
    oom.setSize(d, nO, 0)
    for z in range(nO):
        # oom.tau(z, 0, A[(z+1)*d:(z+2)*d,:].dot(Apinv) / wz[z] )
        oom.tau(z, 0, Ud[(z+1)*d:(z+2)*d,:].dot(Uinv) / wz[z] )
    # oom.sig(A[-1:, :].dot(Apinv))
    oom.sig(Ud[-1:, :].dot(Uinv))
    # ??? oom.w0(Ud.transpose().dot(wY * est.f(Y,e)))
    oom.w0(oom.stationaryState())
    oom.reset()
    return oom
