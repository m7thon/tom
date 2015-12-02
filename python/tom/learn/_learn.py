from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib
#from .. import linalg

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


def identifySubspace(F, W = None, dim = 0):
    if W is None: # We don't use any weights
        U,s,VT = np.linalg.svd(F)
        s = np.sqrt(s[:dim])
        B = U[:,:dim] * s[None,:]; A = s[:,None] * VT[:dim,:]
    elif type(W) in [list, tuple]: # row/column weights
        wI = np.sqrt(W[0]); wJ = np.sqrt(W[1])
        if wI.ndim == 1: wI = wI[:,None]
        if wJ.ndim == 1: wJ = wJ[None,:]
        U,s,VT = np.linalg.svd(wI * F * wJ)
        s = np.sqrt(s[:dim])
        B = 1/wI * U[:,:dim] * s[None,:]
        A = s[:,None] * VT[:dim,:] * 1/wJ
    else: # element-wise weights
        wI = ( 1 / np.sqrt( np.average(1/W, 1) ) )[:,None]
        wJ = ( 1 / np.sqrt( np.average(1/W, 0) ) )[None,:]
        B,A = identifySubspace(F, dim, [wI, wJ])
        B = B.copy(order='F'); A = A.copy(order='F')
        improveWLRA(B, A, F, W, 1e-7, 1000, "LLT")
    return B,A


def spectral(estimator, X, Y, dimension = 0, subspace = None, method = 'Standard', p=-1, q=1, stopConditionWLRA = None):
    """Spectral learning algorithm for (IO)-OOMs.

    This function returns an (IO)-OOM that is estimated using the (weighted) spectral
    learning algorithm. The training data is provided in the form of an
    ``estimator`` together with sets of indicative and characteristic
    sequences given by ``X`` and ``Y``, respectively. The target
    dimension may be specified by ``dimemsion``, otherwise it will be
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
    weightExp : double, optional
        The exponent to apply to the weights for the regression step (default 1).
        
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

    if method not in ['Standard', 'RowColWeighted', 'WLS', 'GLS']:
        raise ValueError('unknown method: ' + method)

    nU = estimator.sequence().nInputSymbols()
    nO = estimator.sequence().nOutputSymbols()
    e = _tomlib.Sequences(1)
    threshold = 1e-12

    if method in ['Standard', 'RowColWeighted'] or subspace in [None, []]:
        if dimension == 0:
            raise ValueError('dimension estimation not yet implemented')
            dimension = estimateDimension(estimator, X, Y)

        if method == 'Standard':
            F = estimator.f(Y,X)
            wI = np.ones((len(Y), 1))
            wJ = np.ones((1, len(X)))
        else:
            F, V = estimator.fv(Y,X)
            W = 1/V
            if p is None:
                F = estimator.f(Y,X)
                wI = estimator.v(Y,e)**q
                wJ = estimator.v(e,X)**q
            else:
                wI = _tomlib.rowwiseMean(W, p)**q
                wJ = _tomlib.colwiseMean(W, p)**q

        U,s,VT = np.linalg.svd(wI * F * wJ, full_matrices=0)
        U = U[:,:dimension]; s = np.sqrt(s[:dimension]); VT = VT[:dimension,:]

        B = 1/wI * U * s[None,:]
        A = s[:,None] * VT * 1/wJ
        if type(subspace) is list:
            subspace.clear()
            subspace.append(B)
            subspace.append(A)

        if method in ['Standard', 'RowColWeighted']:
            for i in range(dimension):
                s[i] = 1 / s[i] if s[i] > threshold else 0
            C = s[:,None] * U.transpose() * wI.transpose()
            Q = wJ.transpose() * VT.transpose() * s[None,:]
            oom = _tomlib.Oom(dimension, nO, nU)
            for u, o in itertools.product(range(max(1,nU)), range(nO)):
                oom.tau(o,u, C.dot(estimator.f(Y, o, u, X)).dot(Q))
            oom.sig( estimator.f(e, X).dot(Q) )
            oom.w0( C.dot(estimator.f(Y, e)) )
            oom.initialize()
            oom.stabilization(preset='default')
            return oom

        del U, s, VT
    else:
        B, A = subspace
        if dimension == 0:
            dimension = B.shape[1]
        if dimension != B.shape[1]:
            raise ValueError('given dimension does not match given subspace')
        F, V = estimator.fv(Y,X)
        W = 1/V

    if stopConditionWLRA is None:
        stopConditionWLRA = _tomlib.StopCondition(100, 1e-5, 1e-7)
    _tomlib.improveWLRA(B, A, F, W, stopCondition = stopConditionWLRA)

    oom = _tomlib.Oom(dimension, nO, nU)
    for u, o in itertools.product(range(max(1,nU)), range(nO)):
        Wz = 1/estimator.v(Y, o, u, X)
        Az = _tomlib.solveLS(B, estimator.f(Y, o, u, X), Wz, transposed=False)
        WAz = _tomlib.transformWeights(Wz, B, covariances = (method == 'GLS'))
        oom.tau(o,u, _tomlib.solveLS(A, Az, WAz, transposed=True))
    oom.sig( _tomlib.solveLS(A, estimator.f(e,X), 1/estimator.v(e,X), transposed=True) )
    oom.w0( _tomlib.solveLS(B, estimator.f(Y,e), 1/estimator.v(Y,e), transposed=False) )

    #oom.w0( oom.stationaryState() )
    oom.initialize()
    oom.stabilization(preset='default')
    return oom



def rowColWeightedTLS(est, Y, X, d, valueErrorBias = 1):
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
