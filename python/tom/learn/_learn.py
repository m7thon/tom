from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib
#from .. import linalg

import numpy as np
import itertools


def rowColWeights(estimator, Y, z, X, p=-0.5, q=1):
    if p is None:
        e = _tomlib.Sequences(1)
        return (estimator.v(Y,e)**-q, estimator.v(e,z,X)**-q)
    else:
        V = estimator.v(Y,z,X)
    if p > 0:
        W = 1/V
        return (_tomlib.rowwiseMean(W, p)**q, _tomlib.colwiseMean(W, p)**q)
    else:
        return (_tomlib.rowwiseMean(V, -p)**-q, _tomlib.colwiseMean(V, -p)**-q)


def estimateDimension(estimator, X, Y, p = -0.5, q = 1, weighted=False, frob=False):
    """Estimate the model dimension for the underlying matrix F.

    Parameters
    ----------
    estimator : tom.Estimator (or equivalent object)
        An object that provides estimates for the function values f(x).
    X : tom.Sequences (or list of tom.Sequence ?)
        The set of indicative sequences, each of type tom.Sequence.
    Y : tom.Sequences (or list of tom.Sequence ?)
        The set of characteristic sequences, each of type tom.Sequence.

    Returns
    -------
    dim : int
    The estimated model dimension.
    """

    e = _tomlib.Sequences(1)
    vP = estimator.regularization()
    estimator.regularization(preset='none')
    F, Vexact = estimator.fv(Y,X)
    estimator.regularization(*vP)

    if weighted:
        wI, wJ = rowColWeights(estimator, Y,e,X, p, q)
    else:
        wI = np.ones((len(Y),1))
        wJ = np.ones((1,len(X)))
    sqrt_wI = np.sqrt(wI)
    sqrt_wJ = np.sqrt(wJ)

    wU, ws, wVT = np.linalg.svd(sqrt_wI * F * sqrt_wJ)
    if frob:
        wef = np.sum(wI * Vexact * wJ)
        s2_tailsum = 0
        d = ws.size
        if d == 0: return d
        while d > 0 and s2_tailsum <= wef:
            s2_tailsum += ws[d-1]**2
            d -= 1
        return d+1
    else:
        we = np.sum(sqrt_wI * np.sqrt(Vexact) * sqrt_wJ) / np.sqrt(F.size)
        d = 0
        while d < ws.size and ws[d] > we:
            d += 1
        return d

def simpleSpectral(estimator, X, Y, dimensions = [0], method ='SPEC', p=-0.5, q=1, vP = (0,0,1,1)):
    if method == 'Standard': method = 'SPEC'
    if method == 'RowColWeighted': method = 'RCW'
    if method not in ['SPEC', 'RCW']:
        raise ValueError('unsupported method: ' + method)

    nU = estimator.sequence().nInputSymbols()
    nO = estimator.sequence().nOutputSymbols()
    e = _tomlib.Sequence(0, nO, nU)
    E = _tomlib.Sequences()
    E.append(e)
    threshold = 1e-15

    if type(dimensions) in [list, tuple]:
        dims = dimensions
    else:
        dims = [dimensions]

    F = estimator.f(Y,X)
    if method == 'SPEC':
        sqrt_wI = np.ones((len(Y), 1))
        sqrt_wJ = np.ones((1, len(X)))
    else: # compute row / col weights and their square roots
        wI, wJ = rowColWeights(estimator, Y,e,X, p, q)
        sqrt_wI = np.sqrt(wI)
        sqrt_wJ = np.sqrt(wJ)

    U,s,VT = np.linalg.svd(sqrt_wI * F * sqrt_wJ, full_matrices=0)

    res = []
    for dim in dims:
        if dim == 0:
            raise ValueError('dimension estimation not yet implemented')
            dim = estimateDimension(estimator, X, Y)

        Ud = U[:, :dim]; sd = np.sqrt(s[:dim]); VdT = VT[:dim, :]
        try:
            for i in range(dim):
                sd[i] = 1 / sd[i] if sd[i] > threshold else 0
        except: # target dimension much too large
            oom = _tomlib.Oom(0, nO, nU)
            oom.stabilization(preset='default')
            res.append(oom)
        else:
            C = sd[:,None] * Ud.transpose() * sqrt_wI.transpose()
            Q = sqrt_wJ.transpose() * VdT.transpose() * sd[None,:]
            oom = _tomlib.Oom(dim, nO, nU)
            for u, o in itertools.product(range(max(1,nU)), range(nO)):
                oom.tau(o,u, C.dot(estimator.f(Y, o, u, X)).dot(Q))
            oom.sig( estimator.f(E, X).dot(Q) )
            oom.w0( C.dot(estimator.f(Y, E)) )
            oom.initialize()
            oom.stabilization(preset='default')
            res.append(oom)
    if type(dimensions) in [list, tuple]:
        return res
    else:
        return res[0]

def rcSpectralFromData(nO, nU, F, Fz, Fsig, Fw0, wI, wJ, dim):
    e = _tomlib.Sequence(0, nO, nU)
    E = _tomlib.Sequences()
    E.append(e)
    threshold = 1e-15

    sqrt_wI, sqrt_wJ = np.sqrt(wI), np.sqrt(wJ)

    U,s,VT = np.linalg.svd(sqrt_wI * F * sqrt_wJ, full_matrices=0)
    U = U[:, :dim]; s = np.sqrt(s[:dim]); VT = VT[:dim, :]

    try:
        for i in range(dim):
            s[i] = 1 / s[i] if s[i] > threshold else 0
    except: # target dimension much too large
        oom = _tomlib.Oom(0, nO, nU)
        oom.stabilization(preset='default')
        return oom
    C = s[:,None] * U.transpose() * sqrt_wI.transpose()
    Q = sqrt_wJ.transpose() * VT.transpose() * s[None,:]
    oom = _tomlib.Oom(dim, nO, nU)
    for z, zid in zip(_tomlib.wordsOverAlphabet(nO, nU), range(10000)):
        oom.tau(z, C.dot(Fz[zid]).dot(Q))
    oom.sig(Fsig.dot(Q))
    oom.w0( C.dot(Fw0) )
    oom.initialize()
    oom.stabilization(preset='default')
    return oom


def spectral(estimator, X, Y, dimension = 0, subspace = None, method = 'SPEC', p=-0.5, q=1, stopConditionWLRA = None):
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
    method : { 'SPEC', 'RCCQ', 'RCW', 'WLS', 'GLS' }
        The (weighted) spectral learning method to use. The default 'SPEC'
        is the standard spectral learning with no weights.
    subspace : [ B, A ], optional
        Provides a method to cache the decomposition F = BA. Safe to ignore.
    p, q : double, optional
        Control the averaging used to obtain row / column weights.
    stopConditionWLRA : tom.util.StopCondition, optional
        Determines the stopping condition for the iterative computation of
        the weighted low-rank approximation of F in the case of methods
        'WLS' or 'GLS'. Safe to leave at the default.

        
    Returns
    -------
    oom : tom.Oom
    The estimated (IO)-OOM.
 
    Notes
    -----
    If the target dimension is not specified (``dimension = 0``), an
    appropriate target dimension is selected based on the numerical rank
    of the matrix F = f(Y,X).

    References
    ----------
    .. [1] M. Thon, H. Jaeger, *Links Between Multiplicity Automata,
    Observable Operator Models and Predictive State Representations -- a
    Unified Learning Framework*, JMLR, 2014, pg. 15"""

    if method == 'Standard': method = 'SPEC'
    if method == 'RowColWeighted': method = 'RCW'
    if method not in ['SPEC', 'RCCQ', 'RCW', 'WLS', 'GLS']:
        raise ValueError('unknown method: ' + method)

    nU = estimator.sequence().nInputSymbols()
    nO = estimator.sequence().nOutputSymbols()
    e = _tomlib.Sequence(0, nO, nU)
    E = _tomlib.Sequences()
    E.append(e)
    threshold = 1e-12

    if method in ['SPEC', 'RCW', 'RCCQ'] or subspace in [None, []]:
        if dimension == 0:
            # raise ValueError('dimension estimation not yet implemented')
            dimension = estimateDimension(estimator, X, Y)

        F = estimator.f(Y,X)

        if method == 'SPEC':
            sqrt_wI = np.ones((len(Y), 1))
            sqrt_wJ = np.ones((1, len(X)))
        else: # compute row / col weights and their square roots
            sqrt_wI, sqrt_wJ = rowColWeights(estimator, Y,e,X, p, 0.5*q)

        U,s,VT = np.linalg.svd(sqrt_wI * F * sqrt_wJ, full_matrices=0)
        U = U[:,:dimension]; s = np.sqrt(s[:dimension]); VT = VT[:dimension,:]

        B = 1/sqrt_wI * U * s[None,:]
        A = s[:,None] * VT * 1/sqrt_wJ
        if type(subspace) is list:
            subspace.clear()
            subspace.append(B)
            subspace.append(A)

        if method in ['SPEC', 'RCCQ']:
            for i in range(dimension):
                s[i] = 1 / s[i] if s[i] > threshold else 0
            C = s[:,None] * U.transpose() * sqrt_wI.transpose()
            Q = sqrt_wJ.transpose() * VT.transpose() * s[None,:]
            oom = _tomlib.Oom(dimension, nO, nU)
            for z in _tomlib.wordsOverAlphabet(nO, nU):
                oom.tau(z, C.dot(estimator.f(Y, z, X)).dot(Q))
            oom.sig( estimator.f(E, X).dot(Q) )
            oom.w0( C.dot(estimator.f(Y, E)) )
            oom.initialize()
            oom.stabilization(preset='default')
            return oom
        elif method in ['RCW']:
            oom = _tomlib.Oom(dimension, nO, nU)
            for z in _tomlib.wordsOverAlphabet(nO, nU):
                wIz, wJz = rowColWeights(estimator, Y, z, X, p, q)
                Az = _tomlib.solveLS(B, estimator.f(Y,z,X), wIz, transposed=False)
                oom.tau(z, _tomlib.solveLS(A, Az, wJz, transposed=True))
            oom.sig( _tomlib.solveLS(A, estimator.f(E,X), estimator.v(E,X)**-q, transposed=True) )
            oom.w0( _tomlib.solveLS(B, estimator.f(Y,E), estimator.v(Y,E)**-q, transposed=False) )
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
    oom.sig( _tomlib.solveLS(A, estimator.f(E,X), 1/estimator.v(E,X), transposed=True) )
    oom.w0( _tomlib.solveLS(B, estimator.f(Y,E), 1/estimator.v(Y,E), transposed=False) )

    #oom.w0( oom.stationaryState() )
    oom.initialize()
    oom.stabilization(preset='default')
    return oom



def RCWTLS(est, Y, X, d, valueErrorBias = 1):
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
