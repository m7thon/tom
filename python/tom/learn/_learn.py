from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib
# from .. import linalg

import numpy as np
try:
    import scipy.linalg as py_linalg
    pinv = py_linalg.pinv2
except ImportError:
    import numpy.linalg as py_linalg
    pinv = py_linalg.pinv
import itertools


def wsvd(sqrt_w_Y, F, sqrt_w_X):
    """Return the singular value decomposition `(U, s, VT)` of the row and column
    weighted matrix `sqrt_w_Y * F * sqrt_w_Y."""
    return py_linalg.svd(sqrt_w_Y * F * sqrt_w_X)


class CachedWSVD:
    """A cached version of `wsvd`."""
    def __init__(self):
        self._F = None
        self._sqrt_w_Y = None
        self._sqrt_w_X = None
        self._wsvd = None

    def __call__(self, F, sqrt_w_Y=1, sqrt_w_X=1):
        redo = False
        if F is not self._F:
            self._F = F
            self._F.setflags(write=False)
            redo = True
        if not (sqrt_w_Y is self._sqrt_w_Y or np.array_equal(sqrt_w_Y, self._sqrt_w_Y)):
            self._sqrt_w_Y = sqrt_w_Y
            self._sqrt_w_Y.setflags(write=False)
            redo = True
        if not (sqrt_w_X is self._sqrt_w_X or np.array_equal(sqrt_w_X, self._sqrt_w_X)):
            self._sqrt_w_X = sqrt_w_X
            self._sqrt_w_X.setflags(write=False)
            redo = True
        if redo is True or self._wsvd is None:
            U, s, VT = wsvd(sqrt_w_Y, F, sqrt_w_X)
            U.setflags(write=False)
            s.setflags(write=False)
            VT.setflags(write=False)
            self._wsvd = U, s, VT
        return self._wsvd


_cached_wsvd = CachedWSVD()


class Data:
    def __init__(self, sequence=None, cache_level=1):
        self._sequence = None
        self._nInputSymbols = self._nOutputSymbols = None
        self._stree = None
        self._regularization = {'vPC': 2, 'vMin': 3}
        self._estimator = None
        self._X = self._Y = None
        self._e = _tomlib.Sequence()
        self._E = _tomlib.Sequences(1)
        self._cache = {}
        if cache_level == 1:
            self.cache(['F_YX', 'V_YX'])
        elif cache_level == 2:
            self.cache(['F_YX', 'F_zYX', 'f_YE', 'f_EX', 'V_YX', 'V_zYX', 'v_YE', 'v_EX'])
        if sequence is not None:
            self.sequence(setTo=sequence)

    def cache(self, keys=None):
        self._cache = {}
        for key in keys:
            if key in ['F_zYX', 'V_zYX']:
                self._cache[key] = {}
            else:
                self._cache[key] = None

    def _reset_cache(self, keys=None):
        if keys is None: keys = self._cache
        else: keys = [k for k in keys if k in self._cache]
        for key in keys:
            if key in ['F_zYX', 'V_zYX']:
                self._cache[key] = {}
            else:
                self._cache[key] = None

    def sequence(self, setTo=None):
        if setTo is None:
            if self._sequence is None: raise ValueError('No sequence has been set.')
            return self._sequence
        else:
            if self._sequence is None or not setTo.hasPrefix(self._sequence, withSameAlphabet=True):
                self._nInputSymbols = setTo.nInputSymbols()
                self._nOutputSymbols = setTo.nOutputSymbols()
                self._stree = _tomlib.STree(setTo)
            else:
                self._stree.extendTo(setTo, False)
            self._sequence = setTo
            self._estimator = _tomlib.Estimator(self._stree)
            self._estimator.regularization(**self._regularization)

            self._X = self._Y = None
            self._reset_cache()

    def nInputSymbols(self, setTo=None):
        if setTo is not None: self._nInputSymbols = setTo
        else:
            if self._nInputSymbols is None: raise ValueError('No sequence has been set.')
            return self._nInputSymbols

    def nOutputSymbols(self, setTo=None):
        if setTo is not None: self._nOutputSymbols = setTo
        if self._nOutputSymbols is None: raise ValueError('No sequence has been set.')
        return self._nOutputSymbols

    def stree(self):
        if self._sequence is None: raise ValueError('No sequence has been set.')
        return self._stree

    def estimator(self):
        if self._sequence is None: raise ValueError('No sequence has been set.')
        return self._estimator

    def regularization(self, setTo=None):
        if setTo is None:
            return self._regularization
        if self._estimator is not None:
            self._estimator.regularization(**setTo)
        self._regularization = setTo
        self._reset_cache(['V_YX', 'V_zYX', 'v_YE', 'v_EX'])

    def X(self, setTo=None):
        if setTo is None:
            if self._X is None: raise ValueError('Indicative words have not been set.')
            return self._X
        self._reset_cache()
        if type(setTo) in [list, tuple]:
            if self._sequence is None: raise ValueError('No sequence has been set.')
            self._X = _tomlib.wordsFromData(self.stree(), *setTo)
        elif type(setTo) is dict:
            if self._sequence is None: raise ValueError('No sequence has been set.')
            self._X = _tomlib.wordsFromData(self.stree(), **setTo)
        else: self._X = setTo

    def Y(self, setTo=None):
        if setTo is None:
            if self._Y is None: raise ValueError('Characteristic words have not been set.')
            return self._Y
        self._reset_cache()
        if type(setTo) in [list, tuple]:
            if self._sequence is None: raise ValueError('No sequence has been set.')
            self._Y = _tomlib.wordsFromData(self.stree(), *setTo)
        elif type(setTo) is dict:
            if self._sequence is None: raise ValueError('No sequence has been set.')
            self._Y = _tomlib.wordsFromData(self.stree(), **setTo)
        else: self._Y = setTo

    def F_YX(self, setTo=None):
        if setTo is not None: self._cache['F_YX'] = setTo
        else:
            try: F = self._cache['F_YX']
            except: return self.estimator().f(self.Y(), self.X())
            if F is None:
                F = self.estimator().f(self.Y(), self.X())
                self._cache['F_YX'] = F
            return F

    def F_zYX(self, z, setTo=None):
        if setTo is not None:
            if 'F_zYX' not in self._cache: self._cache['F_zYX'] = {}
            self._cache['F_zYX'][tuple(z)] = setTo
        else:
            try: F = self._cache['F_zYX'][tuple(z)]
            except: return self._estimator.f(self.Y(), z[0], z[1], self.X())
            if F is None:
                F = self.estimator().f(self.Y(), z[0], z[1], self.X())
                self._cache['F_zYX'][tuple(z)] = F
            return F

    def f_EX(self, setTo=None):
        if setTo is not None: self._cache['f_EX'] = setTo
        else:
            try: f = self._cache['f_EX']
            except: return self.estimator().f(self._E, self.X())
            if f is None:
                f = self.estimator().f(self._E, self.X())
                self._cache['f_EX'] = f
            return f

    def f_YE(self, setTo=None):
        if setTo is not None: self._cache['f_YE'] = setTo
        else:
            try: f = self._cache['f_YE']
            except: return self.estimator().f(self.Y(), self._E)
            if f is None:
                f = self.estimator().f(self.Y(), self._E)
                self._cache['f_YE'] = f
            return f

    def V_YX(self, regularization=None, setTo=None):
        if setTo is not None: self._cache['V_YX'] = setTo
        else:
            if regularization is not None and regularization != self._regularization:
                estimator = _tomlib.Estimator(self.stree())
                estimator.regularization(regularization)
                return estimator.v(self.Y(), self.X())
            try: V = self._cache['V_YX']
            except: return self.estimator().v(self.Y(), self.X())
            if V is None:
                V = self.estimator().v(self.Y(), self.X())
                self._cache['V_YX'] = V
            return V

    def V_zYX(self, z, regularization=None, setTo=None):
        if setTo is not None:
            if 'V_zYX' not in self._cache: self._cache['V_zYX'] = {}
            self._cache['V_zYX'][tuple(z)] = setTo
        else:
            if regularization is not None and regularization != self._regularization:
                estimator = _tomlib.Estimator(self.stree())
                estimator.regularization(regularization)
                return estimator.v(self.Y(), z[0], z[1], self.X())
            try: V = self._cache['V_zYX'][tuple(z)]
            except: return self.estimator().v(self.Y(), z[0], z[1], self.X())
            if V is None:
                V = self.estimator().v(self.Y(), z[0], z[1], self.X())
                self._cache['V_zYX'][tuple(z)] = V
            return V

    def v_EX(self, regularization=None, setTo=None):
        if setTo is not None: self._cache['v_EX'] = setTo
        else:
            if regularization is not None and regularization != self._regularization:
                estimator = _tomlib.Estimator(self.stree())
                estimator.regularization(regularization)
                return estimator.v(self._E, self.X())
            try: v = self._cache['v_EX']
            except: return self.estimator().v(self._E, self.X())
            if v is None:
                v = self.estimator().v(self._E, self.X())
                self._cache['v_EX'] = v
            return v

    def v_YE(self, regularization=None, setTo=None):
        if setTo is not None: self._cache['v_YE'] = setTo
        else:
            if regularization is not None and regularization != self._regularization:
                estimator = _tomlib.Estimator(self.stree())
                estimator.regularization(regularization)
                return estimator.v(self.Y(), self._E)
            try: v = self._cache['v_YE']
            except: return self.estimator().v(self.Y(), self._E)
            if v is None:
                v = self.estimator().v(self.Y(), self._E)
                self._cache['v_YE'] = v
            return v


def v_Y_from_data(data, p=1, q=1, regularization='auto'):
    """Return a column vector (or scalar) v_Y measuring the row precision of the matrix `data.F_YX()`.

    If `p` is `None`, then `data.V_EY()**q` is returned. Otherwise, the precision is computed as the
    rowwise `p`-mean of `data.V_YX()` raised to the power `q`. Setting `q` to be zero will return 1.

    Finally, if `regularization` is `'auto'`, a suitable regularization setting is used for the variance
    estimation ({'vPC': 2, 'vMin': 3} if `p` is `None` and {'vPC': 0, 'vMin': 1e-15} otherwise).
    Otherwise, the provided setting of `regularization` is used (if `None`, this used the regularization
    setting from `data`).
    """
    if q == 0: return 1
    if p is None:
        if regularization == 'auto': regularization = {'vPC': 2, 'vMin': 3}
        return data.v_YE(regularization=regularization)**q
    if regularization == 'auto': regularization = {'vPC': 0, 'vMin': 1e-15}
    V = data.V_YX(regularization=regularization)
    return _tomlib.rowwiseMean(V, p)**q


def v_X_from_data(data, p=1, q=1, regularization='auto'):
    """Return a row vector (or scalar) v_Y measuring the column precision of the matrix `data.F_YX()`.

    If `p` is `None`, then `data.V_XE()**q` is returned. Otherwise, the precision is computed as the
    colwise `p`-mean of `data.V_YX()` raised to the power `q`. Setting `q` to be zero will return 1.

    Finally, if `regularization` is `'auto'`, a suitable regularization setting is used for the variance
    estimation ({'vPC': 2, 'vMin': 3} if `p` is `None` and {'vPC': 0, 'vMin': 1e-15} otherwise).
    Otherwise, the provided setting of `regularization` is used (if `None`, this used the regularization
    setting from `data`).
    """
    if q == 0: return 1
    if p is None:
        if regularization == 'auto': regularization = {'vPC': 2, 'vMin': 3}
        return data.v_EX(regularization=regularization)**q
    if regularization == 'auto': regularization = {'vPC': 0, 'vMin': 1e-15}
    V = data.V_YX(regularization=regularization)
    return _tomlib.colwiseMean(V, p)**q


def v_Y_v_X_from_data(data, p=1, q=1, regularization='auto'):
    """Return a tuple `(v_Y, v_X)` measuring the row and column precision of the matrix `data.F_YX()`,
    respectively, which are 2-dimensional arrays representing a row or column vector, or scalars.

    If `p` is `None`, then `data.V_XE()**q` and `data.V_EY()**q` are returned. Otherwise, the precision
    vectors are computed as the rowwise and columnwise `p`-mean of `data.V_YX()` raised to the power `q`.
    Setting `q` to be zero will return (1,1).

    Finally, if `regularization` is `'auto'`, a suitable regularization setting is used for the variance
    estimation ({'vPC': 2, 'vMin': 3} if `p` is `None` and {'vPC': 0, 'vMin': 1e-15} otherwise).
    Otherwise, the provided setting of `regularization` is used (if `None`, this used the regularization
    setting from `data`).
    """
    if q == 0: return 1, 1
    if p is None:
        if regularization == 'auto': regularization = {'vPC': 2, 'vMin': 3}
        return data.v_EX(regularization=regularization)**q
    if regularization == 'auto': regularization = {'vPC': 0, 'vMin': 1e-15}
    V = data.V_YX(regularization=regularization)
    return _tomlib.rowwiseMean(V, p) ** q, _tomlib.colwiseMean(V, p) ** q


def rank_estimate(F, V, v_Y=1, v_X=1, errorNorm='spec', wsvd=_cached_wsvd):
    """Estimate the numerical rank of the matrix `F`. This uses the element-wise variances given by `V`
    to determine a cutoff for the error on `F` that is measured by the given `errorNorm`, which may
    be 'spec' (default) or 'frob'. """

    w_Y, w_X = 1/v_Y, 1/v_X
    wU, ws, wVT = wsvd(np.sqrt(w_Y), F, np.sqrt(w_X))

    if errorNorm == 'spec':
        e = np.sum(np.sqrt(w_Y * V * w_X)) / np.sqrt(V.size)
        d = 0
        while d < ws.size and ws[d] > e:
            d += 1
        return d
    elif errorNorm == 'frob':
        e = np.sum(w_Y * V * w_X)
        s2_tailsum = 0
        d = ws.size
        if d == 0: return d
        while d > 0 and s2_tailsum <= e:
            s2_tailsum += ws[d-1]**2
            d -= 1
        return d+1
    else:
        raise ValueError('Unknown errorNorm.')


def subspace_from_model(model, Y, v_Y=1, stabilization=None):
    if stabilization is None: stabilization = {'preset': 'none'}
    room = model.reverse(False)
    room.stabilization(**stabilization)
    Pi = np.zeros(len(Y), model.dimension())
    for i, y in enumerate(Y):
        log2_f_y = room.log2_f(y.reverse())
        Pi[i, :] = room.wt().transpose()
        if v_Y is not None: Pi[i, :] *= 2**log2_f_y
    return Pi


def subspace_by_svd(F, dim, v_Y=1, v_X=1, wsvd=_cached_wsvd):
    threshold = 1e-12
    sqrt_v_Y, sqrt_v_X = np.sqrt(v_Y), np.sqrt(v_X)
    U, s, VT = wsvd(1/sqrt_v_Y, F, 1/sqrt_v_X)
    dim = min(dim, len(s))
    ssd = np.sqrt(s[:dim])
    while dim > 1 and ssd[dim-1] < threshold: dim -= 1
    return sqrt_v_Y * U[:, :dim] * ssd[None, :dim]


def subspace_by_alternating_projections(F, dim_subspace, V, stopCondition=_tomlib.StopCondition(100, 1e-5, 1e-7), method='Cholesky'):
    if type(dim_subspace) is int:
        dim_subspace = subspace_by_svd(F, dim_subspace)
    return _tomlib.computeWLRA(F, 1/V, dim_subspace, stopCondition, method)


def subspace_corresponding_to_C_and_v_Y(C, v_Y):
    return v_Y * C.transpose()


def CQ(F_YX, dim_subspace, v_Y=1, v_X=1, wsvd=_cached_wsvd):
    sqrt_w_X = 1/np.sqrt(v_X)
    if type(dim_subspace) is int:
        dim = dim_subspace
        sqrt_w_Y = 1/np.sqrt(v_Y)
        threshold = 1e-12
        U, s, VT = wsvd(sqrt_w_Y, F_YX, sqrt_w_X)
        dim = min(dim, len(s))
        ssdi = np.sqrt(s[:dim])  # ssdi: square-root of Sd inverse
        for i in range(dim):
            if ssdi[i] > threshold: ssdi[i] = 1/ssdi[i]
            else: dim = i; break
        Ud = U[:, :dim]; ssdi = ssdi[:dim]; VdT = VT[:dim, :]
        C = ssdi[:, None] * (sqrt_w_Y * Ud).transpose()
        Q = (VdT * sqrt_w_X).transpose() * ssdi[None, :]
    else:
        subspace = dim_subspace
        C = (1 / v_Y * subspace).transpose()
        Q = np.transpose(sqrt_w_X) * pinv(C.dot(F_YX) * sqrt_w_X)
    return C, Q


def model_by_learning_equations(data, C, Q, CFQ_is_identity=True):
    dim = C.shape[0]
    nU = data.nInputSymbols()
    nO = data.nOutputSymbols()
    CFQ_inv = 1 if CFQ_is_identity is True else pinv(C.dot(data.F_YX()).dot(Q))
    oom = _tomlib.Oom(dim, nO, nU)
    for o, u in itertools.product(range(nO), range(max(1, nU))):
        oom.tau(o, u, C.dot(data.F_zYX([o, u])).dot(Q).dot(CFQ_inv))
    oom.sig(data.f_EX().dot(Q).dot(CFQ_inv))
    oom.w0(C.dot(data.f_YE()))
    oom.initialize()
    return oom


def model_by_weighted_equations(data, subspace, use_covariances=True):
    B = subspace
    dim = B.shape[1]
    A = _tomlib.solveLS(B, data.F_YX, 1/data.V_YX, transposed=False)
    oom = _tomlib.Oom(dim, data.nOutputSymbols(), data.nInputSymbols())
    for u, o in itertools.product(range(max(1, data.nInputSymbols())), range(data.nOutputSymbols())):
        Wz = 1/data.V_zYX([o, u])
        Az = _tomlib.solveLS(B, data.F_zYX([o, u]), Wz, transposed=False)
        WAz = _tomlib.transformWeights(Wz, B, covariances=use_covariances)
        oom.tau(o, u, _tomlib.solveLS(A, Az, WAz, transposed=True))
    oom.sig(_tomlib.solveLS(A, data.f_EX(), 1/data.v_EX(), transposed=True))
    oom.w0(_tomlib.solveLS(B, data.f_YE(), 1/data.v_YE(), transposed=False))
    oom.initialize()
    return oom


def model_estimate(data, dim_subspace=None, method='SPEC', v=None,
                   WLRAstopCondition=_tomlib.StopCondition(100, 1e-5, 1e-7), ES_stabilization=None, return_subspace=False):
    """Estimate a model from the given `data` using a spectral `method`.

    Parameters
    ----------
    data : tom.Data (or equivalent object)
        An representation of the training data that provides at least:
            - nOutputSymbols(), nInputSymbols()
            - F_YX(), F_zYX(z), f_YE(), f_EX()
        Additionally, depending on the requested algorithm, also:
            - V_YX(), V_zYX(z), v_YE(), v_EX()
            - Y()
    dim_subspace : int (or np.array), optional
        The model target dimension (determined numerically by default), or
        the basis of the principal subspace to use.
    method : { 'SPEC', 'RCW', 'ES', 'WLS', 'GLS' }
        The learning method to use:
            - 'SPEC': Standard spectral learning
            - 'RCW' : Row and column weighted spectral learning
            - 'ES'  : (Generalized) efficiency sharpening
            - 'WLS' : Weighted spectral learning using WLS
            - 'GLS' : Weighted spectral learning using GLS
    v : { (v_Y, v_X) }, optional
        The row / column precision to use in the case of 'RCW' or 'ES',
        which defaults to row and column averages of the elementwise
        variances for `data.F_YX()` computed *without* regularization.

        `v_[Y|X]` may each be given as:
            - a scalar or 2D column / row vector
            - a tuple of parameters. In this case `v_[Y|X]` will be
              computed by a call to `v_[Y|X]_from_data(data, *v_[Y|X])`.

        In the case of `ES`, `v_Y` may be given as `None` to mean that the
        row precisions `v_Y` are derived from the previous model estimate,
        which corresponds to using the "reversed characterizer".

        Furthermore, `v_X` may be omitted to indicate that the same settings
        should be used as for `v_Y`. For instance, one may pass
            `v = [[p,q]]`  or  `v = ((p,q), )`  (but *not* `v = ((p,q))`!).
    WLRAstopCondition : tom.util.StopCondition (or None), optional
        Determines the stopping condition for the iterative computation of
        the weighted low-rank approximation of F in the case of methods
        'WLS' or 'GLS'. By default at most 100 iterations are performed.
        This can be set to `None` to perform no iterations (which then uses
        the provided subspace or computes the principal subspace by SVD).
    ES_stabilization : dict of tom.Oom.stabilization() parameters
        Stabilization parameters to use for the computation of Π (defaults
        to no stabilization) in the case of method `ES`.
    return_subspace : bool, optional
        If `True` (default: `False`), return the principal subspace estimate.

    Returns
    -------
    oom : tom.Oom
        The estimated (IO)-OOM.

    Notes
    -----
    If the dimension or subspace is not specified (`dim_subspace = None`), an
    appropriate target dimension is selected based on the numerical rank
    of the matrix `data.F_YX()`.

    References
    ----------
    .. [1] M. Thon, H. Jaeger, *Links Between Multiplicity Automata,
    Observable Operator Models and Predictive State Representations -- a
    Unified Learning Framework*, JMLR, 2014, pg. 15"""

    if method in ['Standard', 'Spectral']: method = 'SPEC'
    elif method in ['RowColWeighted', 'RCCQ']: method = 'RCW'
    elif method == 'EfficienySharpening': method = 'ES'

    if dim_subspace is None:
        raise ValueError('Please provide a target dimension.')

    def _parse_v(data, v):
        if v == 1: return 1, 1
        if type(v) in [list, tuple]:  # row and column weights
            if len(v) == 1:  # same settings for v_Y and v_X
                return v_Y_v_X_from_data(data, *v[0]) if type(v[0]) is list else (v[0], v[0])
            elif len(v) == 2:  # separate settings for v_Y and v_X
                v_Y, v_X = v
                if type(v_Y) in [list, tuple]:
                    v_Y = v_Y_from_data(data, *v_Y)
                if type(v_X) in [list, tuple]:
                    v_X = v_X_from_data(data, *v_X)
                return v_Y, v_X
        raise ValueError('Unrecognized v.')

    if method == 'SPEC':
        v = (1, 1)
        method = 'RCW'

    if method == 'RCW':
        if v is None: v = ((1, 1),)
        v_Y, v_X = _parse_v(data, v)
        C, Q = CQ(data.F_YX(), dim_subspace, v_Y, v_X)
        model = model_by_learning_equations(data, C, Q)
        subspace = subspace_corresponding_to_C_and_v_Y(C, v_Y) if return_subspace and type(dim_subspace) is int else dim_subspace

    elif method == 'ES':
        if v is None: v = [[1, 1]]
        v_Y, v_X = _parse_v(data, v)
        if type(dim_subspace) is int: dim_subspace = model_estimate(data, dim_subspace, method='SPEC')
        subspace = subspace_from_model(dim_subspace, data.Y(), v_Y=v_Y, stabilization=ES_stabilization)
        if v_Y is None: v_Y = 1
        model = model_by_learning_equations(data, *CQ(data.F_YX(), subspace, v_Y, v_X))

    elif method in ['GLS', 'WLS']:
        if WLRAstopCondition is None: WLRAstopCondition = _tomlib.StopCondition(0)
        subspace = subspace_by_alternating_projections(data.F_YX(), dim_subspace, data.V_YX(), stopCondition=WLRAstopCondition)
        model = model_by_weighted_equations(data, subspace, method == 'GLS')

    else:
        raise ValueError('Unknown method: ' + method)

    return (model, subspace) if return_subspace else model
