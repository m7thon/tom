from __future__ import (division, absolute_import, print_function, unicode_literals)
#from .. import _tomlib

import numpy as np
try:
    import scipy.linalg as py_linalg
    pinv = py_linalg.pinv2
except ImportError:
    import numpy.linalg as py_linalg
    pinv = py_linalg.pinv


def spectral_norm_expectation(sqrt_V, n = 5):
    """
    Compute the expectation of the spectral norm of a random matrix with independently normally
    distributed entries with zero mean and standard deviation given by `sqrt_V` (i.e., the
    square root of the element-wise variances) from `n` samples. This returns a tuple of the
    average spectral norm and the estimated standard deviation.
    """
    spec_norm = np.zeros(n)
    for i in range(n):
        spec_norm[i] = np.linalg.norm(sqrt_V * np.random.randn(*sqrt_V.shape), ord=2)
    return np.mean(spec_norm), np.std(spec_norm)


def wsvd(sqrt_w_Y, F, sqrt_w_X):
    """Return the singular value decomposition `(U, s, VT)` of the row and column
    weighted matrix `sqrt_w_Y * F * sqrt_w_Y."""
    return np.linalg.svd(sqrt_w_Y * F * sqrt_w_X, full_matrices=False)


class CachedWSVD:
    """A cached version of `wsvd`."""
    def __init__(self):
        self._F = None
        self._sqrt_w_Y = None
        self._sqrt_w_X = None
        self._wsvd = None
        self._svd = None

    def __call__(self, sqrt_w_Y, F, sqrt_w_X):
        if F is not self._F:
            self._F = F
            self._F.setflags(write=False)
            self._svd = self._wsvd = None
        if np.array_equal(sqrt_w_Y, 1) and np.array_equal(sqrt_w_X, 1):
            if self._svd is None:
                U, s, VT = wsvd(1, F, 1)
                U.setflags(write=False)
                s.setflags(write=False)
                VT.setflags(write=False)
                self._svd = U, s, VT
            return self._svd
        if not (sqrt_w_Y is self._sqrt_w_Y or np.array_equal(sqrt_w_Y, self._sqrt_w_Y)):
            self._sqrt_w_Y = sqrt_w_Y
            try: self._sqrt_w_Y.setflags(write=False)
            except: pass
            self._wsvd = None
        if not (sqrt_w_X is self._sqrt_w_X or np.array_equal(sqrt_w_X, self._sqrt_w_X)):
            self._sqrt_w_X = sqrt_w_X
            try: self._sqrt_w_X.setflags(write=False)
            except: pass
            self._wsvd = None
        if self._wsvd is None:
            U, s, VT = wsvd(sqrt_w_Y, F, sqrt_w_X)
            U.setflags(write=False)
            s.setflags(write=False)
            VT.setflags(write=False)
            self._wsvd = U, s, VT
        return self._wsvd


cached_wsvd = CachedWSVD()
