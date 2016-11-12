from __future__ import (division, absolute_import, print_function, unicode_literals)

import numpy as np

def spectral_norm_expectation(sqrt_V, n = 3):
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
