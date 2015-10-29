from __future__ import (division, absolute_import, print_function, unicode_literals)

import numpy as np

def TLS(A, B, valueErrorBias = 1):
    d = A.shape[0]
    M = np.zeros((d + B.shape[0], A.shape[1]))
    M[:d,:] = A * valueErrorBias
    M[d:,:] = B
    U, S, VT = np.linalg.svd(M)
    # to be continued ...
