from __future__ import division
import numpy as np
import scipy.linalg as linalg
import itertools
from tomlib import *
from tomio import load, save
from tomseqs import generateSequences, stringToSequence, iter
from tomlearn import numericalRank, estimateDimension, identifySubspace, learnSpectral, learnWeightedSpectral
from tomhmm import random_HMM, convert_HMM_to_OOM, learn_EM

def getBasis(oom, co=False, eps_ind = 1e-5, eps_zero=1e-10):
    def addToCL(oom, v, v_seq, cl, sl, co):
        for u, o in itertools.product(range(max(1,oom.nU())), range(oom.nO())):
            v_n = v.dot(oom.tau(o,u)) if co else oom.tau(o, u).dot(v)
            cl.append(v_n)
            if oom.nU() != 0:
                sl.append([u,o] + v_seq) if co else sl.append(v_seq + [u,o])
            else:
                sl.append([o] + v_seq) if co else sl.append(v_seq + [o])
    v = oom.sig() if co else oom.w0()
    v_n = v/linalg.norm(v)
    B = [v_n.reshape(oom.dim())]; S = [[]]
    P = linalg.pinv2(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * linalg.pinv2(np.matrix(B).transpose())
    CL = []; SL = []
    addToCL(oom, v_n, [], CL, SL, co)
    while (len(CL)>0 and len(S) < oom.dim()):
        v = CL.pop(0); s = SL.pop(0)
        if linalg.norm(v) > eps_zero:
            v_n = v/linalg.norm(v)
            p_v_n = v_n.dot(P) if co else P.dot(v_n)
            if linalg.norm(p_v_n - v_n) > eps_ind:
                B.append(v_n.reshape(oom.dim())); S.append(s)
                P = linalg.pinv2(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * linalg.pinv2(np.matrix(B).transpose())
                addToCL(oom, v_n, s, CL, SL, co)
    B = np.matrix(B) if co else np.matrix(B).transpose()
    return B, S

def getGoodBasis(oom, co=False, eps_ind = 1e-5, eps_zero=1e-10, normalizeVectors=True):
    def addToCL(oom, l, L, co):
        for u, o in itertools.product(range(max(1,oom.nU())), range(oom.nO())):
            v = l[0].dot(oom.tau(o,u)) if co else oom.tau(o, u).dot(l[0])
            if linalg.norm(v) > eps_zero:
                vn = v / linalg.norm(v_n) if normalizeVectors else v
                if oom.nU() != 0:
                    seq = [u,o] + l[1] if co else l[1] + [u,o]
                else:
                    seq = [o] + l[1] if co else l[1] + [o]
                L.append([vn, seq, 0])
    def sortCL(L, P, co):
        for l in L:
            vp = l[0].dot(P) if co else P.dot(l[0])
            l[2] = linalg.norm(l[0] - vp)
        L[:] = [l for l in L if l[2] > eps_ind]
        L.sort(key= lambda l: l[2], reverse=True)
        
    v = oom.sig() if co else oom.w0()
    v_n = v / linalg.norm(v) if normalizeVectors else v
    B = [v_n.reshape(oom.dim())]; S = [[]]
    P = linalg.pinv2(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * linalg.pinv2(np.matrix(B).transpose())
    L = []
    addToCL(oom, [v_n, [], 0], L, co)
    sortCL(L, P, co)
    while (len(L)>0 and len(S) < oom.dim()):
        (v,s,x) = L.pop(0)
        B.append(v.reshape(oom.dim())); S.append(s)
        P = linalg.pinv2(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * linalg.pinv2(np.matrix(B).transpose())
        addToCL(oom, (v,s,x), L, co)
        sortCL(L, P, co)
    B = np.matrix(B) if co else np.matrix(B).transpose()
    return B, S

def minimize(oom, eps=1e-5):
    B, S = getBasis(oom, False, eps)
    oom.conjugate(linalg.pinv2(B), B)
    B, S = getBasis(oom, True, eps)
    oom.conjugate(B, linalg.pinv2(B))


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

def TLS(A, B, valueErrorBias = 1):
    d = A.shape[0]
    M = np.zeros((d + B.shape[0], A.shape[1]))
    M[:d,:] = A * valueErrorBias
    M[d:,:] = B
    U, S, VT = linalg.svd(M)
    # to be continued ...
    
def simpleWeightedSpectral(est, Y, X, d, weights = 'average'):
    nO = est.nO_
    e = Sequences(1)
    threshold = 1e-15

        
    # Weights:
    wY = est.v(Y, e)**-0.5
    wX = est.v(e, X)**-0.5

    # Best weighted rank-d approx. to F
    U,S,V_T = linalg.svd(wY * est.f(Y, X) * wX)
    Ud = U[:,:d]; Sd_inv = S[:d]; Vd_T = V_T[:d,:]
    for i in range(d):
        Sd_inv[i] = 1 / Sd_inv[i] if Sd_inv[i] > threshold else 0
    del U, S, V_T
    
    C = Ud.transpose() * wY.transpose()
    Q = wX.transpose() * Vd_T.transpose() * Sd_inv

    # Get OOM:
    oom = Oom()
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
    e = Sequences(1)
    
    # Weights:
    wY = est.v(Y, e)**-0.5
    wX = est.v(e, X)**-0.5
    wz = [(est.f(e, e, z)[0,0]+0.00001)/1.00001 for z in range(nO)]
    wz = [(wz[z] * (1-wz[z]))**-0.5 for z in range(nO)]
    # Use same column weights for Fz
    # wz = [1 for z in range(nO)]

    # Best weighted rank-d approx. to F
    U,S,V_T = linalg.svd(wY * est.f(Y, X) * wX)
    Ud = U[:,:d]; Sd = np.diag(S[:d]); Vd_T = V_T[:d,:]
    del U, S, V_T
    
    A = np.zeros(( (nO+1)*d+1, X.size() ))
    # A[:d,:] = Ud.transpose().dot(wY * est.f(Y,X) * wX)
    A[:d,:] = Sd.dot(Vd_T) * valueErrorBias
    for z in range(nO):
        A[(z+1)*d:(z+2)*d,:] = Ud.transpose().dot(wY * est.f(Y, X, z) * wX * wz[z] )
    A[-1:,:] = 1 * est.f(e, X) * wX # Note that the weight of f(e) should be 1

    U,S,V_T = linalg.svd(A)
    del A
    Ud = U[:,:d]; Sd = np.diag(S[:d]); Vd_T = V_T[:d,:]
    del U, S, V_T
    Uinv = linalg.inv(Ud[:d,:] / valueErrorBias)

    # Apinv = pinv(A[:d,:])
    
    # Get OOM:
    oom = Oom()
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
