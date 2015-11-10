from __future__ import (division, absolute_import, print_function, unicode_literals)
from .. import _tomlib

import numpy as np
import itertools
import bz2, gzip, re, ast, sys

# The following makes the basic tom objects copyable and pickleable
# via their json serialization, thereby improving python integration 
if sys.version_info[0] == 2: import copy_reg as copyreg
else: import copyreg
import copy, pickle

def unpickle_Oom(s):
    oom = _tomlib.Oom(); oom.fromJSON(s); return oom;
def pickle_Oom(oom):
    return unpickle_Oom, (oom.toJSON(),)
def unpickle_Sequence(s):
    seq = _tomlib.Sequence(); seq.fromJSON(s); return seq;
def pickle_Sequence(seq):
    return unpickle_Sequence, (seq.toJSON(),)
def unpickle_Sequences(s):
    return _tomlib.Sequences(s)
def pickle_Sequences(seqs):
    return unpickle_Sequences, ([seq for seq in seqs],)
def unpickle_Policy(s):
    pol = _tomlib.Policy(); pol.fromJSON(s); return pol;
def pickle_Policy(pol):
    return unpickle_Policy, (pol.toJSON(),)
copyreg.pickle(_tomlib.Oom, pickle_Oom)
copyreg.pickle(_tomlib.Sequence, pickle_Sequence)
copyreg.pickle(_tomlib.Sequences, pickle_Sequences)
copyreg.pickle(_tomlib.Policy, pickle_Policy)

def load(filename):
    """Load an Oom, Sequence or Policy from file with given filename.
    If the filename ends in .bz2, assume bzip2 compression."""
    if filename.endswith('.bz2'):
        with bz2.BZ2File(filename) as f:
            content = f.read().decode('UTF-8')
    elif filename.endswith('.gz'):
        with gzip.open(filename) as f:
            content = f.read().decode('UTF-8')
    else:
        with open(filename, 'r') as f:
            content = f.read()

    if content.startswith('{"Type":"OOM"'):
        oom = _tomlib.Oom()
        oom.fromJSON(content)
        return oom
    if content.startswith('{"Type":"SEQUENCE"'):
        seq = _tomlib.Sequence()
        seq.fromJSON(content)
        return seq
    if content.startswith('{"Type":"POLICY"'):
        pol = _tomlib.Policy()
        pol.fromJSON(content)
        return pol    
    if content.startswith('# This file is a POMDP policy'):
        content = re.sub(r'#.*\n', '', content)
        content = re.sub(r' ([a-zA-Z]+) =>', r' "\1":', content)
        content = re.sub(r'\s+', '', content)
        planes = ast.literal_eval(content)['planes']
        policy = _tomlib.Policy()
        for p in planes:
            policy.addPlane(p['action'], p['entries'][0::2], p['entries'][1::2])
        return policy
    if content.startswith('{"Type":"POMDP"'):
        pomdp = ast.literal_eval(content)
        oom = _tomlib.Oom()
        oom.setSize(pomdp['dim'], pomdp['nO'], pomdp['nU'])
        for a, o in itertools.product(range(pomdp['nU']), range(pomdp['nO'])):
            oom.tau(o,a, np.diag(np.array(pomdp['O(a,sp,o)'][a])[:,o]).dot(np.matrix(pomdp['T(a,s,sp)'][a],dtype=np.double).transpose()))
        oom.sig(np.ones((1,pomdp['dim']),dtype=np.double))
        oom.w0(np.array([pomdp['b(s)']],dtype=np.double).transpose())
        oom.minPrediction_ = 0
        oom.normalizationTolerance_ = 10
        oom.maxSetback(0)
        oom.initialize()
        return oom
        
    return None

def save(object, filename):
    """Save an Oom or Sequence object to file with given filename.
    If the filename ends in .bz2 or .gz, compress using bzip2 or gzip."""
    if filename.endswith('.bz2'):
        with bz2.BZ2File(filename, 'w') as f:
            f.write(bytes(object.toJSON(), 'UTF-8'))
    elif filename.endswith('.gz'):
        with gzip.open(filename, 'wb') as f:
            f.write(bytes(object.toJSON(), 'UTF-8'))
    else:
        with open(filename, 'w') as f:
            f.write(object.toJSON())

def getBasis(oom, co=False, eps_ind = 1e-5, eps_zero=1e-10):
    def addToCL(oom, v, v_seq, cl, sl, co):
        for u, o in itertools.product(range(max(1,oom.nInputSymbols())), range(oom.nOutputSymbols())):
            v_n = v.dot(oom.tau(o,u)) if co else oom.tau(o, u).dot(v)
            cl.append(v_n)
            if oom.nInputSymbols() != 0:
                sl.append([u,o] + v_seq) if co else sl.append(v_seq + [u,o])
            else:
                sl.append([o] + v_seq) if co else sl.append(v_seq + [o])
    v = oom.sig() if co else oom.w0()
    v_n = v/np.linalg.norm(v)
    B = [v_n.reshape(oom.dimension())]; S = [[]]
    P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
    CL = []; SL = []
    addToCL(oom, v_n, [], CL, SL, co)
    while (len(CL)>0 and len(S) < oom.dimension()):
        v = CL.pop(0); s = SL.pop(0)
        if np.linalg.norm(v) > eps_zero:
            v_n = v/np.linalg.norm(v)
            p_v_n = v_n.dot(P) if co else P.dot(v_n)
            if np.linalg.norm(p_v_n - v_n) > eps_ind:
                B.append(v_n.reshape(oom.dimension())); S.append(s)
                P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
                addToCL(oom, v_n, s, CL, SL, co)
    B = np.matrix(B) if co else np.matrix(B).transpose()
    return B, S

def getGoodBasis(oom, co=False, eps_ind = 1e-5, eps_zero=1e-10, normalizeVectors=True):
    def addToCL(oom, l, L, co):
        for u, o in itertools.product(range(max(1,oom.nInputSymbols())), range(oom.nOutputSymbols())):
            v = l[0].dot(oom.tau(o,u)) if co else oom.tau(o, u).dot(l[0])
            if np.linalg.norm(v) > eps_zero:
                vn = v / np.linalg.norm(v_n) if normalizeVectors else v
                if oom.nInputSymbols() != 0:
                    seq = [u,o] + l[1] if co else l[1] + [u,o]
                else:
                    seq = [o] + l[1] if co else l[1] + [o]
                L.append([vn, seq, 0])
    def sortCL(L, P, co):
        for l in L:
            vp = l[0].dot(P) if co else P.dot(l[0])
            l[2] = np.linalg.norm(l[0] - vp)
        L[:] = [l for l in L if l[2] > eps_ind]
        L.sort(key= lambda l: l[2], reverse=True)
        
    v = oom.sig() if co else oom.w0()
    v_n = v / np.linalg.norm(v) if normalizeVectors else v
    B = [v_n.reshape(oom.dimension())]; S = [[]]
    P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
    L = []
    addToCL(oom, [v_n, [], 0], L, co)
    sortCL(L, P, co)
    while (len(L)>0 and len(S) < oom.dimension()):
        (v,s,x) = L.pop(0)
        B.append(v.reshape(oom.dimension())); S.append(s)
        P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
        addToCL(oom, (v,s,x), L, co)
        sortCL(L, P, co)
    B = np.matrix(B) if co else np.matrix(B).transpose()
    return B, S

def minimize(oom, eps=1e-5):
    B, S = getBasis(oom, False, eps)
    oom.conjugate(np.linalg.pinv(B), B)
    B, S = getBasis(oom, True, eps)
    oom.conjugate(B, np.linalg.pinv(B))
