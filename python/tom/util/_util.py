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
    elif content.startswith('{"Type":"SEQUENCE"'):
        seq = _tomlib.Sequence()
        seq.fromJSON(content)
        return seq
    elif content.startswith('{"Type":"POLICY"'):
        pol = _tomlib.Policy()
        pol.fromJSON(content)
        return pol
    elif content.startswith('{"Type":"HMM"'):
        hmm = _tomlib.Hmm()
        hmm.fromJSON(content)
        return hmm
    elif content.startswith('# This file is a POMDP policy'):
        content = re.sub(r'#.*\n', '', content)
        content = re.sub(r' ([a-zA-Z]+) =>', r' "\1":', content)
        content = re.sub(r'\s+', '', content)
        planes = ast.literal_eval(content)['planes']
        policy = _tomlib.Policy()
        for p in planes:
            policy.addPlane(p['action'], p['entries'][0::2], p['entries'][1::2])
        return policy
    elif content.startswith('{"Type":"POMDP"'):
        pomdp = ast.literal_eval(content)
        hmm = _tomlib.Hmm(pomdp['dim'], pomdp['nO'], pomdp['nU'], 0)
        for a in range(pomdp['nU']):
            hmm.T(a, np.matrix(pomdp['T(a,s,sp)'][a],dtype=np.double))
        for a, o in itertools.product(range(pomdp['nU']), range(pomdp['nO'])):
            hmm.E(o, a, np.matrix(pomdp['O(a,sp,o)'][a], dtype=np.double)[:,o])
        hmm.pi(np.matrix([pomdp['b(s)']],dtype=np.double).transpose())
        hmm.init()
        return hmm

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

# def getBasis(oom, co=False, eps_ind = 1e-5, eps_zero=1e-10):
#     def addToCL(oom, v, v_seq, cl, sl, co):
#         for u, o in itertools.product(range(max(1,oom.nInputSymbols())), range(oom.nOutputSymbols())):
#             v_n = v.dot(oom.tau(o,u)) if co else oom.tau(o, u).dot(v)
#             cl.append(v_n)
#             if oom.nInputSymbols() != 0:
#                 sl.append([u,o] + v_seq) if co else sl.append(v_seq + [u,o])
#             else:
#                 sl.append([o] + v_seq) if co else sl.append(v_seq + [o])
#     v = oom.sig() if co else oom.w0()
#     v_n = v/np.linalg.norm(v)
#     B = [v_n.reshape(oom.dimension())]; S = [[]]
#     P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
#     CL = []; SL = []
#     addToCL(oom, v_n, [], CL, SL, co)
#     while (len(CL)>0 and len(S) < oom.dimension()):
#         v = CL.pop(0); s = SL.pop(0)
#         if np.linalg.norm(v) > eps_zero:
#             v_n = v/np.linalg.norm(v)
#             p_v_n = v_n.dot(P) if co else P.dot(v_n)
#             if np.linalg.norm(p_v_n - v_n) > eps_ind:
#                 B.append(v_n.reshape(oom.dimension())); S.append(s)
#                 P = np.linalg.pinv(np.matrix(B)) * np.matrix(B) if co else np.matrix(B).transpose() * np.linalg.pinv(np.matrix(B).transpose())
#                 addToCL(oom, v_n, s, CL, SL, co)
#     B = np.matrix(B) if co else np.matrix(B).transpose()
#     return B, S

def getBasis(oom, co=False, eps_independence = 1e-5, eps_zero=1e-12, sort='length'):
    v = oom.sig().transpose() if co else oom.w0()
    v = v / np.linalg.norm(v)
    candidates = [[_tomlib.Sequence(0, oom.nOutputSymbols(), oom.nInputSymbols()), v, 0]]
    basis = np.zeros((oom.dimension(), 0), dtype=np.double)
    basis_words = []
    while (len(candidates) > 0 and len(basis_words) < oom.dimension()):
        # pick next candidate on list
        (word, state, independence) = candidates.pop(0)
        basis = np.hstack((basis, state))
        basis_words.append(word)
        projection = basis.dot(np.linalg.pinv(basis))
        # add new candidates:
        for z in _tomlib.wordsOverAlphabet(oom.nOutputSymbols(), oom.nInputSymbols()):
            candidate_state = oom.tau(z).transpose().dot(state) if co else oom.tau(z).dot(state)
            if np.linalg.norm(candidate_state) > eps_zero:
                candidate_state /= np.linalg.norm(candidate_state)
                candidates.append([z + word if co else word + z, candidate_state, 0])
        # update independence of remaining candidates and sort:
        for candidate in candidates:
            candidate[2] = np.linalg.norm(projection.dot(candidate[1]) - candidate[1])
        candidates = [candidate for candidate in candidates if candidate[2] > eps_independence]
        if sort == 'length':
            candidates.sort(key= lambda c: (-len(c[0]), c[2]), reverse=True)
        elif sort == 'best':
            candidates.sort(key= lambda c: c[2], reverse=True)
    if co: basis = basis.transpose()
    return basis, basis_words

def minimize(oom, eps_independence=1e-5):
    B, words = getBasis(oom, co=False, eps_independence=eps_independence)
    oom.conjugate(np.linalg.pinv(B), B)
    B, words = getBasis(oom, co=True, eps_independence=eps_independence)
    oom.conjugate(B, np.linalg.pinv(B))

def computeStates(oom, words, co = False):
    """
    Return a matrix whose columns are the states of the given `oom` corresponding to the given `words`, or, if `co` is set to true, whose rows are the corresponding co-states.
    """
    room = oom.reverse(False) if co else oom
    S = np.zeros((oom.dimension(), len(words)))
    for i in range(len(words)):
        f = room.f(words[i].reverse() if co else words[i])
        S[:,i] = room.wt()[:,0] * f
    return S.transpose() if co else S

def stringToSequence(str, symbolTable=None):
    if not symbolTable: symbolTable = {}
    char = set()
    for c in str:
        if c not in symbolTable:
            char.add(c)
    char = sorted(list(char), key=lambda c: str.count(c), reverse=True)
    firstFreeSymbolValue = 0 if len(symbolTable) == 0 else max(symbolTable.values()) + 1
    for i in range(len(char)):
        symbolTable[char[i]] = firstFreeSymbolValue + i
    return _tomlib.Sequence([symbolTable[c] for c in str], firstFreeSymbolValue + len(char)), symbolTable
