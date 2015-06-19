from .tomlib import *

import numpy as np
import itertools
import bz2, gzip, re, ast, sys

if sys.version_info[0] == 2: import copy_reg as copyreg
else: import copyreg
import copy, pickle

def unpickle_Oom(s):
    oom = Oom(); oom.fromJSON(s); return oom;
def pickle_Oom(oom):
    return unpickle_Oom, (oom.toJSON(),)
def unpickle_Sequence(s):
    seq = Sequence(); seq.fromJSON(s); return seq;
def pickle_Sequence(seq):
    return unpickle_Sequence, (seq.toJSON(),)
def unpickle_Sequences(s):
    return Sequences(s)
def pickle_Sequences(seqs):
    return unpickle_Sequences, ([seq for seq in seqs],)
def unpickle_Policy(s):
    pol = Policy(); pol.fromJSON(s); return pol;
def pickle_Policy(pol):
    return unpickle_Policy, (pol.toJSON(),)
copy_reg.pickle(Oom, pickle_Oom)
copy_reg.pickle(Sequence, pickle_Sequence)
copy_reg.pickle(Sequences, pickle_Sequences)
copy_reg.pickle(Policy, pickle_Policy)

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
        oom = Oom()
        oom.fromJSON(content)
        return oom
    if content.startswith('{"Type":"SEQUENCE"'):
        seq = Sequence()
        seq.fromJSON(content)
        return seq
    if content.startswith('{"Type":"POLICY"'):
        pol = Policy()
        pol.fromJSON(content)
        return pol
    
    # Keep support for old format for now...
    if content.startswith('OOM'):
        return Oom(content)
    if content.startswith('SEQUENCE'):
        return Sequence(content)
    
    if content.startswith('# This file is a POMDP policy'):
        content = re.sub(r'#.*\n', '', content)
        content = re.sub(r' ([a-zA-Z]+) =>', r' "\1":', content)
        content = re.sub(r'\s+', '', content)
        planes = ast.literal_eval(content)['planes']
        policy = Policy()
        for p in planes:
            policy.addPlane(p['action'], p['entries'][0::2], p['entries'][1::2])
        return policy
    if content.startswith('{"Type":"POMDP"'):
        pomdp = ast.literal_eval(content)
        oom = Oom()
        oom.setSize(pomdp['dim'], pomdp['nO'], pomdp['nU'])
        for a, o in itertools.product(range(pomdp['nU']), range(pomdp['nO'])):
            oom.tau(o,a, np.diag(np.array(pomdp['O(a,sp,o)'][a])[:,o]).dot(np.matrix(pomdp['T(a,s,sp)'][a],dtype=np.double).transpose()))
        oom.sig(np.ones((1,pomdp['dim']),dtype=np.double))
        oom.w0(np.array([pomdp['b(s)']],dtype=np.double).transpose())
        oom.minPrediction_ = 0
        oom.maxPredictionError_ = 10
        oom.maxSetback(0)
        oom.init()
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
