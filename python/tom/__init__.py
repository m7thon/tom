from __future__ import division

from .tomlib import Random, Sequence, Sequences, STree, Hmm, Oom, Policy, Estimator

from .io import load, save
from .sequences import generateSequences, stringToSequence
from .learn import numericalRank, estimateDimension, identifySubspace, learnSpectral, learnWeightedSpectral
from .hmm import random_HMM, convert_HMM_to_OOM, learn_EM
#from .tools import *
