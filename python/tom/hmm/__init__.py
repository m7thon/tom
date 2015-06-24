from ._hmm import Hmm, EMStopCondition, random_HMM, convert_HMM_to_OOM, learn_EM
try:
    from ._hmm import ghmm
except:
    pass
