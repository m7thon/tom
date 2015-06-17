from tomlib import *
import numpy as np
import scipy.linalg as linalg
import warnings
try:
    import ghmm
    ghmm.log.setLevel(ghmm.logging.ERROR)
except:
    warnings.warn('Could not load `ghmm` module. HMM functions will not work.')

def seq_iter(seq):
    for i in range(len(seq)):
        yield seq[i]

def random_HMM(alphabet_size, dim, exponent = 1):
    """Return randomly initialized HMM parameters (T, E, pi).

    Parameters
    ----------
    alphabet_size : int
        the size of the emission alphabet
    dim : int
        the number of hidden states
    exponent : int, double
        exponent to apply to all parameters before normalization

    Returns
    -------
    (T, E, pi) : tuple of
    T : np.array of size `dim` x `dim`
        the state transition matrix. T_{i,j} = P(s_j|s_i)
    E : np.array of size `dim` x `alphabet_size`
        the matrix of emisison probabilities. E_{i,j} = P(x_j|s_i)
    pi : np.array of size `dim`
        the initial state vector. pi_{i} = P(s_i)
    
    Notes
    -----
    All HMM parameters are drawn uniformly and independently from the interval [0,1]. Next, the given
    `exponent` is applied to each entry. Finally, the rows of the transition and emission matrices as
    well as the initial state distribution are normalized to sum to one.
    """

    pi = np.random.rand(dim)**exponent
    pi /= sum(pi)
    T = np.random.rand(dim,dim)**exponent
    for row in T: row /= sum(row)
    E = np.random.rand(dim,alphabet_size)**exponent
    for row in E: row /= sum(row)
    pi_old = np.zeros(dim)
    while linalg.norm(pi - pi_old) > 1e-12: pi_old = pi; pi = pi.dot(T); pi /= sum(pi)
    return (T, E, pi)

def convert_HMM_to_OOM(T, E, pi):
    """Convert a given HMM to an OOM representation.

    Parameters
    ----------
    T : np.array or array of size `dim` x `dim`
        the state transition matrix. T_{i,j} = P(s_j|s_i)
    E : np.array or array of size `dim` x `alphabet_size`
        the matrix of emisison probabilities. E_{i,j} = P(x_j|s_i)
    pi : np.array or list of size `dim`
        the initial state vector. pi_{i} = P(s_i)


    Returns
    -------
    oom : tom.Oom
    An OOM that is equivalent to the given HMM.

    Notes
    -----
    The OOM parameters are computed as
        sig = np.ones(1,dim)
        tau(z) = T.transopse().dot(np.diag(E[:,z]))
        w0 = pi[:,None]

    This implies that all OOM parameters are non-negative.
                    
    """

    T = np.array(T); E = np.array(E); pi = np.array(pi)
    dim = E.shape[0]; nO = E.shape[1]
    oom = Oom()
    oom.setSize(dim, nO, 0)
    for o in range(nO):   oom.tau( o,0, T.transpose() * E[:,o] )
    oom.sig( np.ones((1,dim)) )
    oom.w0( pi[:,None] )
    # oom.w0( oom.stationaryState() )
    oom.init()
    return oom

def learn_EM(train_seq, dim, n_init = 1, init_exponent = 1, max_iter = 100, min_improvement = 0.0001):
    """Train a HMM using the Baum-Welch algorithm.

    Parameters
    ----------
    train_seq : tom.Sequence
        the training sequence
    dim : int
        the number of hidden states for the HMM
    n_init : int
        the number of EM training runs to perform with new initializations (default: 1)
    init_exponent : int, double
        the exponent used for the initialization. See `random_HMM` (default: 2)
    max_iter : int
        the maximum number of EM iterations to perform in each training run (default: 100)
    min_improvement : double
        the EM iteration is terminated if the relative improvement in likelihood is below `min_improvement`
        (default: 0.0001)

    Returns
    -------
    oom : tom.Oom
        the estimated HMM converted to an OOM using `convert_HMM_to_OOM`
        
    Notes
    -----
    To avoid getting trapped in local minima, this funciton performs `n_init` training runs from new
    randomly initialized HMMs and selects the best model based on its likelihood on the training data.
    """

    alphabet = ghmm.IntegerRange(0, train_seq.nO())
    ghmm_train = ghmm.EmissionSequence(alphabet, list(seq_iter(train_seq)))
    (T, E, pi) = random_HMM(train_seq.nO(), dim, init_exponent)
    best_m = ghmm.HMMFromMatrices(alphabet, ghmm.DiscreteDistribution(alphabet), T, E, pi)
    best_m.baumWelch(ghmm_train, max_iter, min_improvement)
    if n_init > 1: best_LL = best_m.loglikelihood(ghmm_train)
    for it in range(n_init - 1):
        (T, E, pi) = random_HMM(train_seq.nO(), dim, init_exponent)
        m = ghmm.HMMFromMatrices(alphabet, ghmm.DiscreteDistribution(alphabet), T, E, pi)
        m.baumWelch(ghmm_train, max_iter, min_improvement)
        m_LL = m.loglikelihood(ghmm_train)
        if m_LL > best_LL: best_LL, best_m = m_LL, m
    m_params = best_m.asMatrices()
    return convert_HMM_to_OOM(m_params[0], m_params[1], m_params[2])

