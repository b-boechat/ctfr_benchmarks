import numpy as np
from scipy.signal import correlate
from libc.math cimport sqrt
from ctfr.utils.arguments_check import _enforce_nonnegative, _enforce_odd_positive_integer
cimport cython

def _baseline_fls_wrapper(X, lk = 21, lm = 11, gamma = 20.0):

    lk = _enforce_odd_positive_integer(lk, 'lk', 21)
    lm = _enforce_odd_positive_integer(lm, 'lm', 11)
    gamma = _enforce_nonnegative(gamma, 'gamma', 20.0)

    return _baseline_fls_cy(X, lk, lm, gamma)

@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.nonecheck(False)
@cython.cdivision(True)
cdef _baseline_fls_cy(double[:,:,::1] X, Py_ssize_t lk, Py_ssize_t lm, double gamma):

    cdef:
        Py_ssize_t P = X.shape[0] # Spectrograms axis.
        Py_ssize_t K = X.shape[1] # Frequency axis.
        Py_ssize_t M = X.shape[2] # Time axis.

        double epsilon = 1e-10 # Small value used to avoid 0 in some computations. Bumped to 1e-10 to avoid numerical issues with cross-correlation.
        double window_size_sqrt = sqrt(<double> lk * lm)

    X_ndarray = np.asarray(X)

    # Local energy containers.
    cdef: 
        double[:,:] local_energy_l1
        double[:,:] local_energy_l2
        double[:,:] local_energy_l1_sqrt

    # Local suitability container.
    suitability_ndarray = np.empty((P, K, M), dtype=np.double)
    cdef double[:,:,:] suitability = suitability_ndarray

    # Containers related to combination.
    combination_weight_ndarray = np.empty((P, K, M), dtype=np.double)
    cdef double[:, :, :] combination_weight = combination_weight_ndarray

    cdef double[:, :] suitability_product

    # Generate the 2D window for local sparsity calculation.
    hamming_window = np.outer(np.hamming(lk), np.hamming(lm))

    ############ Local suitability calculation (using local Hoyer sparsity): {{{

    for p in range(P):
        # Calculate L1 and L2 local energy matrixes and element-wise square root of the L1 matrix.
        # Clipping is used because scipy correlate can return negative values for matrixes with very small positive elements.
        local_energy_l1_ndarray = np.clip(correlate(X_ndarray[p], hamming_window, mode='same'), a_min=epsilon, a_max=None)
        local_energy_l2_ndarray = np.sqrt(np.clip(correlate(X_ndarray[p] * X_ndarray[p], hamming_window * hamming_window, mode='same'), a_min=epsilon, a_max=None))
        local_energy_l1_sqrt_ndarray = np.sqrt(local_energy_l1_ndarray)

        # Point Cython memview to the calculated matrixes.
        local_energy_l1 = local_energy_l1_ndarray
        local_energy_l2 = local_energy_l2_ndarray
        local_energy_l1_sqrt = local_energy_l1_sqrt_ndarray

        # Calculate local suitability.
        for k in range(K):
            for m in range(M):
                suitability[p, k, m] = (window_size_sqrt - local_energy_l1[k, m]/local_energy_l2[k, m])/ \
                                        ((window_size_sqrt - 1) * local_energy_l1_sqrt[k, m]) + epsilon

    ############ }}}

    ############ Spectrograms combination {{{

    suitability_product_ndarray = np.prod(suitability_ndarray, axis=0, dtype=np.double)
    suitability_product = suitability_product_ndarray
    for p in range(P):
        for k in range(K): 
            for m in range(M):
                combination_weight[p, k, m] = (suitability[p, k, m] * suitability[p, k, m] / suitability_product[k, m]) ** gamma
    
    # Calculate spectrogram as a binwise weighted arithmetic mean.
    return np.average(X_ndarray, axis=0, weights=combination_weight_ndarray)

    ############ Spectrograms combination }}}