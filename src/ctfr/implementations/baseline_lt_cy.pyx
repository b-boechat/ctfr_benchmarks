import numpy as np
cimport cython
from libc.math cimport sqrt, pow
from warnings import warn
from ctfr.warning import ArgumentChangeWarning

def _baseline_lt_wrapper(X, lk = 21, lm = 11, eta = 8.0):
    lk = int(lk)
    if lk < 0:
        lk = 21
        warn(f"The 'lk' parameter should be a positive integer. Setting lk = {lk}.", ArgumentChangeWarning)
    if lk % 2 == 0:
        lk += 1
        warn(f"The 'lk' parameter should be an odd integer. Setting lk = {lk}.", ArgumentChangeWarning)

    lm = int(lm)
    if lm < 0:
        lm = 11
        warn(f"The 'lm' parameter should be a positive integer. Setting lm = {lm}.", ArgumentChangeWarning)
    if lm % 2 == 0:
        lm += 1
        warn(f"The 'lm' parameter should be an odd integer. Setting lm = {lm}.", ArgumentChangeWarning)
    
    eta = float(eta)
    if eta < 0.0:
        eta = 8.0
        warn(f"The 'eta' parameter should be a non-negative float. Setting eta = {eta}.", ArgumentChangeWarning)

    return _baseline_lt_cy(X, lk, lm, eta)

@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.nonecheck(False)
@cython.cdivision(True)
cdef _baseline_lt_cy(double[:,:,::1] X_orig, Py_ssize_t lk, Py_ssize_t lm, double eta):

    cdef:
        Py_ssize_t P = X_orig.shape[0]
        Py_ssize_t K = X_orig.shape[1]
        Py_ssize_t M = X_orig.shape[2]
        
        Py_ssize_t lk_lobe = (lk-1)//2
        Py_ssize_t lm_lobe = (lm-1)//2
        Py_ssize_t p, m, k, i, j

        double epsilon = 1e-15
        Py_ssize_t combined_size = lm * lk

    # Zero-pad spectrograms for windowing
    X_ndarray = np.pad(X_orig, ((0, 0), (lk_lobe, lk_lobe), (lm_lobe, lm_lobe)))
    cdef double[:, :, :] X = X_ndarray

    # Container that stores an horizontal segment of a spectrogram, with all frequency bins. Used to calculate smearing.
    calc_vector_ndarray = np.zeros(combined_size, dtype = np.double)
    cdef double[:] calc_vector = calc_vector_ndarray 
 
    # Container that stores the result.
    result_ndarray = np.zeros((K, M), dtype=np.double)
    cdef double[:, :] result = result_ndarray

    # Container that stores the local smearing.
    smearing_ndarray = np.zeros((P, K, M), dtype=np.double)
    cdef double[:,:,:] smearing = smearing_ndarray

    # Variables related to smearing calculation.
    cdef double smearing_numerator, smearing_denominator

    # Variables related to spectrogram combination.
    cdef double weight, weights_sum, result_acc

    ############ Local smearing calculation {{{

    for p in range(P):
        # Iterates through neighborhoods.
        for k in range(lk_lobe, K + lk_lobe):
            for m in range(lm_lobe, M + lm_lobe):

                # Obtain the neighborhood sorted vector.
                for i in range(lk):
                    for j in range(lm):
                        calc_vector[i*lm + j] = X[p, k - lk_lobe + i, m - lm_lobe + j]                                
                calc_vector_ndarray.sort()

                # Compute the smearing function.
                smearing_denominator = 0.0
                smearing_numerator = 0.0
                for i in range(combined_size):
                    smearing_denominator = smearing_denominator + calc_vector[o]
                    smearing_numerator = smearing_numerator + (combined_size-o)*calc_vector[o]
                smearing[p, k - lk_lobe, m - lm_lobe] = smearing_numerator/(sqrt(smearing_denominator) + epsilon)
                
    ############ }}}

    ############ Spectrograms weighted combination {{{

    for k in range(K):
        for m in range(M):
            weights_sum = 0.0
            result_acc = 0.0
            for p in range(P):
                weight = 1./(pow(smearing[p, k, m], eta) + epsilon)
                result_acc = result_acc + weight * X_orig[p, k, m]
                weights_sum = weights_sum + weight
            result[k, m] = result_acc / weights_sum

    ############ }}}

    return result_ndarray