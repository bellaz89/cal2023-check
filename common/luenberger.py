import numpy as np
import numba
from numba import njit

#@njit(fastmath=True)
def luenberger(fwd, prb, hbw, detuning,
             fs, gain_bw, a_threshold):

    bandwidth = hbw * 2
    max_ampsq = np.max(prb[:, 0]**2 + prb[:, 1])
    div_lim = a_threshold * max_ampsq

    samples = fwd.shape[0]
    K = np.zeros((samples, 4, 2))
    x = np.zeros((samples, 4))
    xpred = np.zeros((samples, 4))

    x[0, 2] = -0.5
    x[0, 3] = detuning/bandwidth

    A = np.zeros((samples, 4, 4))

    bandwidth_angle = 2.0 * np.pi * bandwidth
    gain_angle = 2.0 * np.pi * gain_bw

    rot_exp = np.exp(-bandwidth_angle / fs)
    gain_exp = np.exp(-gain_angle / fs)

    A[:, 0, 0] = rot_exp
    A[:, 1, 1] = rot_exp
    A[:, 2, 2] = 1.0
    A[:, 3, 3] = 1.0

    rot_diff = 1-rot_exp
    gain_diff = 1-gain_exp

    B = np.zeros((4, 2))
    B[0, 0] = rot_diff
    B[1, 1] = rot_diff

    C = np.zeros((2, 4))
    C[0, 0] = 1.0
    C[1, 1] = 1.0

    K[:, 0, 0] = -gain_diff + 2.0 * gain_diff
    K[:, 1, 1] = -gain_diff + 2.0 * gain_diff
    K[:, 2:, :] = gain_diff**2 / rot_diff
    K[0, 2:, :] = 0

    for i in range(1, samples):

        # pag 123. S. M. Bozic

        # prediction

        A[i, 0, 2] = - x[i-1, 0] * rot_diff
        A[i, 0, 3] =   x[i-1, 1] * rot_diff
        A[i, 1, 2] = - x[i-1, 1] * rot_diff
        A[i, 1, 3] = - x[i-1, 0] * rot_diff
        
        xpred[i, :] = np.dot(A[i, :, :], x[i-1, :]) + np.dot(B, fwd[i-1, :]) 

        # updating

        x[i, :] = xpred[i, :] + np.dot(K[i-1, :, :], prb[i, :] - xpred[i, :2])

        amp = x[i, 0]**2 + x[i, 1]**2

        #a_recip = 1.0 / div_lim
        a_recip = 0.0

        if div_lim < amp:
            a_recip = 1.0 / amp

        K[i, 2:, :] *= a_recip
        K[i, 2, 0] *= -x[i, 0]
        K[i, 2, 1] *= -x[i, 1]
        K[i, 3, 0] *=  x[i, 1]
        K[i, 3, 1] *= -x[i, 0]

    x[:, 2] = bandwidth * (1.0 + x[:, 2])
    x[:, 3] *= bandwidth
    return x

