import numpy as np
import matplotlib.pyplot as plt
from luenberger import luenberger

plt.ion()

def correct_signals(abcd, vforw, vrefl):
    return (abcd[0] * vforw + abcd[1] * vrefl, abcd[2] * vforw + abcd[3] * vrefl)


def plot_traces(title, fs, probe_cmplx, vforw_cmplx, vrefl_cmplx, abcd, filter_bw=100000, hbw0=0, det0=0):
    dt = 1.0 / fs
    fig, (ax_amp, ax_hbwdet) = plt.subplots(2, 1, sharex=True)
    fig.suptitle(title)
    
    ax_amp.set_xlabel("Time (s)")
    ax_amp.set_ylabel("Amplitude (MV/m)")
    ax_hbwdet.set_xlabel("Time (s)")
    ax_hbwdet.set_ylabel("Detuning/Half bandwidth (Hz)")

    (vforw_cmplx_corr, vrefl_cmplx_corr) = correct_signals(abcd, vforw_cmplx, vrefl_cmplx)

    time_trace = np.linspace(0, probe_cmplx.shape[0]*dt, probe_cmplx.shape[0])

    ax_amp.plot(time_trace, np.abs(probe_cmplx), label="Probe")
    ax_amp.plot(time_trace, np.abs(vforw_cmplx), label="Forward")
    ax_amp.plot(time_trace, np.abs(vforw_cmplx_corr), label="Forward calibrated")

    ax_amp.legend()

    probe = np.column_stack((np.real(probe_cmplx), np.imag(probe_cmplx)))
    vforw = np.column_stack((np.real(vforw_cmplx), np.imag(vforw_cmplx)))
    vforw_corr = np.column_stack((np.real(vforw_cmplx_corr), np.imag(vforw_cmplx_corr)))

    x = luenberger(vforw, probe, hbw0, det0, fs, filter_bw, 0.1)
    x_corr = luenberger(vforw_corr, probe, hbw0, det0, fs, filter_bw, 0.01)

    hbw = x[:, 2]
    det = x[:, 3]

    hbw_corr = x_corr[:, 2]
    det_corr = x_corr[:, 3]

    ax_hbwdet.plot(time_trace, det, label="Detuning")
    ax_hbwdet.plot(time_trace, det_corr, label="Detuning calibrated")

    ax_hbwdet.plot(time_trace, hbw, label="Bandwidth")
    ax_hbwdet.plot(time_trace, hbw_corr, label="Bandwidth calibrated")

    plt.legend()
    plt.show()
