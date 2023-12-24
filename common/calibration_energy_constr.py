import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from calibration_impl import calibrate_energy

ORDER_SAVGOL = 3

#Find the half bandwidth using 
def find_hbw(time_decay, probe_cmplx_decay):
    probe_amp = np.abs(probe_cmplx_decay)

    def decay(t, a, b):
        return a * np.exp(-(t-time_decay[0]) * b)

    popt, pcov = curve_fit(decay,
                           time_decay, probe_amp,
                           p0=[probe_amp[0], np.pi * 130])

    # print("Found cavity half bandwidth is", popt[1]*0.5/np.pi, "Hz")

    return popt[1]


# Calibrate using the energy constrained method
def calibration_energy_constr(fs, probe_cmplx, vforw_cmplx, vrefl_cmplx, 
                              flattop_time, decay_time, transition_guard=100, 
                              start_time=0.0, stop_time=None):

    # Time step
    dt = 1.0/fs

    # Time around the transitions to ignore
    time_guard = transition_guard * dt

    # Filtering steps
    savgol_range = transition_guard + 1 - transition_guard % 2

    trace_length = probe_cmplx.shape[0]
    time_trace = np.linspace(0, trace_length*dt, trace_length)

    if stop_time == None:
        stop_time = dt * trace_length

    # Derivative of the probe squared amplitude
    probe_sq_deriv = savgol_filter(np.abs(probe_cmplx)**2,
                                       savgol_range, ORDER_SAVGOL,
                                       deriv=1,
                                       delta=dt)

    # Indices arrays for filling flattop and decay
    filling_idxs = ((time_trace >= start_time + time_guard) * (time_trace < flattop_time - time_guard))
    flattop_idxs = ((time_trace >= flattop_time + time_guard) * (time_trace < decay_time - time_guard))
    decay_idxs = ((time_trace >= decay_time + time_guard) * (time_trace < stop_time - time_guard)).astype(bool)
    pulse_idxs = (filling_idxs + flattop_idxs + decay_idxs).astype(bool)

    time_trace_decay = time_trace[decay_idxs]
    probe_cmplx_decay = probe_cmplx[decay_idxs]
    vforw_cmplx_decay = vforw_cmplx[decay_idxs]
    vrefl_cmplx_decay = vrefl_cmplx[decay_idxs] 

    time_trace_pulse = time_trace[pulse_idxs]
    probe_cmplx_pulse = probe_cmplx[pulse_idxs]
    vforw_cmplx_pulse = vforw_cmplx[pulse_idxs]
    vrefl_cmplx_pulse = vrefl_cmplx[pulse_idxs] 
    probe_sq_deriv_pulse = probe_sq_deriv[pulse_idxs]

    #Half bandwidth (in angular frequency) computed on the probe decay
    hbw_decay = find_hbw(time_trace_decay, probe_cmplx_decay)

    return calibrate_energy(hbw_decay,
                            probe_cmplx_pulse, vforw_cmplx_pulse, vrefl_cmplx_pulse,
                            probe_sq_deriv_pulse,
                            probe_cmplx_decay, vforw_cmplx_decay, vrefl_cmplx_decay)

