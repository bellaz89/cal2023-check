import numpy as np
from calibration_energy_constr import calibration_energy_constr, find_hbw
from calibration_impl import calibrate_diagonal
from plot_traces import plot_traces
import sys

SET = sys.argv[1]

filename = "data{}.npz".format(SET)
dataset = np.load(filename)
print("Opening file \"", filename, "\"")

#Loads the data from the simulation
fs = dataset["fs"]*1e6
delay_time = dataset["delay"]*1e-6
flattop_time = dataset["filling"]*1e-6 + delay_time
decay_time = dataset["flattop"]*1e-6 + flattop_time

probe_cmplx = dataset["probe_amp"]*np.exp(1.0j*dataset["probe_pha"]*np.pi/180)
vforw_cmplx = dataset["forward_amp"]*np.exp(1.0j*dataset["forward_pha"]*np.pi/180)
vrefl_cmplx = dataset["reflected_amp"]*np.exp(1.0j*dataset["reflected_pha"]*np.pi/180)

[x, _, _, y] = calibrate_diagonal(probe_cmplx, vforw_cmplx, vrefl_cmplx)
print("x", x, "y", y)
vforw_cmplx *= x
vrefl_cmplx *= y

print("Fs:", fs, "Hz")
print("delay:", delay_time, "s")
print("flattop:", flattop_time, "s")
print("decay:", decay_time, "s")
print("pulse length:", probe_cmplx.shape[0]/fs, "s")

time_trace = np.linspace(0, probe_cmplx.shape[0]/fs, probe_cmplx.shape[0])
probe_amp_decay = dataset["probe_amp"][time_trace > decay_time]
time_trace_decay = time_trace[time_trace > decay_time]

hbw0 = find_hbw(time_trace_decay, probe_amp_decay) * 0.5 / np.pi
print("hbw0:", hbw0, "Hz")

# Executes the implementation shown in the paper
[a_andrea, b_andrea, c_andrea, d_andrea] = list(calibration_energy_constr(fs,
    probe_cmplx,
    vforw_cmplx,
    vrefl_cmplx,
    flattop_time,
    decay_time,
    start_time=delay_time))

abcd_andrea = [a_andrea, b_andrea, c_andrea, d_andrea]

print("Andrea-energy based constr.|a:", a_andrea)
print("Andrea-energy based constr.|b:", b_andrea)
print("Andrea-energy based constr.|c:", c_andrea)
print("Andrea-energy based constr.|d:", d_andrea)
print("\n")

plot_traces("Hobicat calibration", fs, probe_cmplx, vforw_cmplx, vrefl_cmplx, abcd_andrea, det0=100, hbw0=hbw0, filter_bw=20000)

input("Press return..")
