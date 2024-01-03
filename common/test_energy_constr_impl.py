import warnings
import numpy as np
import scipy.io
from calibration_energy_constr import calibration_energy_constr
from plot_traces import plot_traces

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    dataset = scipy.io.loadmat("dataset.mat")

#Loads the data from the simulation
hbw0 = float(dataset["w12"]) * 0.5 / np.pi
fs = float(dataset["fs"])
print("Fs:", fs, "Hz")
flattop_time = float(dataset["t1"])/fs
decay_time = float(dataset["t2"])/fs

a_sim = complex(dataset["a"])
b_sim = complex(dataset["b"])
c_sim = complex(dataset["c"])
d_sim = complex(dataset["d"])

a_sven = complex(dataset["a_n"])
b_sven = complex(dataset["b_n"])
c_sven = complex(dataset["c_n"])
d_sven = complex(dataset["d_n"])

x = dataset["x"]
y = dataset["y"]

vforw_cmplx = dataset["FOR_complex"]
vrefl_cmplx = dataset["REFL_complex"]
probe_cmplx = dataset["CAV_complex"]

vforw_cmplx = np.reshape(vforw_cmplx, vforw_cmplx.shape[0]) 
vrefl_cmplx = np.reshape(vrefl_cmplx, vrefl_cmplx.shape[0]) 
probe_cmplx = np.reshape(probe_cmplx, probe_cmplx.shape[0]) 

# Executes the implementation shown in the paper
[a_andrea, b_andrea, c_andrea, d_andrea] = list(calibration_energy_constr(fs,
                                                                          probe_cmplx, 
                                                                          vforw_cmplx, 
                                                                          vrefl_cmplx, 
                                                                          flattop_time,
                                                                          decay_time))
abcd_andrea = [a_andrea, b_andrea, c_andrea, d_andrea]
abcd_sven = [a_sven, b_sven, c_sven, d_sven]

print("Sim/Sven/Andrea-energy based constr.|a:", a_sim, a_sven, a_andrea)
print("Sim/Sven/Andrea-energy based constr.|b:", b_sim, b_sven, b_andrea)
print("Sim/Sven/Andrea-energy based constr.|c:", c_sim, c_sven, c_andrea)
print("Sim/Sven/Andrea-energy based constr.|d:", d_sim, d_sven, d_andrea)
print("\n")
print("Andrea's algorithm diff wrt. ground%|a:", float(np.abs(a_sim - a_andrea)*100/np.abs(x)), "normalized to a+c")
print("Andrea's algorithm diff wrt. ground%|b:", float(np.abs(b_sim - b_andrea)*100/np.abs(y)), "normalized to b+d")
print("Andrea's algorithm diff wrt. ground%|c:", float(np.abs(c_sim - c_andrea)*100/np.abs(x)), "normalized to a+c")
print("Andrea's algorithm diff wrt. ground%|d:", float(np.abs(d_sim - d_andrea)*100/np.abs(y)), "normalized to b+d")
print("\n")
print("Sven's algorithm diff wrt. ground%|a:", float(np.abs(a_sim - a_sven)*100/np.abs(x)), "normalized to a+c")
print("Sven's algorithm diff wrt. ground%|b:", float(np.abs(b_sim - b_sven)*100/np.abs(y)), "normalized to b+d")
print("Sven's algorithm diff wrt. ground%|c:", float(np.abs(c_sim - c_sven)*100/np.abs(x)), "normalized to a+c")
print("Sven's algorithm diff wrt. ground%|d:", float(np.abs(d_sim - d_sven)*100/np.abs(y)), "normalized to b+d")

plot_traces("Andrea's calibration", fs, probe_cmplx, vforw_cmplx, vrefl_cmplx, abcd_andrea, det0=100, hbw0=hbw0)
plot_traces("Sven's calibration", fs, probe_cmplx, vforw_cmplx, vrefl_cmplx, abcd_sven, det0=100, hbw0=hbw0)

input("Press return..")
