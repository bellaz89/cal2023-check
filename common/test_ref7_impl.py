import warnings
import numpy as np
import scipy.io
from calibration_impl import calibrate_from_ref_7_kadd_1

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    dataset = scipy.io.loadmat("dataset.mat")

#Loads the data from the simulation

kadd = complex(dataset["k_add"])

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
vforw_cmplx_decay = dataset["FOR_decay_compl"]
vrefl_cmplx_decay = dataset["REFL_decay_compl"]
probe_cmplx_decay = dataset["CAV_decay_compl"] 

vforw_cmplx = np.reshape(vforw_cmplx, vforw_cmplx.shape[0]) 
vrefl_cmplx = np.reshape(vrefl_cmplx, vrefl_cmplx.shape[0]) 
probe_cmplx = np.reshape(probe_cmplx, probe_cmplx.shape[0]) 
vforw_cmplx_decay = np.reshape(vforw_cmplx_decay, vforw_cmplx_decay.shape[0])
vrefl_cmplx_decay = np.reshape(vrefl_cmplx_decay, vrefl_cmplx_decay.shape[0])
probe_cmplx_decay = np.reshape(probe_cmplx_decay, probe_cmplx_decay.shape[0])

# Executes the implementation shown in the paper
[a_andrea, b_andrea, c_andrea, d_andrea] = list(calibrate_from_ref_7_kadd_1(probe_cmplx, 
                                                                            vforw_cmplx, 
                                                                            vrefl_cmplx, 
                                                                            probe_cmplx_decay, 
                                                                            vforw_cmplx_decay, 
                                                                            vrefl_cmplx_decay, 
                                                                            kadd))

print("Sim/Sven/Andrea|a:", a_sim, a_sven, a_andrea)
print("Sim/Sven/Andrea|b:", b_sim, b_sven, b_andrea)
print("Sim/Sven/Andrea|c:", c_sim, c_sven, c_andrea)
print("Sim/Sven/Andrea|d:", d_sim, d_sven, d_andrea)
print("\n")
print("Algorithm diff %|a:", float(np.abs(a_sven - a_andrea)*100/np.abs(x)), "normalized to a+c")
print("Algorithm diff %|b:", float(np.abs(b_sven - b_andrea)*100/np.abs(y)), "normalized to b+d")
print("Algorithm diff %|c:", float(np.abs(c_sven - c_andrea)*100/np.abs(x)), "normalized to a+c")
print("Algorithm diff %|d:", float(np.abs(d_sven - d_andrea)*100/np.abs(y)), "normalized to b+d")
