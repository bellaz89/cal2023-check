# Scripts to verify the paper 'Calibration of superconducting radio-frequency cavity forward and reflected channels based on stored energy dynamics'

These scripts have the purpouse of checking the correctness of the algorithms contained in the paper [Calibration of superconducting radio-frequency cavity forward and reflected channels based on stored energy dynamics](./Calibration_of_superconducting_radio_frequency_cavity_forward_and_reflected_channels_based_on_stored_energy_dynamics.pdf)

## Before starting

Make sure the following packages are installed:

- python3
- numpy
- scipy
- numba
- matplotlib

If you want to regenerate the simulation data, make sure to have *Matlab* installed.

## test1

This test verifies if Andrea's implementation of [1] with $k_{add} = 1$  Matches Sven's implementation.

### Regenerate the datafile

```shell
  matlab -nodisplay -nosplash -nodesktop -r "run('FORW_REFL_SIMULATION'); exit;"
```

## test2

This test verifies if Andrea's implementation of [1] with $k_{add} = 2.69595 e^{j\frac{118.7687 \pi}{180}}$  Matches Sven's implementation.


### Regenerate the datafile

```shell
  matlab -nodisplay -nosplash -nodesktop -r "run('FORW_REFL_SIMULATION'); exit;"
```

[\[1\] Pfeiffer, Sven, et al. "Virtual cavity probe generation using calibrated forward and reflected signals." MOPWA040, These Proceedings, IPAC 15 \(2015\).](https://accelconf.web.cern.ch/IPAC2015/papers/mopwa040.pdf)

