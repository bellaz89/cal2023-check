# Scripts to verify the paper 'Calibration of superconducting radio-frequency cavity forward and reflected channels based on stored energy dynamics'

These scripts have the purpouse of checking the correctness of the algorithms contained in the paper [Calibration of superconducting radio-frequency cavity forward and reflected channels based on stored energy dynamics](./Calibration_of_superconducting_radio_frequency_cavity_forward_and_reflected_channels_based_on_stored_energy_dynamics.pdf)

## ! DO NOT DIVULGE !

<u>This repository is for the coauthors of the paper only and other member of MSK involved. Please do not share any part of it!</u>.

## Before starting

Make sure the following packages are installed:

- python3
- numpy
- scipy
- numba
- matplotlib

If you want to regenerate the simulation data, make sure to have *Matlab* with the *Control Toolbox* installed.

## test1

This test verifies if Andrea's implementation of [1] with $k_{add} = 1$ matches Sven's implementation. A simulation provided by Sven is used.

### Regenerate the simulation data

```bash
matlab -nodisplay -nosplash -nodesktop -r "run('FORW_REFL_simulatiom.m');exit;"

```

## test2

This test verifies if Andrea's implementation of [1] with $k_{add} = 2.69595 e^{j\frac{118.7687 \pi}{180}}$ matches Sven's implementation. A simulation provided by Sven is used.

### Regenerate the simulation data

```bash
matlab -nodisplay -nosplash -nodesktop -r "run('FORW_REFL_simulatiom.m');exit;"

```

## test3

This test verifies if Andrea's implementation of the *Energy Constrained* algorithm is capable of finding the calibration coefficients. A simulation provided by Sven is used.

### Regenerate the datafile

```shell
  matlab -nodisplay -nosplash -nodesktop -r "run('FORW_REFL_simulation.m'); exit;"

```

[\[1\] Pfeiffer, Sven, et al. "Virtual cavity probe generation using calibrated forward and reflected signals." MOPWA040, These Proceedings, IPAC 15 \(2015\).](https://accelconf.web.cern.ch/IPAC2015/papers/mopwa040.pdf)

