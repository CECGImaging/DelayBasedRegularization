# Delay-Based Regularization for ECG Imaging of Transmembrane Voltages

<https://github.com/CECGImaging/DelayBasedRegularization>

This is an implementation of the delay-based regularization presented at the Computing in Cardiology conference 2019, [see paper](https://github.com/CECGImaging/DelayBasedRegularization/blob/master/CinC2019_DelayBasedRegularization.pdf).

## Running the examples ##

Two exemplary datasets are provided:

- A focal excitation on a spherical geometry simulated using an isotropic monodomain model and the Courtemanche et al. ionic model.
- Three pacings on a realistic ventricular geometry simulated using an anisotropic bidomain model and the ten Tusscher et al. ionic model.

To perform inverse reconstructions for these datasets, go to `MATLAB/delayRegu` and run

- `delayRegu_run_sphere.m` or
- `delayRegu_run_ventricles.m`, respectively.

## License ##

All source code is subject to the terms of the Mozilla Public License, v. 2.0.  
Copyright 2019 Steffen Schuler, Karlsruhe Institute of Technology.

## Contact ##

Steffen Schuler  
Institute of Biomedical Engineering  
Karlsruhe Institute of Technology  
www.ibt.kit.edu
