# NeurovascularDynamics

Extract time series data from widefield and two photon tif images and perform spectral calculations to study neural and vascular dynamics.

Some files in this repository use functions from the [CHRONUX toolbox](http://chronux.org) and/or functions developed by Xiang Ji: [github](https://github.com/xiangjiph/VCRA)

## Wide-field
Extract and save df/f (x,t) from tif images. Compute the peak vasomotor frequency, preform space-frequency SVD, and plot the dominant mode's magnitude and phase. Extract phase at each location along the vasculature, calculate phase gradients and traveling wave properties.

## Penetrating Ves
Process interleaved two photon frame scan tif images. Estimate cross sectional diameter using the [Threshold in Radon Space method](https://journals.sagepub.com/doi/10.1038/jcbfm.2014.67), and perform spectral calculations on the resulting diameter time series. Calculate spectra, phase, and coherence for each trial and plot these results.

## Statistics
Helper functions for various calculations. Perform regression through the origin and calculate the associated uncertainty in slope. Perform [line subtraction](https://direct.mit.edu/neco/article/13/4/717/6503/Sampling-Properties-of-the-Spectrum-and-Coherency) to detect and remove sinusoidal components from complex timeseries data.

## Pial2P
Process and analyze two photon frame scans with two simultaneous imaging channels: Cy5.5 and smooth muscle GCaMP8.1. Calculate GCaMP df/f(x,t), integrated Cy5.5 intensity(x,t), and lumen diameter(x,t) via Full-Width-(Scale factor)-Max. Analyze the relations between these signals using cross correlation and spectral methods.

## VesGraph
Functions required to compute the vessel graph (nodes and links) from the pial vasculature mask. Written by Xiang Ji: [github](https://github.com/xiangjiph/VCRA)

## RBCFlux
Process two photon line scan images with two simultaneous imaging channels: Cy5.5 annd RBCs labeled with CellTrace CFSE. Detect and count RBCs, and calculate vessel diameter using FWHM. Compute flux and diameter as a function of time, and dFlux/dDiameter.

