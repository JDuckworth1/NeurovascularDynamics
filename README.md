# NeurovascularDynamics

Code for processing wide-field and two photon vascular imaging data. Extract time series data from widefield and two photon tif images and perform spectral calculations to study neural and vascular dynamics.

Some files in this repository use functions from the [CHRONUX toolbox](http://chronux.org) and/or functions developed by Xiang Ji: [github](https://github.com/xiangjiph/VCRA)

## Wide-field
Extract and save df/f (x,t) from tif images. Compute the peak vasomotor frequency, preform space-frequency SVD, and plot the dominant mode's magnitude and phase. Extract phase at each location along the vasculature, calculate phase gradients and traveling wave properties.

## Penetrating Ves
Process interleaved two photon frame scan tif images. Estimate cross sectional diameter using the [Threshold in Radon Space method](https://journals.sagepub.com/doi/10.1038/jcbfm.2014.67), and perform spectral calculations on the resulting diameter time series. Calculate spectra, phase, and coherence for each trial and plot these results.

## Statistics
Helper functions for various calculations. Perform regression through the origin and calculate the associated uncertainty in slope. Perform [line subtraction](https://direct.mit.edu/neco/article/13/4/717/6503/Sampling-Properties-of-the-Spectrum-and-Coherency) to detect and remove sinusoidal components from complex timeseries data. Calculate 2D cross correlation using MATLAB's fft2 to detect motion artifacts in two photon data. Plot segmentation of vessels on top of the associated mask or raw image for visualization purposes.
