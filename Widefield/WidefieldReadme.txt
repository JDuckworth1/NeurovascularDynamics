Widefieldreadme.txt

Code to extract df/f(x,t) from individually saved tif images. 
Additional code includes processing (defining vasomotor peak, saving df/f(x,t)) and plotting
(computing the space-freq SVD and mapping magnitude/phase/additional statistics back to the vessel mask).
Also included is code used to create the vessel graph (all "nodes" and "links") from pial vessel masks, and 
calculate df/f phase progression along each vessel link.

load_process_rawdata.m loads tifs and calculates df/f(x,t) for neuronal (Thy1-jRGECO) and SMC GCaMP (GCaMP8.1) fluorescent images.
widefield_VesselSpacefreqSVD.m and widefield_PlotVesselSpaceFreqSVD.m calculats and plots the space-freq SVD magnitude and phase.

widefield_ExtractVesNeuPhaseGrad.m loads each trial's neuronal and vascular df/f, computes the vessel graph, and calculates 
the phase of each signal along each vessel. The phase gradient k is then calculated by fitting phase vs. distance along each vessel. Other calculations such as vessel-neuron coherence are also included.

widefield_CombineVesNeuPhaseGrad.m combines results from widefield_ExtractVesNeuPhaseGrad.m across all animals, trials, and vessels. These combined results are then plotted as summary figures.