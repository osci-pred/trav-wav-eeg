# trav-wav-eeg
Traveling Waves Analysis Tools for EEG


Phase Plane Fit (based onZhang et al. Neuron 2018 / Joshua Jacobs lab):

use as:

wav = tw_fitPlaneEEG(eeg, time, channelLabels, 'Frequency', [7 13]);

where eeg is the raw data as [chan x timepoints x trials]  
time is the corresponding time vector in seconds  
channelLabels is a cell array of channel labels corresponding to first input dimension
