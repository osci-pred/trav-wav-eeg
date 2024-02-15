function out = tw_fitPlaneEEG(raw, time, labels, varargin)
% raw                   eeg input; [nChan x nTime x nTrials];
% time                  time vector in seconds
% labels                cell array of electrode labels, conforming to 10/05
%                       labels (using fieldtrip layout)
%
% Optional Input (as named arg list):
%
% 'Frequency'           temporal freq limits [fmin fmax] used for bp filter
% 'ROI'                 defines limits of electrode region to be used for the fit 
%                       (see fieldtrip layout EEG1005)
% 'WindowSize'          temporal window size (ms) used to average relative phase 
%                       (leave empty if single-point)
% 'MaxCycles'           max. spatial frequency, defined as number of cycles
%                       within the covered electrode space, default 1 cycle
% 'NumStepsSpatFreq'    number of steps for spatial freq [0 MaxCycles],
%                       default 30
% 'NumStepsWaveDir'     number of steps for wave direction [-pi pi],
%                       default 60

p = inputParser();
p.addParameter('Frequency', [7 13]); 
p.addParameter('ROI', {[-0.2 0.2] [-Inf 0.25]});
p.addParameter('WindowSize', []); 
p.addParameter('MaxCycles', 1); 
p.addParameter('NumStepsSpatFreq', 30); 
p.addParameter('NumStepsWaveDir', 60); 

p.parse(varargin{:});


raw = double(raw);
sr = round(1/diff(time(1:2)));

%%

% Get the electrode positions that we want to include in the fit:
[pos, lbl, idxInData] = tw_getElecPos(labels, 'EEG1005', p.Results.ROI);


% Get the indices for the midline electrodes (not used in the fit):
[midLine, midLineLabel, midLineIdxInData] = tw_getMidLine(lbl, idxInData);

% Filter the data and extract phases:
[phi, pow, eegFilt] = tw_getBandpassPhase(raw(idxInData,:,:), sr, p.Results.Frequency);

% Get the phase planes predicted by each fit:
[phiPred, wvdir, a, b, xi] = tw_getPlaneFits(pos,...
    p.Results.NumStepsSpatFreq, p.Results.NumStepsSpatFreq, p.Results.MaxCycles);
% -> first dimension for each is the number of different fits

% Evaluate the fits:
[id, rcc_sq, phi_out] = tw_evalPlaneFits(phi, phiPred, (p.Results.WindowSize/1e3) .* sr);
 


%%
out.t = time;
out.lbl = lbl;
out.pos = pos;
out.xi = xi(id);
out.a = a(id);
out.b = b(id);
out.wavDir = wvdir(id);
out.rcc_sq = rcc_sq;
out.midLineIdx = midLine;
out.midLineIdxInData = midLineIdxInData;
out.midLineLabel = midLineLabel;
out.pow = pow;
out.phi = phi_out;

%%
end


