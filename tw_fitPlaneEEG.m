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
% 'RandShuffleIter'     number of iterations for random
%                       shuffling to estimate the chance level of
%                       rcc_square.
% 'RandShuffleNTrials'  number of trials to perform random shuffling on. If
%                       empty, all trials
% 'DirClassWindowSize'  window size for classification of wave direction (in radians)
% 'RandShuffleTimeStep' time step option for shuffling, in ms. 


p = inputParser();
p.addParameter('Frequency', [7 13]);
p.addParameter('ROI', {[-0.2 0.2] [-Inf 0.25]});
p.addParameter('WindowSize', []);
p.addParameter('MaxCycles', 1);
p.addParameter('NumStepsSpatFreq', 30);
p.addParameter('NumStepsWaveDir', 60);
p.addParameter('RandShuffleIter', 10);
p.addParameter('RandShuffleNTrials', []);
p.addParameter('RandShuffleTimeStep', 200); % in ms
p.addParameter('DirClassWindowSize', 0.5);

p.parse(varargin{:});

%%
raw = double(raw);
nTr = size(raw,3);
nTm = size(raw,2);
sr = round(1/diff(time(1:2)));

% convert moving-avg windows size to samples:
if ~isempty(p.Results.WindowSize)
    movMeanWinSize = (p.Results.WindowSize/1e3) .* sr;
else
    movMeanWinSize = [];
end

if isempty(p.Results.RandShuffleNTrials)
    shuffleTrials = 1:nTr;
else
    shuffleTrials = sort(randsample(1:size(raw,3),...
        min(nTr, p.Results.RandShuffleNTrials)));
end

if ~isempty(p.Results.RandShuffleTimeStep)
    shStep = round((p.Results.RandShuffleTimeStep/1e3) .* sr);
    
    shTm = shStep:shStep:nTm;
else
    shTm = 1:nTm;
end

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
[id, rcc_sq, rcc_sq_rand, phi_out] = tw_evalPlaneFits(phi, phiPred, movMeanWinSize,...
    p.Results.RandShuffleIter, shuffleTrials, shTm);

[fw, bw] = tw_classifyDirection(wvdir(id), p.Results.DirClassWindowSize);


%%
out.t = time;
out.fw = fw;
out.bw = bw;

out.lbl = lbl;
out.pos = pos;
out.xi = xi(id);
out.a = a(id);
out.b = b(id);
out.wavDir = wvdir(id);

out.rcc_sq = rcc_sq;
out.rcc_sq_rand = rcc_sq_rand(:,shuffleTrials,:);
out.rcc_sq_rand = out.rcc_sq_rand(shTm, :,:);
out.rand_trials = shuffleTrials;
out.rand_t = out.t(shTm);
out.rcc_thresh = prctile(out.rcc_sq_rand(:), 99);
out.sig = out.rcc_sq > out.rcc_thresh;

out.pFitFW = mean(out.sig & out.fw,2);
out.pFitBW = mean(out.sig & out.bw,2);

out.midLineIdx = midLine;
out.midLineIdxInData = midLineIdxInData;
out.midLineLabel = midLineLabel;
out.pow = pow;
out.phi = phi_out;

%%
end


