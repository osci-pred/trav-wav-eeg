function out = tw_fitPlaneEEGWavelet(raw, time, labels, varargin)
% Performs the iterative plane fitting procedure based on a continuous
% wavelet transform.
%
% raw                   eeg input; [nChan x nTime x nTrials];
% time                  time vector in seconds
% labels                cell array of electrode labels, conforming to 10/05
%                       labels (using fieldtrip layout)
%
% Optional Input (as named arg list):
%
% 'WaveletArgs'         optional argument list to cwt function passed as cell array. 
% 'ROI'                 defines electrode ROI to be used for the fit
%                       (see fieldtrip layout EEG1005). Can be either in
%                       coordinates as {[ml_min ml_max] [ap_min ap_max]},
%                       or as a simple cell array of electrodes to include
% 'WindowSize'          temporal window size (in osc cycles) used to average relative phase
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
% 'ReturnPowChan'       (optional) list of electrodes to return power
%                       spectra for individually (in addition to the mean power across the ROI)
% 'ReturnWTOnly'        if true, returns after computing the wavelet
%                       transform


p = inputParser();
p.addParameter('WaveletArgs', {'amor', 'FrequencyLimits', [2 30]});
p.addParameter('ROI', {[-0.2 0.2] [-Inf 0.25]});
p.addParameter('WindowSize', 1);
p.addParameter('MaxCycles', 1);
p.addParameter('NumStepsSpatFreq', 30);
p.addParameter('NumStepsWaveDir', 30);
p.addParameter('RandShuffleIter', 10);
p.addParameter('RandShuffleNTrials', []);
p.addParameter('RandShuffleTimeStep', 200); % in ms
p.addParameter('DirClassWindowSize', 0.5);
p.addParameter('ReturnPowChan', {});
p.addParameter('ReturnWTOnly', false);
p.parse(varargin{:});

%%

raw = double(raw);
nTr = size(raw,3);
nTm = size(raw,2);
sr = round(1/diff(time(1:2)));

% check for nans in the data and throw back to user:
assert(~any(isnan(raw(:))), 'Data contains NaNs, cannot run CWT');


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
[pos, lbl, idxInData, powIdx] = tw_getElecPos(labels, 'EEG1005', p.Results.ROI, p.Results.ReturnPowChan);


% Get the indices for the midline electrodes (not used in the fit):
[midLine, midLineLabel, midLineIdxInData] = tw_getMidLine(lbl, idxInData);

% Filter the data and extract phases:
[phi, freq, avgPow, chanPow] = getWT(raw(idxInData,:,:), sr, p.Results.WaveletArgs, powIdx);

if p.Results.ReturnWTOnly
    % Optional return after WT computation (to do initial processing before fitting or re-analysis splits)
    
    [wvdir, a, b, xi, id, rcc_sq] = deal([]);
    
else
    
    % Convert moving mean window sizes to samples using freq vector:
    movMeanWinSize = p.Results.WindowSize .* (sr./freq);
    
    
    % Get the phase planes predicted by each fit:
    [phiPred, wvdir, a, b, xi] = tw_getPlaneFits(pos, ...
        p.Results.NumStepsWaveDir, p.Results.NumStepsSpatFreq, p.Results.MaxCycles);
    % -> first dimension for each is the number of different fits
    
    % Evaluate the fits:
    [id, rcc_sq] = evalPlaneFitsWT(phi, phiPred, movMeanWinSize,...
        p.Results.RandShuffleIter, shuffleTrials, shTm);
    
end

% [fw, bw] = tw_classifyDirection(wvdir(id), p.Results.DirClassWindowSize);


%%
out.t = time;
out.f = freq;
% out.fw = fw;
% out.bw = bw;

out.lbl = lbl;
out.pos = pos;
out.xi = xi(id);
out.a = a(id);
out.b = b(id);
out.wavDir = wvdir(id);

out.rcc_sq = rcc_sq;

% out.rcc_sq_rand = rcc_sq_rand(:,shuffleTrials,:);
% out.rcc_sq_rand = out.rcc_sq_rand(shTm, :,:);
% out.rand_trials = shuffleTrials;
% out.rand_t = out.t(shTm);
% out.rcc_thresh = prctile(out.rcc_sq_rand(:), 99);
% out.sig = out.rcc_sq > out.rcc_thresh;

% out.pFitFW = mean(out.sig & out.fw,2);
% out.pFitBW = mean(out.sig & out.bw,2);

out.midLineIdx = midLine;
out.midLineIdxInData = midLineIdxInData;
out.midLineLabel = midLineLabel;
out.avgPow = avgPow;
out.chanPow = permute(chanPow, [2 3 4 1]);
out.chanPowLbl = lbl(powIdx);

out.param = p.Results;
%%
end
%%
function [phi, freq, avgPow, chanPow] = getWT(raw, sr, wtArgs, powIdx)
% return the CWT
% raw is raw eeg signal as [nChan x nTime x nTrials]
% wtArgs is a cell array of optional additional arguments to the cwt
% function

%%
nTr = size(raw,3);
nTm = size(raw,2);
nElec = size(raw,1);

wtfn = @(x) cwt(x, sr, wtArgs{:}); % shorthand the call to cwt with added arguments


% run once to preallocate output matrix:
[~, freq] = wtfn(raw(1,:,1));

% we flip our freq dimension here and in the output:
freq = flip(freq);

outwt = nan(nElec, numel(freq), nTm, nTr);

wb = waitbar(0, 'Computing continuous wavelet transform...');

for itrial = 1:nTr
    waitbar(itrial/nTr, wb, sprintf('Computing continuous wavelet transform: Trial %d / %d', itrial, nTr));
    for iElec = 1:nElec
        outwt(iElec,:,:,itrial) = wtfn(raw(iElec,:,itrial));
    end
end
delete(wb);

outwt = outwt(:,end:-1:1, :, :); % flip frequency dimension


% return magnitude/power:
avgPow = squeeze(mean(abs(outwt),1));
chanPow = abs(outwt(powIdx, :, :, :));


phi = outwt./abs(outwt);
mphi = mean(phi,1)./abs(mean(phi,1));
phi = phi./mphi; %  correct phase by mean phase across elecs per timepoint
% -> now we can average over time (within short windows)


end
%%
function [id, rcc_sq] = evalPlaneFitsWT(phi, phiPred, movMeanWinSize, nIter, shuffleTrials, shTm)
% compares predicted phases in phiPred to actual phases in phi, after
% applying moving-mean of window size movMeanWinSize (in samples, not ms);
% phases in phi must be relative (e.g. to mean phase across electrode) if
% moving mean is requested
% Returned id is the index of the best fit for each time-point, with the same
% [nTime x nTrial] dimensions as input

%%
[nFreq, nTm, nTr] = size(phi, [2 3 4]);

%%

% if nargin < 4 || isempty(nIter)
%     nIter = 0;
% end

% %% Random shuffling:
% if size(phi,1) >= 10
%     for iter = 1:nIter
%         phi_rand(:,:,:,iter) = phi(randperm(size(phi,1)),:,:);
%     end
% else
%     % for small channel counts (mostly modeling) we make sure that no two shufflings are identical
%     assert(nIter < factorial(size(phi,1)),...
%         'Number of iterations is greater than number of possible permutations')
%     rperm = perms(1:size(phi,1));
%     rperm = rperm(randsample(1:size(rperm,1), nIter),:);
%     
%     for iter = 1:nIter
%         phi_rand(:,:,:,iter) = phi(rperm(iter,:),:,:);
%     end
% end

%%
[id, rcc_sq] = deal(nan(nFreq,nTm,nTr));
% rcc_sq_rand = nan(nTm,nTr,nIter);
% phi_out = nan(size(phi,1), nTm, nTr);


phiPred = repmat(phiPred, 1,1,1, nTr); % expand predicted phases to freq dimension
phiPred = permute(phiPred, [2 3 4 1]); % move fit dimension to the end

wb = waitbar(0, {sprintf('Wave Fit: Timepoint %d/%d', 0,nTm), '~~~ >> ??? << ~~~'});

winHW = ceil(movMeanWinSize(1)/2);
movAvg = zeros(1,nFreq, 2*winHW+1);
movAvg(1,:,:) = (-winHW:winHW) > -movMeanWinSize./2 & (-winHW:winHW) < movMeanWinSize./2;

for iT = 1:nTm
    %     movAvg = zeros(1,nFreq, nTm);
    %     movAvg(1,:,:) = movAvgIdx(iT);
    
    waitbar(iT/nTm,wb, {sprintf('Wave Fit: Timepoint %d/%d', iT, nTm), '~~~ >> ??? << ~~~'});
    
    %     for itrial = 1:nTr
    tIdx = iT + (-winHW:winHW);
    tIdx = tIdx(tIdx >= 1 & tIdx <= nTm);
    
    aphi = squeeze(mean(phi(:,:,tIdx,:).*movAvg(:,:,tIdx-iT+winHW+1),3)); % -> elec x freq
    aphi = aphi./abs(aphi); % discard magnitude of average
    
    [id(:,iT,:), rcc_sq(:,iT,:)] = singleFitEval(aphi, phiPred);
    %         phi_out(:,iT,itrial) = aphi;
    %     end
    
    %     if ismember(itrial, shuffleTrials)
    %         for iT = shTm
    %             for iter = 1:nIter
%                 rphi = squeeze(circ_mean(phi_rand(:,movAvgIdx(iT),itrial,iter), [], 2));
%                 [~, rcc_sq_rand(iT,itrial,iter)] = singleFitEval(rphi, phiPred);
%             end
%         end
%     end
    
end

delete(wb);

% phi_out = phi_out - circ_mean(phi_out, [], 1); % recenter on zero after moving mean

end
%%
function [id, rcc_sq] = singleFitEval(aphi, phiPred)
% gof = squeeze(sqrt(mean(cos(phiPred-aphi)).^2 + mean(sin(phiPred-aphi)).^2)); 

gof = squeeze(abs(mean(phiPred./aphi))); % mean vector length of residuals (-> ie offset is ignored)

[~,id] = max(gof, [], 3); % get best fit ID

bestfit = nan(size(aphi));
for iTr = 1:size(id,2)
    bestfit(:,:,iTr) = squeeze(phiPred(:,1,iTr,id(:,iTr))); % get predicted phases for best fit (freq dim is redundant)
end

bestfit = nan(size(aphi));
for iFreq = 1:size(id,1)
    for iTr = 1:size(id,2)
        bestfit(:,iFreq,iTr) = squeeze(phiPred(:,1,iTr,id(iFreq,iTr))); % get predicted phases for best fit (freq dim is redundant)
    end
end



rcc_sq = sum(cSinVar(aphi) .* cSinVar(bestfit))...
    ./sqrt(...
    sum(cSinVar(aphi).^2) .* sum(cSinVar(bestfit).^2)); % circular correlation for best fit
% 
% rcc_sq = sum(sin(aphi-circ_mean(aphi)) .* sin(bestfit-circ_mean(bestfit)))...
%     ./sqrt(...
%     sum(sin(aphi-circ_mean(aphi)).^2) .* sum(sin(bestfit-circ_mean(bestfit)).^2)); % circular correlation for best fit
rcc_sq = squeeze(rcc_sq.^2);
end
%%
function out = cSinVar(in)
% computes the expression sin(angle(in) - circ_mean(angle(in))) for complex
% input (must all have unit magnitude)
out = in./mean(in);
out = imag(out./abs(out));

end