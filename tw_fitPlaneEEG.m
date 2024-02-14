function out = tw_fitPlaneEEG(raw, time, labels, varargin)
% raw                   eeg input; [nChan x nTime x nTrials];
% time                  time vector in seconds
% labels                cell array of electrode labels, conforming to 10/05
%                       labels (using fieldtrip layout)
%
% Optional Input (as named arg list):
%
% 'Frequency'           temporal freq limits [fmin fmax] used for bp filter
% 'AverageTrials'       bool; if true, average across trials before fitting
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
p.addParameter('AverageTrials', false);
p.addParameter('ROI', {[-0.2 0.2] [-Inf 0.25]});
p.addParameter('WindowSize', []); 
p.addParameter('MaxCycles', 1); 
p.addParameter('NumStepsSpatFreq', 30); 
p.addParameter('NumStepsWaveDir', 60); 

p.parse(varargin{:});


raw = double(raw);

boxLims = p.Results.ROI;


lay = ft_prepare_layout(struct('layout', 'EEG1005'));
box = @(xy) xy(:,1) >= boxLims{1}(1) & xy(:,1) <= boxLims{1}(2) & xy(:,2) >= boxLims{2}(1) & xy(:,2) <= boxLims{2}(2);

inclFromLay = find(box(lay.pos));
idxInData = find(cellfun(@(x) ismember(x, lay.label(inclFromLay)), labels));

idxInLay = nan(1,numel(idxInData));
for iElec = 1:numel(idxInData)
    idxInLay(iElec) = find(strcmpi(lay.label, labels{idxInData(iElec)}));
end

pos = lay.pos(idxInLay,:);
pos = pos - mean(pos); % center on zero/zero
lbl = labels(idxInData);

midLineLabel = {'Iz' 'Oz', 'POz', 'Pz', 'CPz', 'Cz', 'FCz', 'Fz'};
for i = 1:numel(midLineLabel)
    idx = find(strcmpi(lbl, midLineLabel{i}));
    if ~isempty(idx)
        midLine(i) = idx;
        midLineIdxInData(i) = idxInData(idx);
    else
        midLine(i) = nan;
        midLineIdxInData(i) = nan;
    end
end
midLine = midLine(~isnan(midLine));
midLineLabel = midLineLabel(~isnan(midLine));


nTr = size(raw,3);
nTm = size(raw,2);
sr = round(1/diff(time(1:2)));


%% Filter the data
bpFilt = designfilt('bandpassfir', 'StopbandFrequency1', 0.85*p.Results.Frequency(1), 'PassbandFrequency1', p.Results.Frequency(1),...
    'PassbandFrequency2', p.Results.Frequency(2), 'StopbandFrequency2', 1.25*p.Results.Frequency(2),...
    'StopbandAttenuation1', 60, 'PassbandRipple', 1, 'StopbandAttenuation2', 60, 'SampleRate', sr);


eegFilt = nan(numel(idxInData), nTm, nTr);

try
    for itrial = 1:nTr
        for iElec = 1:numel(idxInData)
            eegFilt(iElec,:,itrial) = filtfilt(bpFilt, raw(idxInData(iElec),:,itrial));
        end
    end
catch errmsg
    if strcmpi(errmsg.identifier, 'signal:filtfilt:InvalidDimensionsDataShortForFiltOrder')
        rawPad = cat(2,raw,raw,raw);
        eegFilt = nan(numel(idxInData), 3*nTm, nTr);
        for itrial = 1:nTr
            for iElec = 1:numel(idxInData)
                eegFilt(iElec,:,itrial) = filtfilt(bpFilt, rawPad(idxInData(iElec),:,itrial));
            end
        end
        eegFilt = eegFilt(:,nTm+1:2*nTm,:);
    else
        throw(errmsg)
    end
end

%% average across trials if evoked power is requested
if p.Results.AverageTrials
   eegFilt = mean(eegFilt,3); 
   nTr = 1;
end

%% get phase:
for itrial = 1:nTr
    h = hilbert(eegFilt(:,:,itrial)');
    phi(:,:,itrial) = angle(h)';
    pow(:,itrial) = mean(abs(h),2);
end
phi = phi - circ_mean(phi, [], 1); % subtract mean across elecs per timepoint
% -> now we can average over time (within short windows)

%% Get param space for iterative fitting:

deltapos = sqrt((pos(:,1)-pos(:,1)').^2 + (pos(:,2)-pos(:,2)').^2);
d = max(deltapos(:));

xiMax = (p.Results.MaxCycles*2*pi)/d;

xiStep = xiMax/p.Results.NumStepsSpatFreq;
thetaStep = (2*pi)/p.Results.NumStepsWaveDir;

theta = -pi:thetaStep:pi-thetaStep;
a = (xiStep:xiStep:xiMax)' .* sin(theta);
b = (xiStep:xiStep:xiMax)' .* cos(theta);

a = a(:);
b = b(:);

xi = sqrt(a.^2 + b.^2); % spatial freq (as slope)
wvdir = atan2(b,a); % wave direction 
phiPred = a.*pos(:,1)' + b.*pos(:,2)'; % predicted phases

%%
wsz = (p.Results.WindowSize/1e3) .* sr;
if ~isempty(wsz)
    movAvgIdx = @(t) (1:nTm) > t-wsz/2 & (1:nTm) < t+wsz/2;
else
    movAvgIdx = @(t) t;
end
id = nan(nTm,nTr);
pow = nan(nTm,nTr);
rcc_sq = nan(nTm,nTr);
phi_out = nan(size(pos,1), nTm, nTr);
wb = waitbar(0, {sprintf('Wave Fit: Trial %d/%d', 1,nTr), '~~~ >> ??? << ~~~'});
for itrial = 1:nTr
    waitbar(itrial/nTr,wb, {sprintf('Wave Fit: Trial %d/%d', itrial,nTr), '~~~ >> ??? << ~~~'});
    for iT = 1:nTm
        aphi = squeeze(circ_mean(phi(:,movAvgIdx(iT),itrial), [], 2));
        gof = sqrt(mean(cos(phiPred'-aphi)).^2 + mean(sin(phiPred'-aphi)).^2); % mean vector length of residuals (-> ie offset is ignored)
        [~,id(iT,itrial)] = max(gof); % get best fit ID
        bestfit = phiPred(id(iT,itrial),:)'; % get predicted phases for best fit
        rcc_sq(iT,itrial) = sum(sin(aphi-circ_mean(aphi)) .* sin(bestfit-circ_mean(bestfit)))...
            ./sqrt(...
            sum(sin(aphi-circ_mean(aphi)).^2) .* sum(sin(bestfit-circ_mean(bestfit)).^2)); % circular correlation for best fit
        rcc_sq(iT,itrial) = rcc_sq(iT,itrial)^2; 
        phi_out(:,iT,itrial) = aphi;
    end
end
delete(wb);

phi_out = phi_out - circ_mean(phi_out, [], 1); % recenter on zero after moving mean

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


