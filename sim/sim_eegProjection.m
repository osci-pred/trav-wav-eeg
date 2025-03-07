function [eeg, t, elecLbl] = sim_eegProjection(src, t, varargin)
% project sources to head eeg model;
% input should be [nTm x nSrc x nTr] as signals at 1kHz
% t should be time vector in ms [1 x nTm]

p = inputParser;
p.addParameter('SourceAlignElec', 'Oz'); % electrodes to position sources under
p.addParameter('NumNoiseSources', 2);
p.addParameter('snrRange', [2 5]); % min/max Signal-Source to Noise-Source ratio
p.addParameter('ChansAP', {'I', 'O', 'PO', 'P', 'CP', 'C', 'FC', 'F'});
p.addParameter('ChansML', {'3' '1' 'z' '2' '4'});
p.addParameter('ElecCxDist', 20); % in mm
p.addParameter('SrcRadius', 5); % size of source in mm
p.addParameter('Plot', false);

if nargin==1 && ischar(src) && strcmpi(src, 'help')
    p.parse();
    eeg = p.Results;
    return
else
    p.parse(varargin{:});
end

[dur, nSrc, nTr] = size(src);
sr = round(1/diff(t(1:2)));

% generate noise (not correlated between trials)
for iTr = 1:nTr
    
    for iNs = 1:p.Results.NumNoiseSources
        ns(:,iNs,iTr) = sim_eegNoise(dur, sr);
    end
    
    nsAmp = rand(1,p.Results.NumNoiseSources).*diff(p.Results.snrRange) + p.Results.snrRange(1);
    ns(:,:,iTr) = (ns(:,:,iTr)./rms(ns(:,:,iTr)))./nsAmp;
    src(:,:,iTr) = src(:,:,iTr)./mean(rms(src(:,:,iTr)));
    
end

% get the electrode positions:
elec = ft_read_sens('standard_1020.elc');
elecPos = [];
elecLbl = {};
for iAP = 1:numel(p.Results.ChansAP)
    for iML = 1:numel(p.Results.ChansML)
        idx = find(strcmpi(elec.label, [p.Results.ChansAP{iAP} p.Results.ChansML{iML}]));
        if isempty(idx)
            continue
        end
        elecPos(end+1,:) = elec.chanpos(idx,:);
        elecLbl{end+1} = elec.label{idx};
    end
end


% create a virtual cortical position for every electrode by moving x mm
% towards the center of the head (at [0 0 0])

[r,a,b,c] = sphereFit(elec.chanpos);

[theta, phi] = meshgrid(0:0.05:pi, 0:0.05:2*pi);
theta = theta(:);
phi = phi(:);

sphFun = @(theta,phi,r,a,b,c) [a+r.*sin(theta).*sin(phi) b+r.*cos(theta) c+r.*sin(theta).*cos(phi)];
sphGrid = sphFun(theta,phi,r,a,b,c);

midLine = find(contains(elecLbl, 'z'));
for iElec = 1:size(elecPos,1)
    [~, grdIdx(iElec)] = min(sqrt(sum((elecPos(iElec,:) - sphGrid).^2, 2)));
end

cxPos = sphFun(theta(grdIdx), phi(grdIdx), r-p.Results.ElecCxDist,a,b,c);
% midline pos's aren't excatly zero-centered in the layout, so we fit a
% separate circle:
cxMidLine = @(x) [a+0.*sin(x) b+(r-p.Results.ElecCxDist).*sin(x) c+(r-p.Results.ElecCxDist).*cos(x)];

cxMidLineGrid = (0:0.01:2*pi)';
% bounds are Oz and Fz
[~,cxBounds(1)] =  min(sqrt(sum((cxPos(midLine(2),:) - cxMidLine(cxMidLineGrid)).^2, 2)));
[~,cxBounds(2)] =  min(sqrt(sum((cxPos(midLine(end),:) - cxMidLine(cxMidLineGrid)).^2, 2)));

cxBounds = unwrap(cxMidLineGrid(cxBounds));

srcElec = p.Results.SourceAlignElec;
    
if ~iscell(srcElec)
    srcElec = {srcElec};
end

for iSrc = 1:nSrc
    srcElecIdx(iSrc) = find(strcmpi(elecLbl, srcElec{iSrc}));
end
srcPos = cxPos(srcElecIdx,:);

%% Get source projection weights:

propFun = @(d) max(50-0.4.*d,0)./50; % this is the projection function (intensity as a function of distance)

for iSrc = 1:nSrc
    dist = sqrt(sum((elecPos - srcPos(iSrc,:)).^2, 2));
    dist = dist - p.Results.SrcRadius;
    srcFilt(:,iSrc) = propFun(dist); % spatial filter (in elec space) for this source
end
% 
%% Plot source projections to the electrodes
if p.Results.Plot
    figure
    tiledlayout(1,size(srcFilt,2))
    for iSrc = 1:nSrc
        nexttile
        scatter3(elecPos(:,1), elecPos(:,2), elecPos(:,3), 100, 'k.');
        hold on
        %         scatter3(elecPos(:,1), elecPos(:,2), elecPos(:,3), 0.5+srcFilt(:,iSrc).*4e3, srcFilt(:,iSrc).*4e3, '.');
        axis off
        gpa.surfir(elecPos(2:end,1), elecPos(2:end,2), elecPos(2:end,3), srcFilt(2:end,iSrc));
        alpha 0.8
        gpa.surfir(cxPos(2:end,1), cxPos(2:end,2), cxPos(2:end,3), 0.*srcFilt(2:end,iSrc));
        colormap(brewermap(64, '*RdBu'));
        cxLine = linspace(cxBounds(1), cxBounds(2), 1e3)';
        cxLine = cxMidLine(cxLine);
        plot3(cxLine(:,1), cxLine(:,2), cxLine(:,3), 'w', 'LineWidth', 2)
        for iElec = 1:numel(midLine)
            text(elecPos(midLine(iElec),1), elecPos(midLine(iElec),2), elecPos(midLine(iElec),3), elecLbl{midLine(iElec)})
        end
    end
end

%% Get noise projections:
% pick electrode at random for each noise source
if p.Results.NumNoiseSources > 0
    for iNs = 1:size(ns,2)
        nsPos(iNs,:) = cxPos(randi(size(cxPos,1)),:);
        dist = sqrt(sum((elecPos - nsPos(iNs,:)).^2, 2));
        dist = dist - p.Results.SrcRadius;
        nsFilt(:,iNs) =  propFun(dist);
    end
else
    nsFilt = zeros(size(srcFilt));
end

%% project sources and noise to create eeg signal:
for iTr = 1:nTr
    eeg(:,:,iTr) = ((src(:,:,iTr)*srcFilt') + (ns(:,:,iTr)*nsFilt'))';
end


% 
%% plot midline
if p.Results.Plot
    figure
    tiledlayout(2,1)
    nexttile(1)
    imagesc(t, 1:numel(midLine), eeg(midLine,:,1))
    nexttile(2)
    plot(t, eeg(midLine,:,1)')
end
%%
% keyboard
end
