function [pos, lbl, idxInData, powIdx] = tw_getElecPos(labels, layout, roi, powlabels)

assert(iscell(roi), 'ROI ill-defined. Pass cell array of electrode labels or ROI box coordinates.')
lay = ft_prepare_layout(struct('layout', layout));

hasLabels = ischar(roi{1});


if hasLabels % grab electrodes as requested
    inclFromLay = cellfun(@(x) find(strcmpi(lay.label, x)), roi);
else % select all electrodes inside the box limits
    box = @(xy) xy(:,1) >= roi{1}(1) & xy(:,1) <= roi{1}(2) & xy(:,2) >= roi{2}(1) & xy(:,2) <= roi{2}(2);
    inclFromLay = find(box(lay.pos));
end


idxInData = find(cellfun(@(x) ismember(x, lay.label(inclFromLay)), labels));

idxInLay = nan(1,numel(idxInData));
for iElec = 1:numel(idxInData)
    idxInLay(iElec) = find(strcmpi(lay.label, labels{idxInData(iElec)}));
end

pos = lay.pos(idxInLay,:);
pos = pos - mean(pos); % center on zero/zero
lbl = labels(idxInData);


% get indices for additional elecs we want to return the powerspectra for:
if nargin < 4 || isempty(powlabels)
    powIdx = [];
    return
end

if ~iscell(powlabels)
   powlabels = {powlabels}; 
end

[~, powIdx] = ismember(powlabels, lbl);
powIdx = powIdx(powIdx > 0);

end