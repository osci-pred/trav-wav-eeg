function [pos, lbl, idxInData] = tw_getElecPos(labels, layout, boxLims)


lay = ft_prepare_layout(struct('layout', layout));
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
end