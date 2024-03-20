function [midLine, midLineLabel, midLineIdxInData] = tw_getMidLine(lbl, idxInData)
% get indices for midline electrodes for a given set of electrode labels
% passed as cell array

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

end