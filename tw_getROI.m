function [set, idx] = tw_getROI(labels, defrows, defcols, extra)
% get line-based electrode ROI from row/column definitions.
% this is to support an alternative to the box-defined ROI (in cartesian
% space) which cuts through some columns.
% 
% labels is the full electrode set to match with. If empty, we simply return
% all combinations passed.
% defrows should be a single char or cell array of A/P label definitions ('F' 'C' 'PO' etc.)
% defcols should be an array of integers to include, 0 is used to denote 'z'
% extra is a char or cell array of chars of extra electrodes to add

assert(ischar(defrows) || (iscell(defrows) && ischar(defrows{1})),...
    'Row definition must be char or cell array of chars')

assert(isnumeric(defcols),...
    'Column definition must be integers')


if ~iscell(defrows)
   defrows = {defrows}; 
end

if nargin > 4 && ~iscell(extra)
   extra = {extra}; 
end

defcols = unique(arrayfun(@(x) {num2str(x)}, defcols));
defcols(strcmpi(defcols, '0')) = {'z'};

nRows = numel(defrows);
nCols = numel(defcols);

[rr, cc] = ndgrid(1:nRows, 1:nCols);


set = cellfun(@(x,y) [x y], defrows(rr(:)), defcols(cc(:)), 'UniformOutput', false);

if nargin > 4
    set = unique([set extra]);
end

if ~isempty(labels)
    set = set(ismember(set, labels));
    [~, idx] = ismember(set, labels);
end

end