function varargout = GetFilteredPathsFromXLS(varargin)
%GetFilteredPathsFromXLS Given the path to the xl sheet that has paths to
%   fish behavior data, allows for selection of paths by interactively
%   chosen filters.
% paths = GetFilteredPathsFromXLS();
% paths = GetFilteredPathsFromXLS([]);
% paths = GetFilteredPathsFromXLS(xlPath);
% [paths,paramVals] = GetFilteredPathsFromXL(xlPath);
% Inputs:
% xlPath - Path to xl sheet containing paths. If xlPath = [], then allows
% for interactive selection of file.
% Outputs:
% paths - Paths after filtering selected parameters for selected values
% paramVals = Parameter and value data structure.
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

fileDir = 'S:\Avinash\Ablations and behavior';
currDir = cd;
cd(fileDir)
if nargin==0
    [xlFile,xlPath] = uigetfile('*.xlsx');
    fullPath = fullfile(xlPath,xlFile);
elseif nargin == 1
    fullPath = varargin{1};     
end

if isempty(xlPath)
    [xlFile,xlPath] = uigetfile('*.xlsx');
    fullPath = fullfile(xlPath,xlFile);
end

[~,~,raw] = xlsread(fullPath,1);

varNames = raw(1,:);
nanInds = GetNanInds(varNames);
varNames(nanInds) = [];

disp(varNames)
filtList = input('Select filt inds, e.g. [1 3 4]: ');

rowInds_sub = 1:length(raw);
paramVals = struct;
for ff = 1:length(filtList)
    colInd = filtList(ff);
    valArray = raw(2:end,colInd);
    nanInds = find(GetNanInds(valArray));
    valArray(nanInds) = [];
    raw(nanInds+1,:)=[];
    if ischar(valArray{1})
        blah = unique(valArray);        
        fprintf('\n')
        disp(blah')
        ind = input(['Index of ' num2str(varNames{colInd}) ': ']);
        rowInds = find(strcmpi(valArray,blah(ind)));
    else
        blah = unique(cell2mat(valArray));       
        fprintf('\n')
        disp(blah')
        ind = input(['Index of ' num2str(varNames{colInd}) ': ']);
        rowInds = find(cell2mat(valArray)==blah(ind));
    end 
    if iscell(blah(ind))
        val = blah{ind};
    else
        val = blah(ind);
    end
    paramVals.(varNames{colInd}) = val;
    rowInds_sub = intersect(rowInds_sub,rowInds);
end

pathCol = strcmpi(varNames,'path');
paths = raw(rowInds_sub+1,pathCol);

varargout{1} = paths;
varargout{2} = paramVals;

end

function nanInds = GetNanInds(inputArray)
nanInds = false(length(inputArray),1);
for jj = 1:length(inputArray)
    var = inputArray{jj};
    if any(isnan(var)) || strcmpi(var,'nan')
        nanInds(jj) = true;
    end
end
end



