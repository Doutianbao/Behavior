function varargout = ReadSwimDataFromPaths(varargin)
%ReadSwimDataFromPaths Given the path to an excel sheet formatted in a
%   specific way, reads the stored paths within and extracts data. During
%   execution asks to select a filter based on ablation type.
% data = ReadSwimDataFromPaths(xlPath);
% [data,pathData]  = ReadSwimDataFromPaths(xlPath);
% [data,pathData,dataMat,dimLbls] = ReadSwimDataFromPaths(...);
% Inputs:
% xlPath - Path to xl sheet that in turn contains paths to data locations.
%   The paths must be stored in the 1st sheet.
% Outputs:
% data - Output data of the form data.x.y where x = 'abl' or 'ctrl',
%   standing for ablated or control, and y = 'vib' or 'dark', standing for
%   stimulus type, i.e. vibration or dark flash.
% pathData - Similar to data, but with paths to data instead of actual data
% dataMat = Multidimensional data matrix for ease of calculations
% dimLbls = Cell array with labels for the dimensions of dataMat.
% 
% Avinash Pujala, Koyama lab/HHMI, 2016

if nargin ==0
    [file,path] = uigetfile('*.xls');
    fPath = fullfile(path,file);
else
    fPath = varargin{1};
end

[~,~,raw] = xlsread(fPath,1);
ablationTypes = raw(2:end,1);
currentType = raw(2,1);
delInd = [];
for jj = 2:length(ablationTypes)
    if strcmpi(currentType,ablationTypes(jj)) || strcmpi(ablationTypes(jj),'nan') || sum(isnan(ablationTypes{jj}))
        delInd = [delInd; jj];
    else
        currentType = ablationTypes(jj);
    end
end
ablationTypes(delInd) =[];
[selection,ok] = listdlg('PromptString','Select ablation type...','Name','Select ablation type...',...
    'ListString',ablationTypes);
if ok
    ablationType = ablationTypes(selection);
    disp(['Ablation type is "' ablationType{1} '"'])
else
    ablationType = ablationTypes(selection);
    disp(['Defaulted to "' ablationType{1} '"'])
end

%% Filter raw cells based on ablation type
keepInds  = [];
for jj = 1:size(raw,1)
    if strcmpi(raw(jj,1),ablationType)
        keepInds = [keepInds; jj];
    end
end
raw_abFlt = raw(keepInds,:);

%% Get and store paths in data structure
disp('Reading data...')
pData = struct;
data = struct;
c1 = 0;
c2 = 0;
c3 = 0;
c4 = 0;
tic
for jj = 1:size(raw_abFlt)
    if raw_abFlt{jj,3} == 0
        if strcmpi(raw_abFlt(jj,5),'tap') || strcmpi(raw_abFlt(jj,5),'vib')
            c1 = c1 + 1;
            p = raw_abFlt{jj,6};
            [~,s] = fileparts(p);
            if strcmpi(s,'proc')==0
                p = fullfile(p,'proc');
            end
            disp(p)
            pData.ctrl.vib{c1} = p;
            blah = GetProcData(p);
            data.ctrl.vib{c1} = blah.elicitedSwimInfo;
            toc
        elseif strcmpi(raw_abFlt(jj,5),'dark')
            c2 = c2 + 1;
            p = raw_abFlt{jj,6};
            [~,s] = fileparts(p);
            if strcmpi(s,'proc')==0
                p = fullfile(p,'proc');
            end
            disp(p)
            pData.ctrl.dark{c2} = p;
            blah = GetProcData(p);
            data.ctrl.dark{c2} = blah.elicitedSwimInfo;
            toc
        end
    elseif raw_abFlt{jj,3} ==1
        if strcmpi(raw_abFlt(jj,5),'tap') || strcmpi(raw_abFlt(jj,5),'vib')
            c3 = c3 + 1;
            p = raw_abFlt{jj,6};
            [~,s] = fileparts(p);
            if strcmpi(s,'proc')==0
                p = fullfile(p,'proc');
            end
            disp(p)
            pData.abl.vib{c3} = p;
            blah = GetProcData(p);
            data.abl.vib{c3} = blah.elicitedSwimInfo;
            toc
        elseif strcmpi(raw_abFlt(jj,5),'dark')
            c4 = c4 + 1;
            p = raw_abFlt{jj,6};
            [~,s] = fileparts(p);
            if strcmpi(s,'proc')==0
                p = fullfile(p,'proc');
            end
            disp(p)
            pData.abl.dark{c4} = p;
            blah = GetProcData(p);
            data.abl.dark{c4} = blah.elicitedSwimInfo;
            toc
        end
    end
end
disp('Extracted paths and stored in pData')
disp('Extractd data  and stored in data')
toc

data.ablationType = ablationType{1};
varargout{1} = data;
varargout{2} = pData;

if nargout > 2
    [dataMat, dimLbls] = GetDataMat(data);
    varargout{3} = dataMat;
    varargout{4} = dimLbls;
end

end

function [dataMat, dimLbls] = GetDataMat(data)
%# First need to determine dimensions of multidimensional matrix
grps = {'ctrl','abl'};
dim = zeros(1,6);
for grpNum = 1:length(grps) % Dim 1  
    dim(1) = max(dim(1),grpNum);
    grp = grps{grpNum};
    blah = data.(grp);    
    disp(['Group: ' grp])
    stimTypes = fieldnames(blah);
    for stimNum = 1:length(stimTypes) % Dim 2  
        dim(2)= max(dim(2),stimNum);
        stimType = stimTypes{stimNum};
        disp(['Stim type: ' stimType])
        blah2 = blah.(stimType);
        for fishNum = 1:length(blah2) % Dim 3  
            dim(3) = max(dim(3),fishNum);
            blah3 = blah2{fishNum};
            parNames = fieldnames(blah3);
            for parNum = 1:length(parNames) % Dim 4
                dim(4) = max(dim(4),parNum);
                parName = parNames{parNum};
                blah4 = blah3.(parName);             
                for trlNum = 1:length(blah4) % Dim 5
                    dim(5) = max(dim(5),trlNum);
                    if iscell(blah4(trlNum))
                        blah5 = blah4{trlNum};                       
                    else
                        blah5 = blah4(trlNum);
                    end                    
                   if ~isempty(blah5)                        
                       for pkNum = 1:length(blah5) % Dim 6  
                           dim(6) = max(dim(6),pkNum);                                                 
                       end
                   end
                end                
            end
        end        
    end
end
disp(['Dimentions of data matrix = [ ' num2str(dim) ']']);

%# Now get data matrix
dataMat = nan(dim);
dimLbls  = cell(6,1);
dimLbls{1} = {'Ctrl','Abl'};
dimLbls{2} = {'Dark','Vib'};
dimLbls{3} = {'Fish #'};
dimLbls{4} = {'BendAmp','BendPer','Onset','BendAngVel'};
dimLbls{5} = {'Trl #'};
dimLbls{6} = {'Peak #'};

%##################################
%## ndims(dataMat) = 6 ;
%## Dim1 = # of grps - 'ctrl','abl'
%## Dim2 = # of stim types - 'dark','vib'
%## Dim3 = # of fish
%## Dim4 = # of params - 'bendAmp','bendPer','onset','bendAngVel'
%## Dim5 = # of trls
%## Dim6 = # of pks
%#################################

disp('Getting data matrix...')
for grpNum = 1:length(grps) % Dim 1
    grp = grps{grpNum};
    blah = data.(grp);
    disp(['Group: ' grp])
    stimTypes = fieldnames(blah);
    for stimNum = 1:length(stimTypes) % Dim 2
        stimType = stimTypes{stimNum};
        disp(['Stim type: ' stimType])
        blah2 = blah.(stimType);
        for fishNum = 1:length(blah2) % Dim 3
            blah3 = blah2{fishNum};
            parNames = fieldnames(blah3);
            for parNum = 1:length(parNames) % Dim 4
                parName = parNames{parNum};
                blah4 = blah3.(parName);
                for trlNum = 1:length(blah4) % Dim 5
                    if iscell(blah4(trlNum))
                        blah5 = blah4{trlNum};
                    else
                        blah5 = blah4(trlNum);
                    end
                    if ~isempty(blah5)
                        for pkNum = 1:length(blah5) % Dim 6
                            dataMat(grpNum,stimNum,fishNum,parNum,trlNum,pkNum) = blah5(pkNum);
                        end
                    end
                end
            end
        end
    end
end
disp('Done!')

end

function data = GetProcData(procDataDir)
fNames = GetFilenames(procDataDir,'.mat');
currDir = cd;
if numel(fNames) > 1
    cd(currDir)
    [fName, ~] = uigetfile('*.mat','Select the procData to read...');
    data = OpenMatFile(fullfile(procDataDir,fName));
    cd(currDir)
else
    data = OpenMatFile(fullfile(procDataDir,fNames{1}));
end
end

