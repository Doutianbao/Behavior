function varargout = ReadSwimDataFromPaths(varargin)
%ReadSwimDataFromPaths Given the path to an excel sheet formatted in a
%   specific way, reads the stored paths within and extracts data. During
%   execution asks to select a filter based on ablation type.
% data = ReadSwimDataFromPaths(xlPath);
% [data,pathData]  = ReadSwimDataFromPaths(xlPath);
% Inputs:
% xlPath - Path to xl sheet that in turn contains paths to data locations.
%   The paths must be stored in the 1st sheet.
% Outputs:
% data - Output data of the form data.x.y where x = 'abl' or 'ctrl',
%   standing for ablated or control, and y = 'vib' or 'dark', standing for
%   stimulus type, i.e. vibration or dark flash.
% pathData - Similar to data, but with paths to data instead of actual data
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
            data.ctrl.dark{c1} = blah.elicitedSwimInfo;
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
            data.abl.dark{c3} = blah.elicitedSwimInfo;
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

varargout{1} = data;
varargout{2} = pData;



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
end

