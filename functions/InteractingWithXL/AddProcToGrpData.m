function data = AddProcToGroupData(varargin)
% AddProcToGroupData - Reads paths in xl sheet for specified ablation group
%   finds missing paths from procData pointers and appends procData pointer
%   associated with the missing path
% data = AddProcToGroupData();
% data = AddProcToGroupData(pathToXLSheet);
% data = AddProToGroupData(pathToXLSheet, pathToGroupData);
% data = AddProToGroupData(pathToXLSheet, groupData);
%
% Avinash Pujala, Koyama lab/HHMI, 2016

currPath = cd;
cd('S:\Avinash\Ablations and behavior')

if nargin == 0
    pathsInXl = GetFilteredPathsFromXLS([],1);
    disp('Select mat file with group data...')
    data = OpenMatFile();
    data = data.data;
elseif nargin ==1
    if ~isdir(varargin{1})
        error('First input must be path to xls sheet with paths to procData!')
    end
    pathsInXl = GetFilteredPathsFromXLS(varargin{1},1);
    disp('Select mat file with group data...')
    data = OpenMatFile();
    data = data.data;
elseif nargin == 2
    if ~isdir(varargin{1})
        error('First input must be path to xls sheet with paths to procData!')
    end
    pathsInXl = GetFilteredPathsFromXLS(varargin{1},1);
    if isdir(varargin{2})
        disp('Select mat file with group data...')
        data = OpenMatFile();
        data = data.data;
    else
        data = varargin{2};
    end    
end

pathsInGrp = GetPathsInGrp(data);
pathsInXl = Proctify(pathsInXl);

for jj = 1:length(pathsInXl)
    if ~sum(strcmpi(pathsInGrp,pathsInXl{jj}));
        
    end
    
end

cd(currPath)

end

function pathsInXl = Proctify(pathsInXl)
for jj = 1:length(pathsInXl)
    [~,f] = fileparts(pathsInXl{jj});
    if ~strcmpi(f,'proc')
        pathsInXl{jj} = fullfile(pathsInXl{jj});
    end
end
end

function pathsInGrp = GetPathsInGrp(data)
pathsInGrp = {};
trtmnts = {'ctrl','abl'};
stims = {'vib','dark'};
count = 0;
for ii = 1:length(trtmnts)
    trtmnt = trtmnts{ii};
    for jj = 1:length(stims)
        stim = stims{jj};
        var = data.(trtmnt).(stim);
        for kk = 1:length(var.procData)
            count = count + 1;            
            [pathsInGrp{count},~] = fileparts(var.procData{kk}.Properties.Source);            
        end
    end
end
end
