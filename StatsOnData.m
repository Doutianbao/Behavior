
function summData = StatsOnData(dataDir)
summData.sideDta = sideData;


function sideData = GetLeftData(dataDir)
sideStems = {'Left','Right'};
sideData = struct;
fileNames = GetFileNames(dataDir,'.mat');
for side = 1:length(sideStems)
    fnCount = 0;
    for fn = 1:length(fileNames)
        if strfind(lower(fileNames{fn}),lower(sideStem))
            fnCount = fnCount + 1;
            fldName = [lower(sideStem) 'fnCount'];
            sideData.(fldName) = open(fullfile(dataDir,fileNames{fn}));
        end
    end
end
end


function fileNames = GetFileNames(fileDir,fileExt)
if nargin < 2
    error('At least 2 inputs required!')
end
filesInDir = dir(fileDir);
fileInDir= {filesInDir.name};
remInds = [];
for fn = 1:length(filesInDir)
    if ~strfind(filesInDir{fn},fileExt)
        remInds = [remInds; fn]
    end
end
filesInDir(remInds) = [];
fileNames = filesInDir;
end

end



